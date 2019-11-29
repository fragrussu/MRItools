### Node: voxel-wise ADC or DKI fitting (b-value dependence only, no directional dependence)
#
# Author: Francesco Grussu, University College London
#		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
#
# Code released under BSD Two-Clause license
#
# Copyright (c) 2019 University College London. 
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.

### Load useful modules
import argparse, os, sys
import multiprocessing
import numpy as np
from scipy.optimize import minimize
import nibabel as nib


def signal_gen(mri_seq,fit_dki,tissue_par):
	''' Generate the signal for a diffusion experiment with variable TE according to ADC of DKI
		
		
	    INTERFACE
	    signal = signal_gen(mri_seq,fit_dki,tissue_par)
	    
	    PARAMETERS
	    - mri_seq: list/array of b-values (in s/mm^2)
	    - fit_dki: if 1, signal will be generated according to DKI model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC + (K/6)*(bD)^2)
	               if 0, signal will be generated according to ADC model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC) 
	    - tissue_par: list/array of tissue parameters, in the following order:
                          tissue_par[0] = S0 (T1-weighted proton density)
             		   tissue_par[1] = ADC (apparent diffusion coefficient, in ms)
             		   tissue_par[2] = K (apparent excess kurtosis coefficient, dimensionless; used only if fit_dki=1)
		
	    RETURNS
	    - signal: a numpy array of measurements generated according to the model,
			
		         S = S0*exp(-b ADC + (K/6)*(bD)^2)
		
		       where K is fixed to 0 if fit_dki = 0
		
		
	    Dependencies (Python packages): numpy
	    
	    References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group
	     
	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''
	

	### Handle inputs
	bvals = np.array(mri_seq,'float64')  # MRI sequence	
	s0_value = tissue_par[0]         # S0
	adc_value = tissue_par[1]        # ADC

	### Calculate signal
	with np.errstate(divide='raise',invalid='raise'):
		try:

			if(fit_dki==1):
				akc_value = tissue_par[2]        # Kurtosis
				signal = s0_value * np.exp((-1.0)*bvals*adc_value + (1.0/6.0)*akc_value*(bvals*adc_value)*(bvals*adc_value))
			elif(fit_dki==0):
				signal = s0_value * np.exp((-1.0)*bvals*adc_value)
			else:
				print('')
				print('ERROR: invalid flag fit_dki (set to {}, can be 0 or 1). Exiting with 1.'.format(fit_dki))					 
				print('')
				sys.exit(1)
				

		except FloatingPointError:
			signal = 0.0 * bvals      # Just output zeros when tissue parameters do not make sense			

	### Output signal
	return signal
	

def Fobj(tissue_par,fit_dki,mri_seq,meas):
	''' Fitting objective function for ADC or DKI model		
		
	    INTERFACE
	    fobj = Fobj(tissue_par,fit_dki,mri_seq,meas)
	    
	    PARAMETERS
	    - tissue_par: list/array of tissue parameters, in the following order:
                          tissue_par[0] = S0 (T1-weighted proton density)
             		   tissue_par[1] = ADC (apparent diffusion coefficient, in ms)
             		   tissue_par[2] = K (apparent excess kurtosis coefficient, dimensionless; used only if fit_dki=1)
	    - fit_dki: if 1, signal will be generated according to DKI model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC + (K/6)*(bD)^2)
	               if 0, signal will be generated according to ADC model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC)
	    - mri_seq: list/array indicating the b-values (in s/mm^2) used for the experiment (one measurement per b)
	    - meas: list/array of measurements
		
	    RETURNS
	    - fobj: objective function measured as sum of squared errors between measurements and predictions, i.e.
			
				 fobj = SUM_OVER_n( (prediction - measurement)^2 )
		
		     Above, the prediction are obtained using the signal model implemented by function signal_gen().
	    
	    References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group
	     
	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''
	
	
	### Predict signals given tissue and sequence parameters
	pred = signal_gen(mri_seq,fit_dki,tissue_par)

	### Calculate objective function and return
	fobj = np.sum( (np.array(pred) - np.array(meas))**2 )
	return fobj


def GridSearch(mri_seq,fit_dki,meas):
	''' Grid search for non-linear fitting of ADC of DKI model		
		
	    INTERFACE
	    tissue_estimate, fobj_grid = GridSearch(mri_seq,fit_dki,meas)
	    
	    PARAMETERS
	    - mri_seq: list/array indicating the b-values (in s/mm^2) used for the experiment (one measurement per b-valuee)
	    - meas: list/array of measurements
		
	    RETURNS
	    - tissue_estimate: estimate of tissue parameters that explain the measurements reasonably well. The parameters are
                          tissue_estimate[0] = S0 (T1-weighted proton density)
             		   tissue_estimate[1] = ADC (apparent diffusion coefficient, in ms)
             		   tissue_estimate[2] = K (apparent excess kurtosis coefficient, dimensionless; set to 0.0 if fit_dki=0)
	    - fit_dki: if 1, signal will be modelled according to DKI model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC + (K/6)*(bD)^2)
	               if 0, signal will be modelled according to ADC model, i.e. S = S0*exp(-TE/T2)*exp(-b ADC)
	    - fobj_grid:       value of the objective function when the tissue parameters equal tissue_estimate
	    
	    References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group
	     
	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''
	
	###	###	### CONTINUE FROM HERE
	### Prepare grid for grid search
	adc_grid = 0.001*np.logspace(np.log(0.1),np.log(32),15,base=np.exp(1.0))  # Grid of ADC values in mm^2/s
	adc_grid[0] = 0.0
	s0_grid = np.linspace(0.0,4*np.max(meas),num=15)    # Grid of S0 values
	k_grid = np.linspace(-4.0,9.0,num=15)    # Grid of kurtosis values

	### Initialise objective function to infinity and parameters for grid search
	fobj_best = float('inf')
	s0_best = 0.0
	adc_best = 0.0
	k_best = 0.0
	
	### Run grid search
	for ii in range(0, len(s0_grid)):

		s0_ii =  s0_grid[ii]
		
		for jj in range(0, len(adc_grid)):

			adc_jj =  adc_grid[jj]

			if(fit_dki==1):

				for ww in range(0, len(k_grid)):

					k_ww = k_grid[ww]

					# Tissue parameters 
					params = np.array([s0_ii,adc_jj,k_ww])
					
					# Objective function
					fval = Fobj(params,fit_dki,mri_seq,meas)

					# Check if objective function is smaller than previous value
					if fval<fobj_best:
						fobj_best = fval
						s0_best = s0_ii
						adc_best = adc_jj
						k_best = k_ww

			elif(fit_dki==0):

				# Tissue parameters 
				params = np.array([s0_ii,adc_jj])
					
				# Objective function
				fval = Fobj(params,fit_dki,mri_seq,meas)

				# Check if objective function is smaller than previous value
				if fval<fobj_best:
					fobj_best = fval
					s0_best = s0_ii
					adc_best = adc_jj
					k_best = 0.0

			else:
				print('')
				print('ERROR: invalid flag fit_dki (set to {}, can be 0 or 1). Exiting with 1.'.format(fit_dki))					 
				print('')
				sys.exit(1)
			

	### Return output
	paramsgrid = np.array([s0_best, adc_best, k_best])
	fobjgrid = fobj_best
	return paramsgrid, fobjgrid



def FitSlice(data):
	''' Fit ADC or DKI model to one MRI slice stored as a 3D numpy array (i_coord x j_coord x meas.)  
	    

	    INTERFACE
	    data_out = FitSlice(data)
	     
	    PARAMETERS
	    - data: a list of 7 elements, such that
	            data[0] is a 3D numpy array contaning the data to fit. The first and second dimensions of data[0]
		            are the slice first and second dimensions, whereas the third dimension of data[0] stores
                            measurements
		    data[1] is a numpy monodimensional array storing the b values (s/mm^2) 
		    data[2] is a string describing the fitting algorithm ("linear" or "nonlinear", see FitModel())
		    data[3] is 1 if full DKI model is to fit, 0 if only simple ADC is required 
		    data[4] is a 2D numpy array contaning the fitting mask within the MRI slice (see FitModel())
		    data[5] is a scalar containing the index of the MRI slice in the 3D volume
	    
	    RETURNS
	    - data_out: a list of 4 elements, such that
		    data_out[0] is the parameter S0 (see FitModel()) within the MRI slice
	            data_out[1] is the parameter ADC  (see FitModel()) within the MRI slice
	            data_out[2] is the parameter K  (see FitModel()) within the MRI slice
                    data_out[3] is the exit code of the fitting (see FitModel()) within the MRI slice
		     data_out[4] is the fitting sum of squared errors withint the MRI slice
                    data_out[5] equals data[5]
	
		    Fitted parameters in data_out will be stored as double-precision floating point (FLOAT64)
	    
	    References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group
	     
	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''

	
	### Extract signals and sequence information from the input list
	signal_slice = data[0]      # Signal
	seq_value = data[1]          # Sequence
	fit_algo = data[2]          # fitting algorithm
	fit_dki = data[3]          # fit for DKI? 1 --> yes; 0 --> no
	mask_slice = data[4]        # fitting mask
	idx_slice = data[5]         # Slice index
	slicesize = signal_slice.shape    # Get number of voxels of current MRI slice along each dimension
	seq_value = np.array(seq_value)     # Make sure sequence parameters are stored as an array
	
	### Check whether a sensible algorithm has been requested
	if fit_algo!="linear" and fit_algo!="nonlinear":
		print('')
		print('ERROR: unrecognised fitting algorithm. Exiting with 1.')
		print('')
		sys.exit(1)

	### Allocate output variables
	s0_slice = np.zeros(slicesize[0:2],'float64')
	adc_slice = np.zeros(slicesize[0:2],'float64')
	k_slice = np.zeros(slicesize[0:2],'float64')
	exit_slice = np.zeros(slicesize[0:2],'float64')
	mse_slice = np.zeros(slicesize[0:2],'float64')
	Nmeas = slicesize[2]   # Number of measurements


	### Fit monoexponential decay model in the voxels within the current slice
	for xx in range(0, slicesize[0]):
			for yy in range(0, slicesize[1]):
		
				# Get mask for current voxel
				mask_voxel = mask_slice[xx,yy]           # Fitting mask for current voxel

				# The voxel is not background: fit the signal model				
				if(mask_voxel==1):

					# Get signal and fitting mask
					sig_voxel = signal_slice[xx,yy,:]           # Extract signals for current voxel
					sig_voxel = np.array(sig_voxel)           # Convert to array
						
					## Simplest case: there are only two echo times --> get the solution analytically
					if(Nmeas==2):

						if(fit_dki==1):
							print('')
							print('ERROR: you cannot fit for a full DKI model as you have only two measurements. Exiting with 1.')					 
							print('')
							sys.exit(1)
							
						sig1 = sig_voxel[0]                            # Signal for first b-values
						sig2 = sig_voxel[1] 		                # Signal for second b-value
						b1 = seq_value[0]                              # First b
						b2 = seq_value[1]                              # Second b
						
						# Calculate maps analytically, handling warnings
						with np.errstate(divide='raise',invalid='raise'):	
							try:
								adc_voxel =  np.log( sig1/sig2 ) / (b2 - b1)
								k_voxel = 0.0
								s0_voxel = sig1 / np.exp( (-1.0)*b1*adc_voxel )
								exit_voxel = 1

								# Check whether the solution is plausible
								if adc_voxel<0:
									s0_voxel = np.mean(sig_voxel)
									adc_voxel = 0.001*32.0    # Maximum possible ADC in mm2/s
									exit_voxel = -1
								if s0_voxel<0:
									s0_voxel = 0.0
									exit_voxel = -1

								mse_voxel = Fobj([s0_voxel,adc_voxel,k_voxel],fit_dki,seq_value,sig_voxel)   # Error (0 if  adc > 0 ad s0 > 0 at the first attempt) 
								
								
							except FloatingPointError:
								s0_voxel = 0.0
								adc_voxel = 0.0
								k_voxel = 0.0
								exit_voxel = -1
								mse_voxel = 0.0

					## General case: there are more than two echo times --> get the solution minimising an objective function
					else:

						# Perform linear fitting as first thing - if non-linear fitting is required, the linear fitting will be used to initialise the non-linear optimisation afterwards
						seq_column = np.reshape(seq_value,(Nmeas,1))    # Store TE values as a column array						
						sig_voxel_column = np.reshape(sig_voxel,(Nmeas,1))   # Reshape measurements as column array

						# Calculate linear regression coefficients as ( W * Q )^-1 * (W * m), while handling warnings
						with np.errstate(divide='raise',invalid='raise'):
							try:
								# Create matrices and arrays to be combinted via matrix multiplication
								Yvals = np.log(sig_voxel)           # Independent variable of linearised model
								Xvals = (-1.0)*seq_column            # Dependent variable of linearised model
								Xvals2 = seq_column*seq_column            # Dependent variable of linearised model
								allones = np.ones([Nmeas,1])        # Column of ones						
								if(fit_dki==1):
									Qmat = np.concatenate((allones,Xvals,Xvals2),axis=1)    # Design matrix Q
								else:
									Qmat = np.concatenate((allones,Xvals),axis=1)    # Design matrix Q
								Wmat = np.diag(sig_voxel)                        # Matrix of weights W
									
								# Calculate coefficients via matrix multiplication
								coeffs = np.matmul( np.linalg.pinv( np.matmul(Wmat,Qmat) ) , np.matmul(Wmat,Yvals) )
									
								# Retrieve signal model parameters from linear regression coefficients
								s0_voxel = np.exp(coeffs[0])
								adc_voxel = coeffs[1]
								if(fit_dki==1):
									k_voxel = 6.0*coeffs[2] / (adc_voxel*adc_voxel)
								else:
									k_voxel = 0.0
								exit_voxel = 1	
								
								# Check whether the solution is plausible: if not, declare fitting failed
								if adc_voxel<0:
									s0_voxel = np.mean(sig_voxel)
									adc_voxel = 0.001*32.0    # We fix the maximum possible value
									k_voxel = 0.0
									exit_voxel = -1
								if s0_voxel<0:
									s0_voxel = 0.0
									exit_voxel = -1	

								mse_voxel = Fobj([s0_voxel,adc_voxel,k_voxel],fit_dki,seq_value,sig_voxel)   # Measure of quality of fit
								
								
							except FloatingPointError:
								s0_voxel = 0.0
								adc_voxel = 0.0
								k_voxel = 0.0
								exit_voxel = -1
								mse_voxel = 0.0
								
						
						# Refine the results from linear with non-linear optimisation if the selected algorithm is "nonlinear"
						if fit_algo=="nonlinear":

							# Check whether linear fitting has failed
							if exit_voxel==-1:
								param_init, fobj_init = GridSearch(seq_value,fit_dki,sig_voxel)   # Linear fitting has failed: run a grid search
								if(fit_dki==0):
									param_init = [param_init[0],param_init[1]] 

							else:
								if(fit_dki==1):
									param_init = [s0_voxel,adc_voxel,k_voxel]   # Linear fitting 
								else:
									param_init = [s0_voxel,adc_voxel]   # Linear fitting

								fobj_init = mse_voxel               
							
							# Minimise the objective function numerically
							if(fit_dki==1):
								param_bound = ((0,2*s0_voxel),(0,0.001*32),(-4.0,9.0),)  # Range for parameters						
							else:
								param_bound = ((0,2*s0_voxel),(0,0.001*32),)  # Range for parameters						
							modelfit = minimize(Fobj, param_init, method='L-BFGS-B', args=tuple([fit_dki,seq_value,sig_voxel]), bounds=param_bound)
							fit_exit = modelfit.success
							fobj_fit = modelfit.fun

							# Get fitting output if non-linear optimisation was successful and if succeeded in providing a smaller value of the objective function as compared to the grid search
							if fit_exit==True and fobj_fit<fobj_init:
								param_fit = modelfit.x
								s0_voxel = param_fit[0]
								adc_voxel = param_fit[1]
								if(fit_dki==1):
									k_voxel = param_fit[2]
								else:
									k_voxel = 0.0
								exit_voxel = 1
								mse_voxel = fobj_fit

							# Otherwise, output the best we could find with linear fitting or, when linear fitting fails, with grid search (note that grid search cannot fail by implementation)
							else:
								s0_voxel = param_init[0]
								adc_voxel = param_init[1]
								if(fit_dki==1):
									k_voxel = param_init[2]
								else:
									k_voxel = 0.0
								exit_voxel = -1
								mse_voxel = fobj_init
							
						
							
				# The voxel is background
				else:
					s0_voxel = 0.0
					adc_voxel = 0.0
					k_voxel = 0.0
					exit_voxel = 0
					mse_voxel = 0.0
				
				# Store fitting results for current voxel
				s0_slice[xx,yy] = s0_voxel
				adc_slice[xx,yy] = adc_voxel
				k_slice[xx,yy] = k_voxel
				exit_slice[xx,yy] = exit_voxel
				mse_slice[xx,yy] = mse_voxel

	### Create output list storing the fitted parameters and then return
	data_out = [s0_slice, adc_slice, k_slice, exit_slice, mse_slice, idx_slice]
	return data_out
	



def FitModel(*argv):
	''' Fit ADC or DKI to multi b-value data 
	    

	    INTERFACES
	    FitModel(mb_nifti, b_text, fit_dki, output_basename, algo, ncpu)
	    FitModel(mb_nifti, b_text, fit_dki, output_basename, algo, ncpu, mask_nifti)
	    FitModel(mb_nifti, b_text, fit_dki, output_basename, algo, ncpu, te_file, t2_nifti)
	    FitModel(mb_nifti, b_text, fit_dki, output_basename, algo, ncpu, mask_nifti, te_file, t2_nifti)
	     
	    PARAMETERS
	    - mb_nifti: path of a Nifti file storing the multi b-value data as 4D data.
	    - b_text: path of a text file storing the b-values (s/mm^2) used to acquire the data.
	    - fit_dki: flag to specify whether to fit DKI (if set to 1) or simple ADC (if set to 0)
	    - output_basename: base name of output files. Output files will end in 
                            "_S0.nii"   --> T1-weighted proton density, with receiver coil field bias
		            "_ADC.nii"  --> ADC map (mm^2/2)
		            "_K.nii"  --> Kurtosis excess map (dimensionless, if requested)
			    "_Exit.nii" --> exit code (1: successful fitting; 0 background; -1: unsuccessful fitting)
			    "_SSE.nii"  --> fitting sum of squared errors
			    
			    Note that in the background and where fitting fails, metrics are set to 0.0
			    Output files will be stored as double-precision floating point (FLOAT64)
			
	    - algo: fitting algorithm ("linear" or "nonlinear")
	    - ncpu: number of processors to be used for computation
	    - mask_nifti: path of a Nifti file storing a binary mask, where 1 flgas voxels where the 
			  signal model needs to be fitted, and 0 otherwise
	    - te_file: text file storing TEs if TE varies across measurements
	    - t2_nifti: NIFTI file storting a quantitative T2 map to account for measurement-dependent TE
	    
	    References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group
	     
	    Dependencies: numpy, nibabel, scipy (other than standard library)

	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''
	

	### Get input parametrs
	Nargv = len(argv)
	sig_nifti = argv[0]
	seq_text = argv[1]
	fit_dki = argv[2]
	output_rootname = argv[3]
	algo = argv[4]
	ncpu = argv[5]
	ncpu_physical = multiprocessing.cpu_count()
	if ncpu>ncpu_physical:
		print('')
		print('WARNING: {} CPUs were requested. Using {} instead (all available CPUs)...'.format(ncpu,ncpu_physical))					 
		print('')
		ncpu = ncpu_physical     # Do not open more workers than the physical number of CPUs

	### Check whether the requested fitting algorithm makes sense or not
	if algo!="linear" and algo!="nonlinear":
		print('')
		print('ERROR: unrecognised fitting algorithm. Exiting with 1.')
		print('')
		sys.exit(1)

	### Load MRI data
	print('    ... loading input data')
	
	# Make sure MRI data exists
	try:
		sig_obj = nib.load(sig_nifti)
	except:
		print('')
		print('ERROR: the 4D input NIFTI file {} does not exist or is not in NIFTI format. Exiting with 1.'.format(me_nifti))					 
		print('')
		sys.exit(1)
	
	# Get image dimensions and convert to float64
	sig_data = sig_obj.get_data()
	imgsize = sig_data.shape
	sig_data = np.array(sig_data,'float64')
	imgsize = np.array(imgsize)
	
	# Make sure that the text file with sequence parameters exists and makes sense
	try:
		seqarray = np.loadtxt(seq_text)
		seqarray = np.array(seqarray,'float64')
		seqarray_size = seqarray.size
	except:
		print('')
		print('ERROR: the b-value file {} does not exist or is not a numeric text file. Exiting with 1.'.format(seq_text))					 
		print('')
		sys.exit(1)
			
	# Check consistency of sequence parameter file and number of measurements
	if imgsize.size!=4:
		print('')
		print('ERROR: the input file {} is not a 4D nifti. Exiting with 1.'.format(sig_nifti))					 
		print('')
		sys.exit(1)
	if seqarray_size!=imgsize[3]:
		print('')
		print('ERROR: the number of measurements in {} does not match the number of b-values in {}. Exiting with 1.'.format(sig_nifti,seq_text))					 
		print('')
		sys.exit(1)
	seq = seqarray
	mask_data = np.ones(imgsize[0:3],'float64')

	### Deal with optional arguments: mask
	if (Nargv==7) or (Nargv==9):
		mask_nifti = argv[6]
		try:
			mask_obj = nib.load(mask_nifti)
		except:
			print('')
			print('ERROR: the mask file {} does not exist or is not in NIFTI format. Exiting with 1.'.format(mask_nifti))					 
			print('')
			sys.exit(1)

		
		# Make sure that the mask has consistent header with the input data containing the measurements
		sig_header = sig_obj.header
		sig_affine = sig_header.get_best_affine()
		sig_dims = sig_obj.shape
		mask_dims = mask_obj.shape		
		mask_header = mask_obj.header
		mask_affine = mask_header.get_best_affine()	
		
		# Make sure the mask is a 3D file
		mask_data = mask_obj.get_data()
		masksize = mask_data.shape
		masksize = np.array(masksize)
		if masksize.size!=3:
			print('')
			print('WARNING: the mask file {} is not a 3D Nifti file. Ignoring mask...'.format(mask_nifti))				 
			print('')
			mask_data = np.ones(imgsize[0:3],'float64')
		elif ( (np.sum(sig_affine==mask_affine)!=16) or (sig_dims[0]!=mask_dims[0]) or (sig_dims[1]!=mask_dims[1]) or (sig_dims[2]!=mask_dims[2]) ):
			print('')
			print('WARNING: the geometry of the mask file {} does not match that of the input data. Ignoring mask...'.format(mask_nifti))					 
			print('')
			mask_data = np.ones(imgsize[0:3],'float64')
		else:
			mask_data = np.array(mask_data,'float64')
			# Make sure mask data is a numpy array
			mask_data[mask_data>0] = 1
			mask_data[mask_data<=0] = 0


	### Deal with optional arguments: TE and T2 files
	t2_data = np.ones(imgsize[0:3],'float64')
	tearray = 0.0*seqarray
	if Nargv==8:

		try:
			tearray = np.loadtxt(argv[6])
			tearray = np.array(tearray,'float64')
			tearray_size = tearray.size
		except:
			print('')
			print('ERROR: the TE file {} does not exist or is not a numeric text file. Exiting with 1.'.format(argv[6]))					 
			print('')
			sys.exit(1)

		if tearray_size!=imgsize[3]:
			print('')
			print('ERROR: the number of measurements in {} does not match the number of echo times in {}. Exiting with 1.'.format(sig_nifti,argv[6]))					 
			print('')
			sys.exit(1)


		try:
			t2_obj = nib.load(argv[7])
			t2_niftifile = argv[7]
		except:
			print('')
			print('ERROR: the T2 map file {} does not exist or is not in NIFTI format. Exiting with 1.'.format(argv[7]))					 
			print('')
			sys.exit(1)


	if Nargv==9:

		try:
			tearray = np.loadtxt(argv[7])
			tearray = np.array(tearray,'float64')
			tearray_size = tearray.size
		except:
			print('')
			print('ERROR: the TE file {} does not exist or is not a numeric text file. Exiting with 1.'.format(argv[7]))					 
			print('')
			sys.exit(1)

		if tearray_size!=imgsize[3]:
			print('')
			print('ERROR: the number of measurements in {} does not match the number of echo times in {}. Exiting with 1.'.format(sig_nifti,argv[7]))					 
			print('')
			sys.exit(1)


		try:
			t2_obj = nib.load(argv[8])
			t2_niftifile = argv[8]
		except:
			print('')
			print('ERROR: the T2 map file {} does not exist or is not in NIFTI format. Exiting with 1.'.format(argv[8]))					 
			print('')
			sys.exit(1)


	if ((Nargv==8) or (Nargv==9)):
	
		# Make sure that the T2 map has consistent header with the input data containing the measurements
		sig_header = sig_obj.header
		sig_affine = sig_header.get_best_affine()
		sig_dims = sig_obj.shape
		t2_dims = t2_obj.shape		
		t2_header = t2_obj.header
		t2_affine = t2_header.get_best_affine()	
		
		# Make sure the mask is a 3D file
		t2_data = t2_obj.get_data()
		t2size = t2_data.shape
		t2size = np.array(t2size)
		if t2size.size!=3:
			print('')
			print('WARNING: the T2 map file {} is not a 3D Nifti file. Ignoring mask...'.format(t2_niftifile))				 
			print('')
			t2_data = np.ones(imgsize[0:3],'float64')
		elif ( (np.sum(sig_affine==t2_affine)!=16) or (sig_dims[0]!=t2_dims[0]) or (sig_dims[1]!=t2_dims[1]) or (sig_dims[2]!=t2_dims[2]) ):
			print('')
			print('WARNING: the geometry of the T2 map file {} does not match that of the input data. Ignoring mask...'.format(t2_niftifile))					 
			print('')
			t2_data = np.ones(imgsize[0:3],'float64')
		else:
			t2_data = np.array(t2_data,'float64')
			# Make sure T2 data is a numpy array
			t2_data[t2_data<=0] = 3000.0  # in ms

		# Correct for measaurement-dependent TE
		for uu in range(0, imgsize[0]):
			for vv in range(0, imgsize[1]):
				for ww in range(0, imgsize[2]):
					sig_copy = np.copy(sig_data[uu,vv,ww,:])
					sig_data[uu,vv,ww,:] = sig_copy / np.exp((-1.0)*tearray/t2_data[uu,vv,ww])

		sig_data[np.isnan(sig_data)]=0.0
		sig_data[np.isinf(sig_data)]=0.0 


	### Allocate memory for outputs
	s0_data = np.zeros(imgsize[0:3],'float64')	       # T1-weighted proton density with receiver field bias (double-precision floating point)
	adc_data = np.zeros(imgsize[0:3],'float64')	       # ADC (double-precision floating point)
	k_data = np.zeros(imgsize[0:3],'float64')	       # K (double-precision floating point)
	exit_data = np.zeros(imgsize[0:3],'float64')           # Exit code (double-precision floating point)
	mse_data = np.zeros(imgsize[0:3],'float64')            # Fitting sum of squared errors (MSE) (double-precision floating point)


	#### Fitting
	print('    ... diffusion analysis')
	# Create the list of input data
	inputlist = [] 
	for zz in range(0, imgsize[2]):
		sliceinfo = [sig_data[:,:,zz,:],seq,algo,fit_dki,mask_data[:,:,zz],zz]  # List of information relative to the zz-th MRI slice
		inputlist.append(sliceinfo)     # Append each slice list and create a longer list of MRI slices whose processing will run in parallel

	# Clear some memory
	del sig_data, mask_data 
	
	# Call a pool of workers to run the fitting in parallel if parallel processing is required (and if the the number of slices is > 1)
	if ncpu>1 and imgsize[2]>1:

		# Create the parallel pool and give jobs to the workers
		fitpool = multiprocessing.Pool(processes=ncpu)  # Create parallel processes
		fitpool_pids_initial = [proc.pid for proc in fitpool._pool]  # Get initial process identifications (PIDs)
		fitresults = fitpool.map_async(FitSlice,inputlist)      # Give jobs to the parallel processes
		
		# Busy-waiting: until work is done, check whether any worker dies (in that case, PIDs would change!)
		while not fitresults.ready():
			fitpool_pids_new = [proc.pid for proc in fitpool._pool]  # Get process IDs again
			if fitpool_pids_new!=fitpool_pids_initial:               # Check whether the IDs have changed from the initial values
				print('')					 # Yes, they changed: at least one worker has died! Exit with error
				print('ERROR: some processes died during parallel fitting. Exiting with 1.')					 
				print('')
				sys.exit(1)
		
		# Work done: get results
		fitlist = fitresults.get()

		# Collect fitting output and re-assemble MRI slices		
		for kk in range(0, imgsize[2]):					
			fitslice = fitlist[kk]    # Fitting output relative to kk-th element in the list
			slicepos = fitslice[5]    # Spatial position of kk-th MRI slice
			s0_data[:,:,slicepos] = fitslice[0]    # Parameter S0 
			adc_data[:,:,slicepos] = fitslice[1]   # Parameter ADC
			k_data[:,:,slicepos] = fitslice[2]     # Parameter Kurtosis
			exit_data[:,:,slicepos] = fitslice[3]  # Exit code
			mse_data[:,:,slicepos] = fitslice[4]   # Sum of Squared Errors	


	# Run serial fitting as no parallel processing is required (it can take up to 1 hour per brain)
	else:
		for kk in range(0, imgsize[2]):
			fitslice = FitSlice(inputlist[kk])   # Fitting output relative to kk-th element in the list
			slicepos = fitslice[5]    # Spatial position of kk-th MRI slice
			s0_data[:,:,slicepos] = fitslice[0]    # Parameter S0 
			adc_data[:,:,slicepos] = fitslice[1]   # Parameter ADC
			k_data[:,:,slicepos] = fitslice[2]     # Parameter Kurtosis
			exit_data[:,:,slicepos] = fitslice[3]  # Exit code
			mse_data[:,:,slicepos] = fitslice[4]   # Sum of Squared Errors


	### Save the output maps
	print('    ... saving output files')

	if(fit_dki==1):
		buffer_string=''
		seq_string = (output_rootname,'_K.nii')
		k_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (output_rootname,'_ADC.nii')
	adc_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (output_rootname,'_S0.nii')
	s0_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (output_rootname,'_Exit.nii')
	exit_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (output_rootname,'_SSE.nii')
	mse_outfile = buffer_string.join(seq_string)

	buffer_header = sig_obj.header
	buffer_header.set_data_dtype('float64')   # Make sure we save quantitative maps as float64, even if input header indicates a different data type

	if(fit_dki==1):
		k_obj = nib.Nifti1Image(k_data,sig_obj.affine,buffer_header)
		nib.save(k_obj, k_outfile)

	adc_obj = nib.Nifti1Image(1000.0*adc_data,sig_obj.affine,buffer_header)  # Convert ADC in um^2/ms from mm2/s when saving
	nib.save(adc_obj, adc_outfile)

	s0_obj = nib.Nifti1Image(s0_data,sig_obj.affine,buffer_header)
	nib.save(s0_obj, s0_outfile)

	exit_obj = nib.Nifti1Image(exit_data,sig_obj.affine,buffer_header)
	nib.save(exit_obj, exit_outfile)

	mse_obj = nib.Nifti1Image(mse_data,sig_obj.affine,buffer_header)
	nib.save(mse_obj, mse_outfile)

	### Done
	print('')




# Run the module as a script when required
if __name__ == "__main__":

	
	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Voxel-wise fitting of isotropic (monodimensional) ADC or DKI model from multi-b-value MRI magnitude data already corrected for motion. Dependencies (Python packages): numpy, nibabel, scipy (other than standard library). References: "Quantitative MRI of the brain", 2nd edition, Tofts, Cercignani and Dowell editors, Taylor and Francis Group. Author: Francesco Grussu, University College London. Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>.')
	parser.add_argument('dwi_file', help='4D Nifti file of multi-b-value magnitude images from a diffusion MRI experiment')
	parser.add_argument('b_file', help='text file of b-values used to acquire the images (in s/mm^2; b-values separated by spaces)')
	parser.add_argument('out_root', help='root of output file names, to which file-specific strings will be added; output files will be double-precision floating point (FLOAT64) and will end in "_S0.nii" (T1-weighted proton density, with receiver coil field bias); "_ADC.nii" (ADC map in um^2/ms); "_K.nii" (dimensionless excess kurtosis map, if requested with flag --fitdki); "_Exit.nii" (exit code: 1 for successful non-linear fitting; 0 background; -1 for failing of non-linear fitting, with results from a grid search/linear fitting provided instead); "_SSE.nii" (fitting sum of squared errors).')
	parser.add_argument('--mask', metavar='<file>', help='mask in Nifti format where 1 flags voxels where fitting is required, 0 where is not')
	parser.add_argument('--algo', metavar='<type>', default='linear', help='fitting algorithm; choose among "linear" and "nonlinear" (default: "linear")')
	parser.add_argument('--ncpu', metavar='<N>', help='number of CPUs to be used for computation (default: half of available CPUs)')
	parser.add_argument('--fitdki', metavar='<flag>', default='0', help='flag indicating whether to fit ADC (if set to 0) or DKI (if set to 1; defalut 0 --> fit ADC). Any value different from 0 or 1 will be treated as 0.')
	parser.add_argument('--tefile', metavar='<file>', help='text file of echo times TEs if different echo times have been used for different TEs (TEs in ms, separated by a space)')
	parser.add_argument('--t2map', metavar='<file>', help='path of a 3D NIFTI file storing a T2 map (expressed in ms) to account for measurement-dependent TE (compulsory if a TE file is specified with option --tefile)')
	args = parser.parse_args()

	### Get input arguments
	sigfile = args.dwi_file
	seqfile = args.b_file
	outroot = args.out_root
	maskfile = args.mask
	fittype = args.algo
	nprocess = args.ncpu
	tefile = args.tefile
	t2map = args.t2map
	fitdki = int(args.fitdki)

	### Deal with optional arguments
	if isinstance(maskfile, str)==1:
		# A mask for fitting has been provided
		maskrequest = True
	else:
		# A mask for fitting has not been provided
		maskrequest = False

	
	if isinstance(nprocess, str)==1:
		# A precise number of CPUs has been requested
		nprocess = int(float(nprocess))
	else:
		# No precise number of CPUs has been requested: use 50% of available CPUs
		nprocess = multiprocessing.cpu_count()
		nprocess = int(float(nprocess)/2)


	if isinstance(tefile, str)==1:
		t2request = True
		if isinstance(t2map, str)==0:
			print('')
			print('ERROR! You have specified a TE text file but not a T2 map! Exiting with 1.')
			print('')
			sys.exit(1)
	else:
		t2request = False

	if isinstance(t2map, str)==1:
		if isinstance(tefile, str)==0:
			print('')
			print('ERROR! You have specified a T2 map but not a TE text file! Exiting with 1.')
			print('')
			sys.exit(1)


	### Sort out things to print
	buffer_str=''
	seq_str = (outroot,'_ADC.nii')
	adc_out = buffer_str.join(seq_str)

	buffer_str=''
	seq_str = (outroot,'_K.nii')
	k_out = buffer_str.join(seq_str)
		
	buffer_str=''
	seq_str = (outroot,'_S0.nii')
	s0_out = buffer_str.join(seq_str)

	buffer_str=''
	seq_str = (outroot,'_Exit.nii')
	exit_out = buffer_str.join(seq_str)

	buffer_str=''
	seq_str = (outroot,'_SSE.nii')
	mse_out = buffer_str.join(seq_str)


	print('')
	print('********************************************************************************')
	if(fitdki==1):
		print('                 Fitting of monodimensional DKI                 ')
	else:
		print('                 Fitting of monodimensional ADC                 ')
	print('********************************************************************************')
	print('')
	print('Called on 4D Nifti file: {}'.format(sigfile))
	print('b-value file: {}'.format(seqfile))
	if(fitdki==1):
		print('Output files: {}, {}, {}, {}'.format(adc_out,k_out,s0_out,exit_out,mse_out))
	else:
		print('Output files: {}, {}, {}'.format(adc_out,s0_out,exit_out,mse_out))
	print('')

	### Call fitting routine
	# The entry point of the parallel pool has to be protected with if(__name__=='__main__') (for Windows): 
	if(__name__=='__main__'):
		if (maskrequest==False):
			if(t2request==False):
				FitModel(sigfile, seqfile, fitdki, outroot, fittype, nprocess)
			else:
				FitModel(sigfile, seqfile, fitdki, outroot, fittype, nprocess, tefile, t2map)
		else:
			if(t2request==False):
				FitModel(sigfile, seqfile, fitdki, outroot, fittype, nprocess, maskfile)
			else:
				FitModel(sigfile, seqfile, fitdki, outroot, fittype, nprocess, maskfile, tefile, t2map)
	
	### Done
	print('Processing completed.')
	print('')
	sys.exit(0)



