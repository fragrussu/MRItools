# Author: Francesco Grussu, University College London
#		   <f.grussu@ucl.ac.uk>  <francegrussu@gmail.com>
#
# Code released under BSD Two-Clause license
#
# Copyright (c) 2020 University College London. 
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in 
#    the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.
# 

### Load useful modules
import argparse, os, sys
import nibabel as nib
import numpy as np
from scipy import ndimage
import warnings
import math
from scipy.special import hyp1f1        # Confluent hypergeometric function needed for channel-dependent bias correction
from scipy.special import factorial     # Factorial
from scipy.special import factorial2	 # Double factorial
import multiprocessing                 # For parallel computing
import mpdenoise as mppca              # From https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py


def correction_factor(snrval,Nch):
	# Calculate the channel-dependent and SNR-dependent correction factor
	CONSTFACTOR = math.sqrt(0.5*math.pi)  # Constant 
	snrval_squared = snrval*snrval;      # SNR squared
	confluent = hyp1f1(-0.5,float(Nch),-0.5*snrval_squared)    # confluent hypergeometric function             
	confluent_squared = confluent*confluent                    # squared confluent hypergeometric function
	Nch_minus_one = Nch - 1                                    # Number of channel reduced by one
	betaval = float(factorial2(2*Nch - 1)) / float( ( float(factorial(Nch - 1))* float(2**Nch_minus_one) ) )   # Factor beta in equation 18, Koay and Basser, JMR (2006), 179, 317-322
	factorval = CONSTFACTOR*float(betaval)                        # One more scaling required
	factorval_squared = factorval*factorval                       # Square the factor
	corrfact = 2*float(Nch) + snrval_squared - factorval_squared*confluent_squared    # Final correction factor
	return corrfact



def moment_ratio_lowerbound(N):
	# Minimum plausible value of ratio mean_signal / std_signal for a given number of coil channels
	lowestval = math.sqrt( (2.0*float(N))/correction_factor(0.0,N) - float(1.0)  )
	return lowestval



def snrmap(snr,meansig,stdsig,Ncoil):
	# Calculate SNR map g(theta) in Koay and Basser, Journal of Magnetic Resonance (2006), 179, 317-322 (equation 19) 
	mapval = math.sqrt( correction_factor(snr,Ncoil)*( 1 + ( meansig / stdsig )*( meansig / stdsig ) ) - 2*float(Ncoil) )
	return mapval



def iterate_snr_map(meansignal,stdsignal,Nchannel,tolerance):
	# Maximum number of map iteration
	NITERMAX = 80
	# Iteration counter
	niter = 1
	# Lower bound for the ratio of moments for the given number of coil channels
	minval = moment_ratio_lowerbound(Nchannel)
	# Get ratio of moments
	momratio = float(meansignal)/float(stdsignal)

	# Check whether the ratio is plausible or not
	if  momratio<=minval:
		theta = 0.0
		return theta, niter

	# Initialise the iterative map
	theta = momratio - minval	
	g_of_theta = snrmap(theta,meansignal,stdsignal,Nchannel)
        
	# Iterate map
	while niter<NITERMAX:
		tolval = math.fabs(theta - g_of_theta)
		if tolval<tolerance:
			break
		else:
			theta = g_of_theta
			g_of_theta = snrmap(theta,meansignal,stdsignal,Nchannel)
		
		niter = niter + 1

	# Return the SNR value
	return theta, niter



def MethodMoments(data):

	den_data = data[0]     # Estimate of the sample mean of the signal 
	sigma_data = data[1]   # Estimate of the sample standard deviation of the signal
	Navg = data[2]         # Number of signal averaging performed at acquisition stage before MP-PCA denoising
	myslice = data[3]      # Position of current volume in the stack of slices
	Nch = 1           # Number of coil channels (set to 1 as this leads to ideal Rician distribution of noisy MRI measurements)
	tolval = 1e-08    # Tolerance for iterative SNR estimation

	imgsize = den_data.shape   # Dimensions of 4D input
	
	# Allocate outputs
	den_data_unbias = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3])) # den_data with noise floor mitigation
	iteration_map = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))   # number of iterations
	sgmmap_corr = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))     # corrected sigma of noise
	codes_map = np.ones((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))        # output code

	for ii in range(0,imgsize[0]):
		for jj in range(0,imgsize[1]):
			for kk in range(0,imgsize[2]):
				for mm in range(0,imgsize[3]):
		
											
					# Get value of first moment for current voxel (i.e. denoised measurement)
					meanimgval = den_data[ii,jj,kk,mm]
					# Get noise standard deviation estimate
					sigmaval = np.sqrt(Navg)*sigma_data[ii,jj,kk]						
					# If sigmaval is 0 (i.e. we are in the image border), set it to a very small value
					if sigmaval==0:
						sigmaval = float(1e-12)
							
					try:
						# Find iteratively an estimate of the SNR from the ratio of signal moments
						snr_estimate, numberiter = iterate_snr_map(meanimgval,sigmaval,Nch,tolval)
						# Get the correction factor knowing SNR and number of channels
						epsilonval = correction_factor(snr_estimate,Nch)
						# Mitigate noise floor
						newval = math.sqrt( math.fabs( meanimgval*meanimgval + (1 - 2*float(Nch)/float(epsilonval))*sigmaval*sigmaval )  )
						# Correct std of noise estimate to account for noise floor
						sgmmap_corr[ii,jj,kk,mm] = np.sqrt(sigmaval*sigmaval/epsilonval) / np.sqrt(Navg);				
						# Save the unbiased and denoised measurement in a 4D matrix 
						den_data_unbias[ii,jj,kk,mm] = newval
						# Save the number of iterations for current voxel and current measurement in a 4D matrix
						iteration_map[ii,jj,kk,mm] = float(numberiter)
					except ValueError:
						sgmmap_corr[ii,jj,kk,mm] = np.nan
						iteration_map[ii,jj,kk,mm] = 0.0
						den_data_unbias[ii,jj,kk,mm] = den_data[ii,jj,kk,mm]
						codes_map[ii,jj,kk,mm] = -1


	dataout = [den_data_unbias, sgmmap_corr, iteration_map, codes_map, myslice]
	return dataout



def runMPPCA(scanin,rootout,slwindow,ricecorr=True,nsa=1,ncpu=1):
	''' MP-PCA denoising with a freely available python implementation followed by optional Rician bias mitigation  
	    

	    INTERFACES
	    runMPPCA(scanin,rootout,slwindow,ricecorr=True,nsa=1,ncpu=1)

	     
	    PARAMETERS
	    - scanin:    4D NIFTI file storing a multi-contrast MRI scan
	    - rootout:   root name for output files. The following strings will be appended:
	    			*den_den.nii (scan denoised with MP-PCA)
				*den_res.nii (residuals: denoised scan - input scan)
				*den_sigma.nii (noise standard deviation estimated by MP-PCA)
				*den_npar.nii (number of PCA components above the noise)
				*den_den_unbias.nii (denoised scan corrected for Rician bias)
				*den_sigma_unbias.nii (noise standard deviation corrected for Rician bias)
				*den_exitcode.nii (exit code for Rician bias mitigation: 1 for success, -1 for failure)
				*den_iterations.nii (number of iterations in the Rician bias mitigation)
	    - slwindow:  string indicating the sliding window kernel size in voxels (e.g. "7,7,5", "9,9,9", etc)
	    - ricecorr:  a boolean value indicating whether Rician noise floor mitigation is required
	                 (default True when function called as runMPPCA(scanin,rootout,kernel))
            - nsa:       number of signal averages performed at acquisition stage
	    - ncpu:      number of threads to be used for Rician bias mitigation with the method of moments
	                 (default 1)
	    
	    Dependencies (Python packages): mpdenoise.py, nibabel, numpy, scipy, tqdm, joblib 
	                                    (other than standard library)
	    
	    References on MP-PCA denoising: 
			
		        Veraart J, Fieremans E, Novikov DS. Diffusion MRI noise mapping
                        using random matrix theory. Magnetic Resonance in Medicine 2016, 76(5):1582-1593, 
                        doi:  10.1002/mrm.26059 

		        Veraart J , Novikov DS , Christiaens D, Ades-Aron Benjamin, Sijbers Jan, Fieremans E.
                        Denoising of diffusion MRI using random matrix theory. NeuroImage 2016, 142:394-406,
                        doi: 10.1016/j.neuroimage.2016.08.016 
                        
           Reference on noise floor mitigation (method of moments):
           
           	         Koay CG and Basser PJ, Journal of Magnetic Resonance (2006), 179, 317-322
                        doi: 10.1016/j.jmr.2006.01.016
		        
	   Freely available MP-PCA python implementation: 
           https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py

	     
	    Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''
	
	
	### Load 4D MRI scan to denoise
	print('    ... loading input data')
	try:
		scan_obj = nib.load(scanin)
	except:
		print('')
		raise RuntimeError('ERROR: the 4D input scan {} does not exist or is not in NIFTI format.'.format(scanin)) 
	scan_data = scan_obj.get_fdata()
	imgsize = scan_data.shape
	imgsize = np.array(imgsize)
	if imgsize.size!=4:
		print('')
		raise RuntimeError('ERROR: the 4D input scan {} is not a 3D NIFTI.'.format(scanin))			 
		print('')
		sys.exit(1)
	
	### Convert sliding window to a tuple
	slwindow_array = slwindow.split(',')
	slwindow_array = np.array( list(map( int, slwindow_array )) )

	### Run denoising
	print('    ... denoising with MP-PCA')
	den_data, sigma_data, npar_data = mppca.denoise( scan_data, kernel=slwindow )  # Run MP-PCA on 4D data
	res_map = den_data - scan_data   # Get residuals
	
	### Run noise floor mitigation if required
	if(ricecorr):
		print('    ... mitigating noise floor with the method of moments')

		# Smooth the 'noisy' sigma of noise map
		sigma_data_smooth = ndimage.median_filter(sigma_data, size=slwindow_array)
	
		# Allocate memory for output	
		den_data_moments = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))
		iteration_moments = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))
		sgmmap_moments = np.zeros((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))
		codes_moments = np.ones((imgsize[0],imgsize[1],imgsize[2],imgsize[3]))
	

		# Create list to run parallel processing
		inputlist = [] 
		for zz in range(0, imgsize[2]):
			buffer_den = np.zeros((imgsize[0],imgsize[1],1,imgsize[3]))
			buffer_den[:,:,0,:] = den_data[:,:,zz,:]
			buffer_sgm = np.zeros((imgsize[0],imgsize[1],1))
			buffer_sgm[:,:,0] = sigma_data_smooth[:,:,zz]
			volinfo = [buffer_den,buffer_sgm,float(nsa),zz]  # List of information relative to the zz-th MRI slice
			inputlist.append(volinfo)     # Append each slice list to create a longer list to send for parallel processing

		# Clear some memory
		del buffer_den, buffer_sgm 
	
		# Call a pool of workers whe parallel processing is required (and if the the number of slices is > 1)
		if ncpu>1 and imgsize[2]>1:

			# Make sure that the required number of threads does not exceed the maximum available
			maxthreads = multiprocessing.cpu_count()
			if(ncpu>=maxthreads):
				ncpu = maxthreads - 1

			# Create the parallel pool and give jobs to the workers
			fitpool = multiprocessing.Pool(processes=ncpu)  # Create parallel processes
			fitpool_pids_initial = [proc.pid for proc in fitpool._pool]  # Get initial process identifications (PIDs)
			fitresults = fitpool.map_async(MethodMoments,inputlist)      # Give jobs to the parallel processes
			
			# Busy-waiting: until work is done, check whether any worker dies (in that case, PIDs would change!)
			while not fitresults.ready():
				fitpool_pids_new = [proc.pid for proc in fitpool._pool]  # Get process IDs again
				if fitpool_pids_new!=fitpool_pids_initial:   # Check whether the IDs have changed from the initial values
					print('')        # Yes, they changed: at least one worker has died! Exit with error
					raise RuntimeError('ERROR: some processes died during the analysis!')			 
		
			# Work done: get results
			fitlist = fitresults.get()

			# Collect output and re-assemble MRI slices		
			for kk in range(0, imgsize[2]):					
				fitslice = fitlist[kk]    # Output relative to kk-th element in the list
				slicepos = fitslice[4]    # Spatial position of kk-th MRI slice
				den_data_moments[:,:,slicepos,:] = fitslice[0][:,:,0,:]    # Denoised data with mitigated noise floor
				sgmmap_moments[:,:,slicepos,:] = fitslice[1][:,:,0,:]      # Unbiased noise standard deviation
				iteration_moments[:,:,slicepos] = fitslice[2][:,:,0]       # Number of iterations
				codes_moments[:,:,slicepos] = fitslice[3][:,:,0]           # Exit code

		# Run serial processing as no parallel processing is required
		else:
			for kk in range(0, imgsize[2]):
				fitslice = MethodMoments(inputlist[kk])   # Fitting output relative to kk-th element in the list
				slicepos = fitslice[4]    # Spatial position of kk-th MRI slice
				den_data_moments[:,:,slicepos,:] = fitslice[0] [:,:,0,:]   # Denoised data with mitigated noise floor
				sgmmap_moments[:,:,slicepos,:] = fitslice[1][:,:,0,:]      # Unbiased noise standard deviation
				iteration_moments[:,:,slicepos] = fitslice[2][:,:,0]       # Number of iterations
				codes_moments[:,:,slicepos] = fitslice[3][:,:,0]           # Exit code
		
		
		# Get final estimate of unbiased noise standard deviation
		sgmmap_moments_mean = np.nanmean(sgmmap_moments,axis=3)
	
	
	### Store MP-PCA output
	print('    ... saving output files')
	scan_header = scan_obj.header
	scan_header.set_data_dtype('float64')   # Make sure we save output as a float64
	scan_affine = scan_obj.affine
	
	den_file = '{}_den.nii'.format(rootout)
	sigma_file = '{}_sigma.nii'.format(rootout)
	npar_file = '{}_npar.nii'.format(rootout)
	res_file = '{}_res.nii'.format(rootout)
	
	den_obj = nib.Nifti1Image(den_data,scan_affine,scan_header)
	nib.save(den_obj, den_file)
	
	sigma_obj = nib.Nifti1Image(sigma_data,scan_affine,scan_header)
	nib.save(sigma_obj, sigma_file)
	
	npar_obj = nib.Nifti1Image(npar_data,scan_affine,scan_header)
	nib.save(npar_obj, npar_file)
	
	res_obj = nib.Nifti1Image(res_map,scan_affine,scan_header)
	nib.save(res_obj, res_file)
	
	### Store noise floor mitigation output if required
	if(ricecorr):
		denunbias_file = '{}_den_unbias.nii'.format(rootout)
		sigmaunbias_file = '{}_sigma_unbias.nii'.format(rootout)
		iter_file = '{}_iterations.nii'.format(rootout)
		exit_file = '{}_exitcode.nii'.format(rootout)

		denunbias_obj = nib.Nifti1Image(den_data_moments,scan_affine,scan_header)
		nib.save(denunbias_obj, denunbias_file)
		
		sigmaunbias_obj = nib.Nifti1Image(sgmmap_moments_mean,scan_affine,scan_header)
		nib.save(sigmaunbias_obj, sigmaunbias_file)
		
		niter_obj = nib.Nifti1Image(iteration_moments,scan_affine,scan_header)
		nib.save(niter_obj, iter_file)
		
		exit_obj = nib.Nifti1Image(codes_moments,scan_affine,scan_header)
		nib.save(exit_obj, exit_file)


	### Done
	print('')




# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='** MP-PCA denoising wrapper with optional Rician bias mitigation **. This program provides a command-line interface to perform MP-PCA denoising with the freely available MP-PCA python implementation (https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py), and mitigates the noise floor bias on the MP-PCA output with a custom-implementation of the method of moments. Dependencies: mpdenoise.py, nibabel, numpy, scipy, tqdm, joblib. References on MP-PCA denoising:  Veraart J, Fieremans E, Novikov DS. Diffusion MRI noise mapping using random matrix theory. Magnetic Resonance in Medicine 2016, 76(5):1582-1593, doi:  10.1002/mrm.26059.  Veraart J , Novikov DS , Christiaens D, Ades-Aron Benjamin, Sijbers Jan, Fieremans E. Denoising of diffusion MRI using random matrix theory. NeuroImage 2016, 142:394-406, doi: 10.1016/j.neuroimage.2016.08.016. Reference on noise floor mitigation (method of moments): Koay CG and Basser PJ, Journal of Magnetic Resonance (2006), 179, 317-322 doi: 10.1016/j.jmr.2006.01.016. Author: Francesco Grussu, University College London (<f.grussu@ucl.ac.uk><francegrussu@gmail.com>). Code released under BSD Two-Clause license. Copyright (c) 2020 University College London. All rights reserved.')
	parser.add_argument('scanin', help='4D Nifti file storing the multi-contrast MRI scan to denoise')
	parser.add_argument('rootout', help='root name of output files. The following strings will be appended: *den_den.nii (scan denoised with MP-PCA); *den_res.nii (residuals: denoised scan - input scan); *den_sigma.nii (noise standard deviation estimated by MP-PCA); *den_npar.nii (number of PCA components above the noise); if Rician bias mitigation is required, also *den_den_unbias.nii (denoised scan corrected for Rician bias); *den_sigma_unbias.nii (noise standard deviation corrected for Rician bias); *den_exitcode.nii (exit code for Rician bias mitigation: 1 for success, -1 for failure); *den_iterations.nii (number of iterations in the Rician bias mitigation)')
	parser.add_argument('--kernel', metavar='<size>', default='5,5,5', help='string indicating the sliding window kernel size in voxels (e.g. "7,7,5", "9,9,9", etc; default="5,5,5"; set 1 along the slice dimension for 2D slice-wise denoising)')
	parser.add_argument('--riceunbias', metavar='<size>', default='1', help='flag controlling optional Rician noise floor mitigation with the method of moments (1: perform noise floor mitigation; 0: do not perform; default 1; anything different from 1 will be treated as 1)')
	parser.add_argument('--nsa', metavar='<N>', default='1', help='number of signal averages performed at acquisition stage, used if noise floor mitigation is required (default: --nsa equal to 1)')
	parser.add_argument('--nthread', metavar='<N>', default='1', help='number of threads to be used for noise floor mitigation, if required by --riceunbias flag (default: --nthread equal to 1)')
	args = parser.parse_args()

	### Get input arguments
	scanin = args.scanin
	rootoutstr = args.rootout
	kernelstr = args.kernel
	try:
		riceflag = np.bool( int(args.riceunbias) )
	except:
		riceflag = True
	try:
		nsigav = int(args.nsa)
	except:
		nsigav = 1
	try:
		nthr = int(args.nthread)
	except:
		nthr = 1


	print('')
	print('********************************************************************')
	print('    MP-PCA DENOISING FOLLOWED BY OPTIONAL RICIAN BIAS MITIGATION    ')
	print('********************************************************************')
	print('')
	print('Called on 4D Nifti files: {}'.format(scanin))
	print('Output root: {}'.format(rootoutstr))
	print('Sliding window kernel: {} voxels'.format(kernelstr))
	if(riceflag):
		print('Noise floor mitigation: REQUIRED - {} thread(s) will be used'.format(nthr))
		print('                                 - {} magnitude signal averages used for acquisition'.format(nsigav))	
	else:
		print('Noise floor mitigation: NOT REQUIRED')
	
	print('')

	runMPPCA(scanin, rootoutstr, kernelstr, riceflag, nsigav, nthr)

	### Done
	print('Processing completed.')
	print('')
	sys.exit(0)
	
	

