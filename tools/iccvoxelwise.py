### Calculate ICC voxel-wise
#
# Author: Francesco Grussu, UCL Queen Square Institute of Neurology
#                           UCL Centre for Medical Image Computing
#         <f.grussu@ucl.ac.uk>
#
# Copyright (c) 2019, 2020 University College London. 
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


# Import modules
import os, argparse, sys, argparse, math, shutil, math
import nibabel as nib
import numpy as np
from nipype.algorithms.icc import ICC_rep_anova as getICC

# Print help and parse arguments
parser = argparse.ArgumentParser(description='Get voxel-wise ICC and COV from a set of Nifti files storing an MRI metric and named as subj1_meas1.nii, subj2_meas1.nii, ..., subj1_meas2.nii, subj2_meas2.nii, ... . The script relies on the formalism of Shrout and Fleiss, Psychological Bulletin 1979, 86(2):420-8 "Intraclass correlations: Uses in assessing rater reliability", DOI: 10.1037/0033-2909.86.2.420. The script fits a one-way random effects model and a two-way mixed effect model voxel-by-voxel, both under the assumption of individual rating (i.e. ICC11 and ICC31 according to the formalism of Shrout and Fleiss). ICC31 is calculated using the ICC class (http://github.com/nipy/nipype/blob/master/nipype/algorithms/icc.py) in Nipype (http://nipype.readthedocs.io/en/latest), while ICC11 is calculated by estimating variances "by hand", as shown in Appendix B of Grussu F et al, NeuroImage 2015, 111:590-601 "Neurite orientation dispersion and density imaging of the healthy cervical spinal cord in vivo", DOI: 10.1016/j.neuroimage.2015.01.045.')
parser.add_argument('Folder', help='Folder containing the voxel-wise quantitative maps as Nifti (.nii).')
parser.add_argument('Nsubj', help='Number of subjects')
parser.add_argument('Nmeas', help='Number of MRI measurements per subject (note that the script does not deal with missing data, so you need exactly the same number of scans for each subject)')
parser.add_argument('Out', help='Output string for output files. Output file names for the two-may mixed effects model will append the following to the output string: *_icc31.nii (ICC values), *_bms31.nii, *_ems31.nii, *_dfc31.nii, *_dfe31.nii, *_effectF31.nii (EMS, BMS, dFC, dFE and session effect from ANOVA), *_TotVar31.nii (total variance), *_covTotVar31.nii (coefficient of variation with respect to total variance). Output file names for one-way random effects model will append the following to the output string: *_icc11.nii (ICC values), *_BVar11.nii (between-subject variance), *_WVar11.nii (within-subject variance), *_TotVar11.nii (total variance), *_covTotVar11.nii (COV with respect to total variance), *_covWVar11.nii (COV with respect to within-subject variance), *_popMean11.nii (population mean). Note that metrics for two-way mixed effects model are calculated using the ICC class in Nipype (http://github.com/nipy/nipype/blob/master/nipype/algorithms/icc.py), while metrics for the one-way random effects model are evaluated by estimating variances "by hand" as shown in Appendix B of Grussu F et al, NeuroImage 2015, 111:590-601, DOI: 10.1016/j.neuroimage.2015.01.045.')
parser.add_argument('--mask', metavar='<FILE>', help='NIFTI file containing a binary mask flagging voxels to include/exclude in the ICC analysis (1/0 in voxels to include/exclude)')
args = parser.parse_args()

# Get file names of input and output
datadir = args.Folder
Nsubj = int(args.Nsubj)
Nmeas = int(args.Nmeas)
outstring = args.Out
maskid = args.mask

# Load first measurement of first subject to get image size and reference header
refname = '{}/subj1_meas1.nii'.format(datadir)
refnifti = nib.load(refname)    # Reference NIFTI
refmat = refnifti.get_fdata()   
imgsize = refmat.shape          # Reference image size
fullmatsize = tuple([Nsubj]) + tuple([Nmeas]) + imgsize   # To create a matrix Nsubj x Nmeas x Nvox1 x Nvox2 x Novx3 storing all the measurements from all the subjects
fulldata = np.zeros(fullmatsize)
ref_header = refnifti.header     # Reference header
ref_affine = ref_header.get_best_affine()

# Loop through subjects and measurements: load data
print("")
print("****************************************************************************")
print("              ICC analysis in Python 3 with iccvoxelwise.py                 ")
print("****************************************************************************")
print("")
print("  Data folder: {}".format(datadir))
print("  Number of subjects: {}".format(Nsubj))
print("  Number of measurements per subject: {}".format(Nmeas))
print("  Output string: {}".format(outstring))
print("")
print("")
print("      Loading data for ICC calculation. Please wait...")

for ss in range(1, Nsubj+1):

	for mm in range(1, Nmeas+1):
		
		myfilename = '{}/subj{}_meas{}.nii'.format(datadir,str(ss),str(mm))          # Name of Nifti file with current map
		try:
			mynifti = nib.load(myfilename)            # Load current map
		except:
			print('')
			print('ERROR: the parametric map {} does not exist or is not in NIFTI format. Exiting with 1.'.format(myfilename))	  			 
			print('')
			sys.exit(1)
		mydata = mynifti.get_fdata()               # Extract voxel-values
		my_header = mynifti.header                # Header 
		my_affine = my_header.get_best_affine()   # Geometric information of header
		my_dims = mynifti.shape                   # Image size
		if ( (np.sum(my_affine==ref_affine)!=16) or (my_dims[0]!=imgsize[0]) or (my_dims[1]!=imgsize[1]) or (my_dims[2]!=imgsize[2]) ):
			print('')
			print('ERROR: the header of the parametric map {} does not match that of the previously loaded parametric maps. Exiting with 1.'.format(myfilename))	  			 
			print('')
			sys.exit(1)
		fulldata[ss-1,mm-1,:,:,:] = mydata        # Store values in the global data set

# Load mask if necessary
if isinstance(maskid, str)==1:
	# The user provided a mask: load it
	try:
		maskobj = nib.load(maskid)
	except:
		print('')
		print('ERROR: the mask file {} does not exist or is not in NIFTI format. Exiting with 1.'.format(maskid))	  			 
		print('')
		sys.exit(1)
	maskvals = maskobj.get_fdata()
	mask_header = maskobj.header                  # Header 
	mask_affine = mask_header.get_best_affine()   # Geometric information of header
	mask_dims = maskobj.shape                     # Image size
	if ( (np.sum(mask_affine==ref_affine)!=16) or (mask_dims[0]!=imgsize[0]) or (mask_dims[1]!=imgsize[1]) or (mask_dims[2]!=imgsize[2]) ):
		print('')
		print('ERROR: the header of the mask file {} does not match that of the parametric maps. Exiting with 1.'.format(maskid))	  			 
		print('')
		sys.exit(1)
else:
        # No mask provided: flag that you want to analyse all voxels
	maskvals = np.ones(imgsize)

# Allocate matrices for output files
icc_mat = np.zeros(imgsize)     # Matrix of ICC values for ICC(3,1) calculation
bms_mat = np.zeros(imgsize)     # Matrix of BMS values for ICC(3,1) calculation
ems_mat = np.zeros(imgsize)     # Matrix of EMS values for ICC(3,1) calculation
totvar_mat = np.zeros(imgsize)     # Matrix of e_var values for ICC(3,1) calculation
session_effect_F_mat = np.zeros(imgsize)     # Matrix of session_effect_F values for ICC(3,1) calculation
dfc_mat = np.zeros(imgsize)     # Matrix of dfc values for ICC(3,1) calculation
dfe_mat = np.zeros(imgsize)     # Matrix of dfe values for ICC(3,1) calculation
cov_mat = np.zeros(imgsize)     # Matrix of COV values obtained using total variance as estimated with the ICC(3,1) framework

icc11_mat = np.zeros(imgsize)     # Matrix of ICC values for ICC(1,1) calculation
varB11_mat = np.zeros(imgsize)    # Matrix of between-subject variance values for ICC(1,1) calculation
varW11_mat = np.zeros(imgsize)    # Matrix of within-subject variance values for ICC(1,1) calculation
varT11_mat = np.zeros(imgsize)    # Matrix of total variance values for ICC(1,1) calculation
covT11_mat = np.zeros(imgsize)    # Matrix of values obtained using total variance as estimated with the ICC(1,1) framework
covW11_mat = np.zeros(imgsize)    # Matrix of values obtained using within-subject variance as estimated with the ICC(1,1) framework
popmean11_mat = np.zeros(imgsize) # Matrix of population means as estimated with the ICC(1,1) framework

# Loop through voxels: calculate actual ICC
print("      Looping through voxels for ICC calculation. Please wait...")

for xx in range(0, imgsize[0]):
	
	for yy in range(0, imgsize[1]):
		
		for zz in range(0, imgsize[2]):
			
			# Get value of mask
			maskflag = maskvals[xx,yy,zz]
			
                        # If mask does not carry 0 for current voxel, calculate ICC; otherwise go to next voxel
			if maskflag!=0:
				
				# Calculate ICC31 using Nipype utilities (ICC(3,1) in Shrout & Fleiss, Phsy Bull (1979), 86(2): 420-428)
				data_for_anova = fulldata[:,:,xx,yy,zz]
				icc_val, r_var_val, e_var_val, session_effect_F_val, dfc_val, dfe_val = getICC(data_for_anova)
				
				# Get relevant nomenclautre
				EMS =  e_var_val                             # Residual sum of squares
				BMS = float(Nmeas)*r_var_val - EMS           # Between-target mean squares
				tot_var_val = ( BMS + (float(Nmeas) - 1)*EMS )/float(Nmeas)    # Total variance
				tot_var_val = math.fabs(tot_var_val)

				# Calculate COV as the ratio between the estimate of the population mean and the standard deviation of the total variance (ICC31)
				subjectwise_mean = np.mean(data_for_anova,axis=1)   # Array of subject-wise mean value
				cohort_mean = np.mean(subjectwise_mean)             # Cohort mean value
				tot_std_val = np.sqrt(tot_var_val)                  # Square root of total variance
				cov_val = tot_std_val/math.fabs(cohort_mean)        # COV as the ratio between std or error variance and true value

				# Store values provided by ICC/COV calculation in the location of the current voxel (ICC31)			
				icc_mat[xx,yy,zz] = icc_val
				bms_mat[xx,yy,zz] = BMS
				ems_mat[xx,yy,zz] = EMS
				session_effect_F_mat[xx,yy,zz] = session_effect_F_val
				dfc_mat[xx,yy,zz] = dfc_val
				dfe_mat[xx,yy,zz] = dfe_val
				cov_mat[xx,yy,zz] = cov_val
				totvar_mat[xx,yy,zz] = tot_var_val

				# Calculate ICC11 as Grussu et al, NeuroImage 2015, Appendix B (custom implementation of equation 1 of Shrout & Fleiss, Phsy Bull (1979), 86(2): 420-428)
				data_for_anova_mat = np.matrix(data_for_anova)
				b_of_subjects11 = np.mean(data_for_anova_mat,axis=1)        # Estimates of the value of the metric for each subject
				mu11 = np.mean(b_of_subjects11)                             # Estimate of the true population mean
				b_of_subjects_zeromean11 = b_of_subjects11 - mu11           # This is the array of the between-subject variability in the form of zero-mean
				sigmaB11_val = np.var(b_of_subjects_zeromean11)             # Estimate of the between-subject variance
				noise11 = data_for_anova_mat - np.repeat(b_of_subjects11,Nmeas,axis=1)   # Matrix of random numbers chacterising the measurement errors: it is the within-subject variability
				shape_noise11 = noise11.shape
				Nelems_11 = shape_noise11[0]*shape_noise11[1]
				Nsize_11 = tuple([1]) + tuple([Nelems_11])
				noise11_array = np.reshape(noise11,Nsize_11)                     # Matrix of random numbers reshaped as an array
				sigmaW11_val = np.var(noise11_array)                             # Estimate of the within-subject variance
				stdW11_val = np.sqrt(sigmaW11_val)                               # Standard deviation of the within-subject variance
				sigmaT11_val = sigmaW11_val + sigmaB11_val                       # Estimate of the total variance
				stdT11_val = np.sqrt(sigmaT11_val)                               # Standard deviation of the total variance
				icc11_val = sigmaB11_val/sigmaT11_val                            # Estimate of ICC11
				covT11_val = stdT11_val/math.fabs(mu11)                          # CoV from the population mean calcualated with the total variance 
				covW11_val = stdW11_val/math.fabs(mu11)                          # CoV from the population mean calcualated with the within-subject variance

				# Store values provided by ICC/COV calculation in the location of the current voxel (ICC11)
				icc11_mat[xx,yy,zz] = icc11_val
				varB11_mat[xx,yy,zz] = sigmaB11_val
				varW11_mat[xx,yy,zz] = sigmaW11_val
				varT11_mat[xx,yy,zz] = sigmaT11_val
				covT11_mat[xx,yy,zz] = covT11_val
				covW11_mat[xx,yy,zz] = covW11_val
				popmean11_mat[xx,yy,zz] = mu11

# Save the riproducibility maps as Nifti files
print("      Saving voxel-wise maps as Nifti files. Please wait...")

# Output files for ICC31
icc_file = '{}_icc31.nii'.format(outstring)
bms_file = '{}_bms31.nii'.format(outstring)
ems_file = '{}_ems31.nii'.format(outstring)
session_effect_F_file = '{}_effectF31.nii'.format(outstring)
dfc_file = '{}_dfc31.nii'.format(outstring)
dfe_file = '{}_dfe31.nii'.format(outstring)
cov_file = '{}_covTotVar31.nii'.format(outstring)
totvar_file = '{}_TotVar31.nii'.format(outstring)

# Output files for ICC11
icc11_file = '{}_icc11.nii'.format(outstring)
varB11_file = '{}_BVar11.nii'.format(outstring)
varW11_file = '{}_WVar11.nii'.format(outstring)
varT11_file = '{}_TotVar11.nii'.format(outstring)
covT11_file = '{}_covTotVar11.nii'.format(outstring)
covW11_file = '{}_covWVar11.nii'.format(outstring)
popmean11_file = '{}_popMean11.nii'.format(outstring)

# Save files for ICC31
buffer_header = refnifti.header
buffer_header.set_data_dtype('float64')   # Make sure we save quantitative maps as float64, even if input header indicates a different data type

buffernifti = nib.Nifti1Image(icc_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, icc_file)

buffernifti = nib.Nifti1Image(bms_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, bms_file)

buffernifti = nib.Nifti1Image(ems_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, ems_file)

buffernifti = nib.Nifti1Image(session_effect_F_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, session_effect_F_file)

buffernifti = nib.Nifti1Image(dfc_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, dfc_file)

buffernifti = nib.Nifti1Image(dfe_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, dfe_file)

buffernifti = nib.Nifti1Image(cov_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, cov_file)

buffernifti = nib.Nifti1Image(totvar_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, totvar_file)


# Save files for ICC11
buffernifti = nib.Nifti1Image(icc11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, icc11_file)

buffernifti = nib.Nifti1Image(varB11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, varB11_file)

buffernifti = nib.Nifti1Image(varW11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, varW11_file)

buffernifti = nib.Nifti1Image(varT11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, varT11_file)

buffernifti = nib.Nifti1Image(covT11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, covT11_file)

buffernifti = nib.Nifti1Image(covW11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, covW11_file)

buffernifti = nib.Nifti1Image(popmean11_mat,refnifti.affine,buffer_header)
nib.save(buffernifti, popmean11_file)



print("                 Done.")
print("")

sys.exit(0)

