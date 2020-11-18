# Code released under BSD Two-Clause license
#
# Copyright (c) 2020 University College London. 
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.
import nibabel as nib
import numpy as np
import sys, argparse
from scipy import ndimage
from scipy import signal


def rephaseDC(img_real,img_imag,out_base,kernel):
	''' Rephase complex MR data in image space
	
	INTERFACE
	rephaseDC(img_real,img_imag,out_root,kernel)

	PARAMETERS
	img_real: path of a NIFTI file storing a 3D or 4D image (real channel)
	img_imag: path of a NIFTI file storing a 3D or 4D image (imaginary channel)
	out_base: base name of output files (the output files contain the real and imaginary channels 
		  after noise decorrelation and rephasing; these will end in *_RealReph.nii 
		  (real channel rephased), *_RealRephThresh.nii (real channel rephased with outlier 
		  detection), *_ImagReph.nii (imaginary channel rephased), *_ImagRephThresh.nii 
		  (imaginary channel rephased with outlier detection), *_OutlierDetected.nii 
		  (flagging with 1 outliers), *_PhaseOriginal.nii (storing the original phase),
		  *_PhaseBackground.nii (storing the estimated background phase), *_PhaseRephased.nii 
		  (storing the phase after rephasing using the background phase), *_PhaseRephasedOutliers.nii 
	          (storing the original phase after rephasing where outliers are set to zero phase).
	          Note that the imaginary channel after rephasing should contain mostly noise and 
		  negligible true signal information. 
	kernel:   string of the 2D kernel to use for decorrenation (choose among "B3", "B5", "G3F1", "G5F2",
	          "G3F1H", "G5F2H", "Opt3", "Opt5"; see Sprenger T et al, MRM 2017, 77:559–570 for more 
		  information about the kernels.)

	DESCRIPTION
	The function implements noise decorrelation and rephasing algorithm presented in Sprenger T et al, 
	MRM 2017, 77:559-570. The function works with 3D and 4D NIFTI files (in the latter case,
	each volume of the 4D NIFTI is treated independently).

	References: "Real valued diffusion-weighted imaging using decorrelated 
		     phase filtering", Sprenger T et al, Magnetic Resonance 
		     in Medicine (2017), 77:559-570
	     
	Author: Francesco Grussu, University College London
		<f.grussu@ucl.ac.uk> <francegrussu@gmail.com>

	Code released under BSD Two-Clause license. 
	Copyright (c) 2020 University College London. All rights reserved.'''


	# Load real MRI
	try:
		imgR_obj = nib.load(img_real)
	except:
		print('')
		print('ERROR: the file storing the real channel {} does not exist or is not in NIFTI format. Exiting with 1.'.format(img_real))					 
		print('')
		sys.exit(1)
	
	imgR_data = imgR_obj.get_fdata()
	imgR_size = imgR_data.shape
	imgR_size = np.array(imgR_size)
	imgR_ndim = imgR_size.size
	imgR_data = np.array(imgR_data,'float64')
	imgR_header = imgR_obj.header
	imgR_affine = imgR_header.get_best_affine()	

	# Load imaginary MRI
	try:
		imgI_obj = nib.load(img_imag)
	except:
		print('')
		print('ERROR: the file storing the imaginary channel {} does not exist or is not in NIFTI format. Exiting with 1.'.format(img_imag))					 
		print('')
		sys.exit(1)

	imgI_data = imgI_obj.get_fdata()
	imgI_size = imgI_data.shape
	imgI_size = np.array(imgI_size)
	imgI_ndim = imgI_size.size
	imgI_data = np.array(imgI_data,'float64')
	imgI_header = imgI_obj.header
	imgI_affine = imgI_header.get_best_affine()

	# Check consistency of real and imaginay MRIs
	if ((imgR_ndim>4) or (imgR_ndim<2) or (imgI_ndim>4) or (imgI_ndim<2)):
		print('')
		print('ERROR: the input files {} and {} cannot have more than 4 dimensions and less than 2. Exiting with 1.'.format(img_real,img_imag))					 
		print('')
		sys.exit(1)

	if imgR_ndim!=imgI_ndim:
		print('')
		print('ERROR: the input files {} is {}D while the input file {} is {}D. Exiting with 1.'.format(img_real,imgR_ndim,img_imag,imgI_ndim))					 
		print('')
		sys.exit(1)
	
	if imgR_ndim==4:
		if imgR_size[3]!=imgI_size[3]:
			print('')
			print('ERROR: the input files {} and {} store a different number of measurements along the 4th dimension. Exiting with 1.'.format(img_real,img_imag))					 
			print('')
			sys.exit(1)

	if ( (np.sum(imgI_affine==imgR_affine)!=16) or (imgI_size[0]!=imgR_size[0]) or (imgI_size[1]!=imgR_size[1]) ):
		print('')
		print('ERROR: the geometry of the input files {} and {} do not match. Exiting with 1.'.format(img_real,img_imag))					 
		print('')
		sys.exit(1)

	if imgR_ndim>2:
		if imgI_size[2]!=imgR_size[2]:
			print('')
			print('ERROR: the geometry of the input files {} and {} do not match. Exiting with 1.'.format(img_real,img_imag))					 
			print('')
			sys.exit(1)


	# Load kernel
	if kernel=='B3':

		# Boxcar 3x3
		kernel_weights = np.array([[1.0/9.0, 1.0/9.0, 1.0/9.0],
					   [1.0/9.0, 1.0/9.0, 1.0/9.0],
					   [1.0/9.0, 1.0/9.0, 1.0/9.0]],'float64')

	elif kernel=='B5':

		# Boxcar 5x5
		kernel_weights = np.array([[1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0],
					   [1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0],
					   [1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0],
					   [1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0],
					   [1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0, 1.0/25.0]],'float64') 

	elif kernel=='G3F1':

		# Gaussian 3x3 with sigma = 1 voxel
		kernel_weights = np.array([[0.075113607954111, 0.123841403152974, 0.075113607954111],
   					   [0.123841403152974, 0.204179955571658, 0.123841403152974],
   					   [0.075113607954111, 0.123841403152974, 0.075113607954111]],'float64') 

	elif kernel=='G5F2':

		# Gaussian 5x5 with sigma = 2 voxels
		kernel_weights = np.array([[0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 0.023246839878294],
   					   [0.033823952439922, 0.049213560408541, 0.055766269846849, 0.049213560408541, 0.033823952439922],
   					   [0.038327559383904, 0.055766269846849, 0.063191462410265, 0.055766269846849, 0.038327559383904],
   					   [0.033823952439922, 0.049213560408541, 0.055766269846849, 0.049213560408541, 0.033823952439922],
   					   [0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 0.023246839878294]],'float64')

	elif kernel=='G3F1H':

		# Gaussian 3x3 with sigma = 1 voxel and center coefficient equal to 0
		kernel_weights = np.array([[0.075113607954111, 0.123841403152974, 0.075113607954111],
   					   [0.123841403152974, 0.0, 0.123841403152974],
   					   [0.075113607954111, 0.123841403152974, 0.075113607954111]],'float64') 

	elif kernel=='G5F2H':

		# Gaussian 5x5 with sigma = 2 voxels and center coefficient equal to 0
		kernel_weights = np.array([[0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 0.023246839878294],
   					   [0.033823952439922, 0.049213560408541, 0.055766269846849, 0.049213560408541, 0.033823952439922],
   					   [0.038327559383904, 0.055766269846849, 0.0, 0.055766269846849, 0.038327559383904],
   					   [0.033823952439922, 0.049213560408541, 0.055766269846849, 0.049213560408541, 0.033823952439922],
   					   [0.023246839878294, 0.033823952439922, 0.038327559383904, 0.033823952439922, 0.023246839878294]],'float64')

	elif kernel=='Opt3':

		# Optimised kernel 3x3
		kernel_weights = np.array([[0.107235538162453, 0.142764461837547, 0.107235538162453],
   					   [0.142764461837547, 0.0, 0.142764461837547],
   					   [0.107235538162453, 0.142764461837547, 0.107235538162453]],'float64') 

	elif kernel=='Opt5':

		# Optimised kernel 5x5
		kernel_weights = np.array([[0.025441320175391, 0.037016902431746, 0.041945645727859, 0.037016902431746, 0.025441320175391],
  	 				   [0.037016902431746, 0.053859275233950, 0.054719953999307, 0.053859275233950, 0.037016902431746],
   					   [0.041945645727859, 0.054719953999307, 0.0, 0.054719953999307, 0.041945645727859],
   					   [0.037016902431746, 0.053859275233950, 0.054719953999307, 0.053859275233950, 0.037016902431746],
   					   [0.025441320175391, 0.037016902431746, 0.041945645727859, 0.037016902431746, 0.025441320175391]],'float64') 

	else:
		print('')
		print('ERROR: the kernel {} is not supported. Exiting with 1.'.format(kernel))					 
		print('')
		sys.exit(1)

	### Filter data with the specified kernel
	if imgR_ndim==2:

		# Filter real and imaginary channels independently
		imgR_data_filt = ndimage.convolve(imgR_data, kernel_weights, mode='constant', cval=0.0)
		imgI_data_filt = ndimage.convolve(imgI_data, kernel_weights, mode='constant', cval=0.0)
		
	elif imgR_ndim==3:
		
		# Filter real and imaginary channels independently
		imgI_data_filt = np.zeros(imgR_size,'float64')
		imgR_data_filt = np.zeros(imgR_size,'float64')
		for zz in range(0, imgR_size[2]):
			imgI_data_filt[:,:,zz] = ndimage.convolve(imgI_data[:,:,zz], kernel_weights, mode='constant', cval=0.0)
			imgR_data_filt[:,:,zz] = ndimage.convolve(imgR_data[:,:,zz], kernel_weights, mode='constant', cval=0.0)

	elif imgR_ndim==4:

		# Filter real and imaginary channels independently
		imgI_data_filt = np.zeros(imgR_size,'float64')
		imgR_data_filt = np.zeros(imgR_size,'float64')
		for vv in range(0, imgR_size[3]):
			for zz in range(0, imgR_size[2]):
				imgI_data_filt[:,:,zz,vv] = ndimage.convolve(imgI_data[:,:,zz,vv], kernel_weights, mode='constant', cval=0.0)
				imgR_data_filt[:,:,zz,vv] = ndimage.convolve(imgR_data[:,:,zz,vv], kernel_weights, mode='constant', cval=0.0)

	
	### Get phase of complex data after filtering
	phase_data_orig = np.angle(imgR_data + imgI_data*1j)
	phase_data_filt = np.angle(imgR_data_filt + imgI_data_filt*1j)
		
	### Rephase measured signals so that the true information is in the real channel only; for the rephasing, use the phase of the signal after filtering
	rephased_data_complex = (imgR_data + 1j*imgI_data)*(np.exp(-1j*phase_data_filt))    # Rephase signals
	rephased_data_R = np.real(rephased_data_complex)    # Get real channel of rephased signals (this should theoretically contain only true information)
	rephased_data_I = np.imag(rephased_data_complex)    # Get imaginary channel of rephased signals (this should theoretically contain only Gaussian noise)
	rephased_data_M = np.sqrt(rephased_data_R*rephased_data_R + rephased_data_I*rephased_data_I)   # Get magnitude of rephased signals (when this differs too much from rephased_data_R, then the rephasing has probably gone wrong)
	rephased_data_deltaMR = np.abs(rephased_data_M - rephased_data_R)   # Difference between magnitude and real channel
	phase_data_new = np.angle(rephased_data_R + rephased_data_I*1j)     # Phase after rephasing
	

	### Clear some memory
	del imgI_data, imgR_data

	### Calculate noise level and remove outliers (look at MAD within a window the same size as the kernels)
	rephased_data_R_thresh = rephased_data_M
	rephased_data_I_thresh = np.zeros(imgR_size,'float64')
	outliers_flag = np.ones(imgR_size,'float64')
	if imgR_ndim==2:

		absdev = np.abs(rephased_data_I - signal.medfilt(rephased_data_I,kernel_weights.shape))   # Absolute deviation of imaginary channel within kernel window
		medabsdev = signal.medfilt(absdev,kernel_weights.shape)    # Median absolute deviation of imaginary channel within kernel window
		thresh = 2.5000*1.4826*medabsdev   # Local threhsold
		rephased_data_R_thresh[rephased_data_deltaMR<thresh] = rephased_data_R[rephased_data_deltaMR<thresh]
		rephased_data_I_thresh[rephased_data_deltaMR<thresh] = rephased_data_I[rephased_data_deltaMR<thresh]
		outliers_flag[rephased_data_deltaMR<thresh] = 0.0

	elif imgR_ndim==3:
	
		thresh = np.zeros(imgR_size,'float64')	
		for zz in range(0, imgR_size[2]):

			absdev_slice = np.abs(rephased_data_I[:,:,zz] - signal.medfilt(rephased_data_I[:,:,zz],kernel_weights.shape))   # Absolute deviation of imaginary channel within kernel window
			medabsdev_slice = signal.medfilt(absdev_slice,kernel_weights.shape)    # Median absolute deviation of imaginary channel within kernel window
			thresh_slice = 2.5000*1.4826*medabsdev_slice   # Local threhsold
			thresh[:,:,zz] = thresh_slice
		 
		rephased_data_R_thresh[rephased_data_deltaMR<thresh] = rephased_data_R[rephased_data_deltaMR<thresh]
		rephased_data_I_thresh[rephased_data_deltaMR<thresh] = rephased_data_I[rephased_data_deltaMR<thresh]
		outliers_flag[rephased_data_deltaMR<thresh] = 0.0


	elif imgR_ndim==4:

		thresh = np.zeros(imgR_size,'float64')	
		for vv in range(0, imgR_size[3]):
			for zz in range(0, imgR_size[2]):

				absdev_vol_slice = np.abs(rephased_data_I[:,:,zz,vv] - signal.medfilt(rephased_data_I[:,:,zz,vv],kernel_weights.shape))   # Absolute deviation of imaginary channel within kernel window
				medabsdev_vol_slice = signal.medfilt(absdev_vol_slice,kernel_weights.shape)    # Median absolute deviation of imaginary channel within kernel window
				thresh_vol_slice = 2.5000*1.4826*medabsdev_vol_slice   # Local threhsold
				thresh[:,:,zz,vv] = thresh_vol_slice
		 
		rephased_data_R_thresh[rephased_data_deltaMR<thresh] = rephased_data_R[rephased_data_deltaMR<thresh]
		rephased_data_I_thresh[rephased_data_deltaMR<thresh] = rephased_data_I[rephased_data_deltaMR<thresh]
		outliers_flag[rephased_data_deltaMR<thresh] = 0.0

	phase_data_new_thresh = np.angle(rephased_data_R_thresh + rephased_data_I_thresh*1j)

	### Save as real and imaginary channels after rephasing and after rephasing + outlier detection as NIFTI
	
	# Create file names
	buffer_string=''
	seq_string = (out_base,'_RealReph.nii')
	rephased_R_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_RealRephThresh.nii')
	rephased_R_thresh_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_ImagReph.nii')
	rephased_I_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_ImagRephThresh.nii')
	rephased_I_thresh_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_OutlierDetected.nii')
	flag_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseOriginal.nii')
	phaseorig_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseRephased.nii')
	phasenew_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseBackground.nii')
	phaseest_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseRephasedOutliers.nii')
	phasenewthresh_outfile = buffer_string.join(seq_string)

	# Create header
	buffer_header = imgR_header
	buffer_header.set_data_dtype('float64')   # Make sure we save output data as float64, even if input header indicates a different data type
	
	# Save files
	rephased_obj_R = nib.Nifti1Image(rephased_data_R,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_R, rephased_R_outfile)
	
	rephased_obj_R_thresh = nib.Nifti1Image(rephased_data_R_thresh,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_R_thresh, rephased_R_thresh_outfile)

	rephased_obj_I = nib.Nifti1Image(rephased_data_I,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_I, rephased_I_outfile)
	
	rephased_obj_I_thresh = nib.Nifti1Image(rephased_data_I_thresh,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_I_thresh, rephased_I_thresh_outfile)

	flag_outobj = nib.Nifti1Image(outliers_flag,imgR_obj.affine,buffer_header)
	nib.save(flag_outobj, flag_outfile)	

	phaseorig_obj = nib.Nifti1Image(phase_data_orig,imgR_obj.affine,buffer_header)
	nib.save(phaseorig_obj, phaseorig_outfile)

	phasenew_obj = nib.Nifti1Image(phase_data_new,imgR_obj.affine,buffer_header)
	nib.save(phasenew_obj, phasenew_outfile)

	phasenewthresh_obj = nib.Nifti1Image(phase_data_new_thresh,imgR_obj.affine,buffer_header)
	nib.save(phasenewthresh_obj, phasenewthresh_outfile)

	phaseest_obj = nib.Nifti1Image(phase_data_filt,imgR_obj.affine,buffer_header)
	nib.save(phaseest_obj, phaseest_outfile)



# Run the module as a script when required
if __name__ == "__main__":

	### Parse arguments or print help
	parser = argparse.ArgumentParser(description='Rephasing of complex MR images with noise decorrelation according to Sprenger T et al, MRM 2017, 77:559-570. Author: Francesco Grussu, University College London. Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>. Code released under BSD Two-Clause license. Copyright (c) 2020 University College London. All rights reserved.')
	parser.add_argument('img_real', help='3D or 4D Nifti file storing the real channel in image space')
	parser.add_argument('img_imag', help='3D or 4D Nifti file storing the imaginary channel in image space')
	parser.add_argument('out_base', help='base name of the output files (output files are: *_RealReph.nii, storing the real channel rephased; *_RealRephThresh.nii, storing the real channel rephased with outlier detection; *_ImagReph.nii, storing the imaginary channel rephased; *_ImagRephThresh.nii, storing the real channel rephased with outlier detection; *_OutlierDetected.nii, flagging with 1 outliers; *_PhaseOriginal.nii, storing the original phase; *_PhaseBackground.nii, storing the estimated background phase; *_PhaseRephased.nii, storing the phase after rephasing using the background phase; *_PhaseRephasedOutliers.nii, storing the original phase after rephasing where outliers are set to zero phase)')
	parser.add_argument('kernel', help='kernel for decorrelation filers (choose among B3, B5, G3F1, G5F2, G3F1H, G5F2H, Opt3 and Opt5; see Sprenger T et al, MRM 2017, 77:559-570)')
	args = parser.parse_args()

	### Rephase the data
	rephaseDC(args.img_real,args.img_imag,args.out_base,args.kernel)  
    
