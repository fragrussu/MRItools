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
from skimage.restoration import denoise_tv_bregman

def rephaseTV(img_real,img_imag,out_base,lam_value):
	''' Rephase complex MR data in image space
	
	INTERFACE
	rephaseTV(img_real,img_imag,out_root,lam_value)

	PARAMETERS
	img_real: path of a NIFTI file storing a 3D or 4D image (real channel)
	img_imag: path of a NIFTI file storing a 3D or 4D image (imaginary channel)
	out_base: base name of output files (the output files contain the real and imaginary channels 
		  after noise decorrelation and rephasing; these will end in *_RealReph.nii 
		  (real channel rephased), *_ImagReph.nii (imaginary channel rephased), 
	          *_PhaseOriginal.nii (storing the original phase), _*PhaseBackground.nii (storing the estimated 
		  background phase), *_PhaseRephased.nii (storing the phase after rephasing using the 
		  background phase). Note that the imaginary channel after rephasing should contain 
		  mostly noise and negligible true signal information. 
	lam_value: weight of total variation regulariser (see Eichner C et al, NeuroImage 2015, 122:373-384).

	DESCRIPTION
	The function the algorithm based on total variation denoising of Eichner C et al, 
	NeuroImage 2015, 122:373-384. The function works with 3D and 4D NIFTI files (in the latter case,
	each volume of the 4D NIFTI is treated independently).

	References: "Real diffusion-weighted MRI enabling true signal averaging and 
		     increased diffusion contrast.", Eichner C et al, NeuroImage 2015, 122:373-384
	     
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


	# Weight of regulariser in TV denoising
	try:
		lam_value = float(lam_value)
	except:
		print('')
		print('ERROR: the requested weight of TV regulariser does not appear to be numeric. Exiting with 1.')					 
		print('')
		sys.exit(1)

	### Filter data with the specified kernel
	imgC_data = imgR_data + imgI_data*1j
	imgP_data = np.angle(imgC_data)
	imgM_data = np.sqrt(imgR_data*imgR_data + imgI_data*imgI_data) 
	phase_data_filt = np.zeros(imgR_size,'float64')
	
	if imgR_ndim==2:

		# Use multichannel TV denoising
		buffimg = denoise_tv_bregman(np.dstack((imgM_data,imgP_data)), 2.0*lam_value)
		phase_data_filt = buffimg[:,:,1]
		
	elif imgR_ndim==3:
		
		# Use multichannel TV denoising
		for zz in range(0, imgR_size[2]):
			buffimg = denoise_tv_bregman(np.dstack((imgM_data[:,:,zz],imgP_data[:,:,zz])), 2.0*lam_value)
			phase_data_filt[:,:,zz] = buffimg[:,:,1]

	elif imgR_ndim==4:

		# Use multichannel TV denoising
		for vv in range(0, imgR_size[3]):
			for zz in range(0, imgR_size[2]):
				buffimg = denoise_tv_bregman(np.dstack((imgM_data[:,:,zz,vv],imgP_data[:,:,zz,vv])), 2.0*lam_value)
				phase_data_filt[:,:,zz,vv] = buffimg[:,:,1]

		
	### Rephase measured signals so that the true information is in the real channel only; for the rephasing, use the phase of the signal after filtering
	rephased_data_complex = imgC_data*(np.exp(-1j*phase_data_filt))    # Rephase signals
	rephased_data_R = np.real(rephased_data_complex)    # Get real channel of rephased signals (this should theoretically contain only true information)
	rephased_data_I = np.imag(rephased_data_complex)    # Get imaginary channel of rephased signals (this should theoretically contain only Gaussian noise)
	phase_data_new = np.angle(rephased_data_R + rephased_data_I*1j)     # Phase after rephasing
	

	### Clear some memory
	del imgI_data, imgR_data

	### Save as real and imaginary channels after rephasing and after rephasing + outlier detection as NIFTI
	
	# Create file names
	buffer_string=''
	seq_string = (out_base,'_RealReph.nii')
	rephased_R_outfile = buffer_string.join(seq_string)


	buffer_string=''
	seq_string = (out_base,'_ImagReph.nii')
	rephased_I_outfile = buffer_string.join(seq_string)


	buffer_string=''
	seq_string = (out_base,'_PhaseOriginal.nii')
	phaseorig_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseRephased.nii')
	phasenew_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (out_base,'_PhaseBackground.nii')
	phaseest_outfile = buffer_string.join(seq_string)


	# Create header
	buffer_header = imgR_header
	buffer_header.set_data_dtype('float64')   # Make sure we save output data as float64, even if input header indicates a different data type
	
	# Save files
	rephased_obj_R = nib.Nifti1Image(rephased_data_R,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_R, rephased_R_outfile)
	
	rephased_obj_I = nib.Nifti1Image(rephased_data_I,imgR_obj.affine,buffer_header)
	nib.save(rephased_obj_I, rephased_I_outfile)

	phaseorig_obj = nib.Nifti1Image(imgP_data,imgR_obj.affine,buffer_header)
	nib.save(phaseorig_obj, phaseorig_outfile)

	phasenew_obj = nib.Nifti1Image(phase_data_new,imgR_obj.affine,buffer_header)
	nib.save(phasenew_obj, phasenew_outfile)

	phaseest_obj = nib.Nifti1Image(phase_data_filt,imgR_obj.affine,buffer_header)
	nib.save(phaseest_obj, phaseest_outfile)



# Run the module as a script when required
if __name__ == "__main__":

	### Parse arguments or print help
	parser = argparse.ArgumentParser(description='Rephasing of complex MR images using total variation denoising according to Eichner C et al, NeuroImage 2015, 122:373-384. Author: Francesco Grussu, University College London. Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>. Code released under BSD Two-Clause license. Copyright (c) 2020 University College London. All rights reserved.')
	parser.add_argument('img_real', help='3D or 4D Nifti file storing the real channel in image space')
	parser.add_argument('img_imag', help='3D or 4D Nifti file storing the imaginary channel in image space')
	parser.add_argument('out_base', help='base name of the output files (output files are: *_RealReph.nii, storing the real channel rephased; *_ImagReph.nii, storing the imaginary channel rephased; *_PhaseOriginal.nii, storing the original phase; *_PhaseBackground.nii, storing the estimated background phase; *_PhaseRephased.nii, storing the phase after rephasing with the background phase)')
	parser.add_argument('--w', metavar='<value>', default='6.0', help='weight of total variation regulariser for split-Bregman total variation denoising (default 6.0; see Eichner C et al, NeuroImage 2015, 122:373-384)')
	args = parser.parse_args()

	### Rephase the data
	rephaseTV(args.img_real,args.img_imag,args.out_base,args.w)  
	sys.exit(0)
	
