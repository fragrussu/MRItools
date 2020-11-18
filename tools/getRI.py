### Get real and imaginary MRI signals from magnitude and phase
#
#
# Author: Francesco Grussu, UCL Queen Square Institute of Neurology
#                           UCL Centre for Medical Image Computing
#         <f.grussu@ucl.ac.uk>
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
#

import nibabel as nib
import numpy as np
import sys, argparse

def imgRI(phase_nifti,mag_nifti,real_nifti,imag_nifti):
	'''Calculate phase and magnitude of complex MR images given real and imaginary images
	
	   INTERFACE
	   imgRI(phase_nifti,mag_nifti,real_nifti,imag_nifti)
	
	   PARAMETERS
	   phase_nifti: 2D, 3D or 4D Nifti storing the phase in image space
           mag_nifti: 2D, 3D or 4D Nifti storing the magnitude in image space
	   real_nifti: path of output file containing the real channel
	   imag_nifti: path of output file containing the imaginary channel
	   
	   NOTES
	   - phase_nifti and mag_nifti must have the same size and same header geometry
	   - real_nifti and imag_nifti will store an image with the same dimensions of mag_nifti and
	     phase_nifti (voxel-wise calculation)
           - phase assumed to be in radians
	
	   Author: Francesco Grussu, University College London
		  <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''


	### Load input data
	
	# Load phase
	try:
		imgP_obj = nib.load(phase_nifti)
	except:
		print('')
		print('ERROR: the file storing the phase {} does not exist or is not in NIFTI format. Exiting with 1.'.format(phase_nifti))
		print('')
		sys.exit(1)
	
	imgP_data = imgP_obj.get_fdata()
	imgP_size = imgP_data.shape
	imgP_size = np.array(imgP_size)
	imgP_ndim = imgP_size.size
	imgP_data = np.array(imgP_data,'float64')
	imgP_header = imgP_obj.header
	imgP_affine = imgP_header.get_best_affine()	

	# Load magnitude
	try:
		imgM_obj = nib.load(mag_nifti)
	except:
		print('')
		print('ERROR: the file storing the magnitude {} does not exist or is not in NIFTI format. Exiting with 1.'.format(mag_nifti))
		print('')
		sys.exit(1)

	imgM_data = imgM_obj.get_fdata()
	imgM_size = imgM_data.shape
	imgM_size = np.array(imgM_size)
	imgM_ndim = imgM_size.size
	imgM_data = np.array(imgM_data,'float64')
	imgM_header = imgM_obj.header
	imgM_affine = imgM_header.get_best_affine()

	# Check consistency of real and imaginay MRIs
	if ((imgP_ndim>4) or (imgP_ndim<2) or (imgM_ndim>4) or (imgM_ndim<2)):
		print('')
		print('ERROR: the input files {} and {} cannot have more than 4 dimensions and less than 2. Exiting with 1.'.format(phase_nifti,mag_nifti))	
		print('')
		sys.exit(1)

	if imgP_ndim!=imgM_ndim:
		print('')
		print('ERROR: the input files {} is {}D while the input file {} is {}D. Exiting with 1.'.format(phase_nifti,imgP_ndim,mag_nifti,imgP_ndim))
		print('')
		sys.exit(1)
	
	if imgP_ndim==4:
		if imgP_size[3]!=imgM_size[3]:
			print('')
			print('ERROR: the input files {} and {} store a different number of measurements along the 4th dimension. Exiting with 1.'.format(phase_nifti,mag_nifti))
			print('')
			sys.exit(1)

	if ( (np.sum(imgP_affine==imgM_affine)!=16) or (imgP_size[0]!=imgM_size[0]) or (imgP_size[1]!=imgM_size[1]) ):
		print('')
		print('ERROR: the geometry of the input files {} and {} do not match. Exiting with 1.'.format(phase_nifti,mag_nifti))
		print('')

	if imgP_ndim>2:
		if imgM_size[2]!=imgP_size[2]:
			print('')
			print('ERROR: the geometry of the input files {} and {} do not match. Exiting with 1.'.format(phase_nifti,mag_nifti))
			print('')

	### Get real and imaginary parts
	complex_data = imgM_data*np.exp(1j*imgP_data)
	real_data = np.real(complex_data)
	imag_data = np.imag(complex_data)

	### Save output	

	# Create header
	buffer_header = imgP_header
	buffer_header.set_data_dtype('float64')   # Make sure we save output data as float64, even if input header indicates a different data type
	# Save files
	out_obj1 = nib.Nifti1Image(real_data,imgP_obj.affine,buffer_header)
	nib.save(out_obj1, real_nifti)
	out_obj2 = nib.Nifti1Image(imag_data,imgP_obj.affine,buffer_header)
	nib.save(out_obj2, imag_nifti)




# Run the module as a script when required
if __name__ == "__main__":

	### Parse arguments or print help
	parser = argparse.ArgumentParser(description='Calculate real and imaginary parts of complex-valued MR data in image space. Author: Francesco Grussu, University College London. Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>.')
	parser.add_argument('img_phase', help='2D, 3D or 4D Nifti file storing the phase in image space')
	parser.add_argument('img_mag', help='2D, 3D or 4D Nifti file storing the magnitude in image space')
	parser.add_argument('img_real', help='path of the output file storing the real channel')
	parser.add_argument('img_imag', help='path of the output file storing the imaginary channel')
	args = parser.parse_args()

	### Get real and imaginary
	imgRI(args.img_phase,args.img_mag,args.img_real,args.img_imag)
	sys.exit(0)
	
