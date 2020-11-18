### Get MR phase from real and imaginary MR images
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

def imgphase(real_nifti,imag_nifti,phase_nifti):
	'''Calculate phase of complex MR images given real and imaginary images
	
	   INTERFACE
	   imgphase(real_nifti,imag_nifti,phase_nifti)
	
	   PARAMETERS
	   real_nifti: 2D, 3D or 4D Nifti storing the real channel in image space
           imag_nifti: 2D, 3D or 4D Nifti storing the imaginary channel in image space
	   phase_nifti: path of output file containing phase
	
	   NOTES
	   - imag_nifti and real_nifti must have the same size and same header geometry
	   - phase_nifti will store an image with the same dimensions of real_nifti and
	     imag_nifti (voxel-wise phase calculation)
	
	   Author: Francesco Grussu, University College London, July 2017
		<f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''


	### Load input data
	
	# Load real MRI
	try:
		imgR_obj = nib.load(real_nifti)
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
		imgI_obj = nib.load(imag_nifti)
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

	if imgR_ndim>2:
		if imgI_size[2]!=imgR_size[2]:
			print('')
			print('ERROR: the geometry of the input files {} and {} do not match. Exiting with 1.'.format(img_real,img_imag))					 
			print('')

	### Get phase of complex data
	phase_data = np.angle(imgR_data + imgI_data*1j)

	### Save output	
	# Create header
	buffer_header = imgR_header
	buffer_header.set_data_dtype('float64')   # Make sure we save output data as float64, even if input header indicates a different data type
	# Save files
	out_obj = nib.Nifti1Image(phase_data,imgR_obj.affine,buffer_header)
	nib.save(out_obj, phase_nifti)




# Run the module as a script when required
if __name__ == "__main__":

	### Parse arguments or print help
	parser = argparse.ArgumentParser(description='Calculate phase of complex MR data in image space. Author: Francesco Grussu, University College London. Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>.')
	parser.add_argument('img_real', help='2D, 3D or 4D Nifti file storing the real channel in image space')
	parser.add_argument('img_imag', help='2D, 3D or 4D Nifti file storing the imaginary channel in image space')
	parser.add_argument('img_phase', help='path of the output file storing the phase')
	args = parser.parse_args()

	### Rephase the data
	imgphase(args.img_real,args.img_imag,args.img_phase)
	sys.exit(0)
	
