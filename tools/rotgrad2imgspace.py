### Rotate diffusion gradient directions from scanner space to image space
#   The code attempts to replicate the behaviour of the freely available dcm2niix 
#   DICOM to NIFTI converter (http://github.com/rordenlab/dcm2niix)
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

import numpy as np
import nibabel as nb
import argparse
import math
import sys
from nibabel import quaternions
from nibabel import affines

### Print help and parse arguments
parser = argparse.ArgumentParser(description='Rotate diffusion gradient directions from scanner space to image space. The code attempts to replicate the behaviour of the freely available dcm2niix DICOM to NIFTI converter (http://github.com/rordenlab/dcm2niix). Author: Francesco Grussu, University College London <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>')
parser.add_argument('grads_scanner', help='Text file containing the list of gradient directions in scanner space. Components along lines, different directions along columns (the file must contain 3 lines and as many columns as number of directions, with each entry separated by a space)')
parser.add_argument('input_nifti', help='Nifti file to the voxel space of which gradient directions will be rotated')
parser.add_argument('grads_img', help='Path of the output text file that will store the gradient directions rotated to the image space, using the FSL format')
parser.add_argument('--meth', metavar='METHOD', default='ROT', help='Method to be employed to rotate the gradients (ROT: rotation matrices; QUAT: quaternion multiplication; AFF: header affine transformation; default is ROT)')
args = parser.parse_args()

### Print information
print('***********************************')
print('        rotgrad2imgspace.py        ')
print('***********************************')
print('')
print('Gradients in scanner space to be rotated: {}'.format(args.grads_scanner))
print('Reference NIFTI file: {}'.format(args.input_nifti))
print('Output file with gradients rotated in image space: {}'.format(args.grads_img))
print('Rotation method: {}'.format(args.meth))
print('')

### Get input parameters and load input data
g_scanner = np.loadtxt(args.grads_scanner)
g_scanner = np.array(g_scanner,'float64')
img = nb.load(args.input_nifti)
g_imgfile = args.grads_img

### Rotate gradients according to the specified algorithm

# Rotate gradients using rotation matricies constructed from the quaterion elements
if args.meth=='ROT':

	# Construct rotation matrix from header
	img_header = img.header
	img_header_raw = img_header.structarr
	pixdim = img_header_raw['pixdim'];
	qfac = float(pixdim[0])
	if qfac==0:
	    qfac = 1.0
	qb = float(img_header_raw['quatern_b'])
	qc = float(img_header_raw['quatern_c'])
	qd = float(img_header_raw['quatern_d'])
	qa = math.sqrt(1 - (qb*qb  +  qc*qc  +  qd*qd))
	R11 = qa*qa + qb*qb - qc*qc - qd*qd
	R12 = 2*qb*qc - 2*qa*qd
	R13 = 2*qb*qd + 2*qa*qc
	R21 = 2*qb*qc + 2*qa*qd
	R22 = qa*qa + qc*qc - qb*qb - qd*qd
	R23 = 2*qc*qd - 2*qa*qb
	R31 = 2*qb*qd - 2*qa*qc
	R32 = 2*qc*qd + 2*qa*qb
	R33 = qa*qa + qd*qd - qc*qc - qb*qb
	R = np.array([[R11, R12, R13],
		      [R21, R22, R23],
		      [R31, R32, R33]],'float64') 
	Rt = np.transpose(R)

	# Apply rotation matrix to rotate gradients from scanner to image space
	g_dims = g_scanner.shape
	g_img = np.zeros(g_dims)
	Ngrads = g_img.shape[1]

	for nn in range(0,Ngrads):
	    g_buff = np.array([g_scanner[0,nn], g_scanner[1,nn], g_scanner[2,nn]/qfac]) 
	    g_nn_scanner = np.reshape(g_buff,(3,1))
	    g_nn_img = np.matmul(Rt,g_nn_scanner)
	    g_img[0,nn] = g_nn_img[0]
	    g_img[1,nn] = g_nn_img[1]
	    g_img[2,nn] = g_nn_img[2]/qfac

# Rotate gradients using quaterion algebra
elif args.meth=='QUAT':

	# Construct quaterion from header
	img_header = img.header
	img_header_raw = img_header.structarr
	pixdim = img_header_raw['pixdim'];
	qfac = np.array(pixdim[0],'float64')
	if qfac==0:
	    qfac = 1.0
	qb = float(img_header_raw['quatern_b'])
	qc = float(img_header_raw['quatern_c'])
	qd = float(img_header_raw['quatern_d'])
	qa = math.sqrt(1 - (qb*qb  +  qc*qc  +  qd*qd))
	qinv = quaternions.inverse([qa, qb, qc, qd]) 

	# Apply quaterion to rotate gradients from scanner to image space
	g_dims = g_scanner.shape
	g_img = np.zeros(g_dims)
	Ngrads = g_img.shape[1]

	for nn in range(0,Ngrads):
	    q_nn_scanner = [0.0, g_scanner[0,nn], g_scanner[1,nn], g_scanner[2,nn]/qfac] 
	    qleft_nn = [qinv[0], qinv[1], qinv[2], qinv[3]]
	    qright_nn = [qinv[0], (-1.0)*qinv[1], (-1.0)*qinv[2], (-1.0)*qinv[3]]
	    q_nn_image = quaternions.mult(qleft_nn,quaternions.mult(q_nn_scanner,qright_nn))
	    g_img[0,nn] = q_nn_image[1]
	    g_img[1,nn] = q_nn_image[2]
	    g_img[2,nn] = q_nn_image[3]/qfac

# Rotate gradients using the affine transformation stored in the header
elif args.meth=='AFF':

	# Construct affine transformation from header
	img_header = img.header
	img_header_raw = img_header.structarr
	pixdim = img_header_raw['pixdim'];
	qfac = np.array(pixdim[0],'float64')
	if qfac==0:
	    qfac = 1.0
	mymat = img.affine
	Aff = mymat[0:3,0:3]
	invAff = np.linalg.inv(Aff)

	# Apply affine transformations to rotate gradients from scanner to image space
	g_dims = g_scanner.shape
	g_img = np.zeros(g_dims)
	Ngrads = g_img.shape[1]

	for nn in range(0,Ngrads):
	    g_buff = np.array([g_scanner[0,nn], g_scanner[1,nn], g_scanner[2,nn]/qfac]) 
	    g_nn_scanner = np.reshape(g_buff,(3,1))
	    tmp_nn_img = np.matmul(invAff,g_nn_scanner)
	    g_img[0,nn] = tmp_nn_img[0]*pixdim[1]
	    g_img[1,nn] = tmp_nn_img[1]*pixdim[2]
	    g_img[2,nn] = tmp_nn_img[2]*pixdim[3]


else:
	print('')
	print('ERROR: the method {} is not supported. Exiting with 1.'.format(args.meth))					 
	print('')
	sys.exit(1)


# Save gradient directions rotated in image space to the output text file
np.savetxt(g_imgfile,g_img,fmt='%.6f')

# Exit with no error
print('Done')
print('')
sys.exit(0)

