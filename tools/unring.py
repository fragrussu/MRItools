# Author: Francesco Grussu, University College London
#		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
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
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.


### Load useful modules
import argparse, os, sys
import nibabel as nib
import numpy as np
import warnings
from pathlib import Path as pt
sys.path.insert(0, os.path.dirname(pt(__file__).absolute()) )
import gibbs_removal as gbr   # From: https://github.com/RafaelNH/gibbs-removal/blob/master/gibbs_removal.py



# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Wrapper enabling the application of MRI Gibbs unringing from command line, based on the implementation by Neto Henriques R (https://github.com/RafaelNH/gibbs-removal) of the algorithm by Kellner et al. References: Kellner et al. Gibbs‐ringing artifact removal based on local subvoxel‐shifts (2016), Magn Reson Med 76:1574–1581; Neto Henriques, Advanced Methods for Diffusion MRI Data Analysis and their Application to the Healthy Ageing Brain (2018, Doctoral thesis, DOI: 10.17863/CAM.29356). Email: <francegrussu@gmail.com> <f.grussu@ucl.ac.uk>.')
	parser.add_argument('infile', help='Nifti file storing the volume to unring (it can be 2D, 3D or 4D)')
	parser.add_argument('outfile', help='Nifti file storing the unringed output volume (same size as input)')
	parser.add_argument('--np', metavar='<value>', default=9, help='Number of neighbouring points for local total variation calculation (default: 9)')
	parser.add_argument('--axis', metavar='<value>', default=2, help='Axis along which slices are acquired (choose among 0, 1, 2 if data is 3D or among 0, 1, 2, 3 if data is 4D; default: 2, i.e. slices stacked along third dimension)')
	args = parser.parse_args()

	### Print some info
	print('')
	print('********************************************************************')
	print('                                Unring                              ')
	print('********************************************************************')
	print('')
	print('* Input data: {}'.format(args.infile))
	print('* Output data: {}'.format(args.outfile))
	print('* Neighbourhood size for TV calculation: {}'.format(args.np))
	print('* Slice axis: {}'.format(args.axis))
	print('')


	### Get input MRI data
	print('    ... loading input data')
	try:
		mri_obj = nib.load(args.infile)
	except:
		print('')
		raise RuntimeError('ERROR: the input MRI file {} does not exist or is not in NIFTI format.'.format(infile))	 
	mri_data = mri_obj.get_fdata()
	mri_data[np.isnan(mri_data)] = 0.0    # Remove NaNs
	mri_data[np.isinf(mri_data)] = 0.0    # Remove Infs

	### Unring
	print('    ... unringing')
	np = int(args.np)
	axdir = int(args.axis)
	mri_data_corr = gbr.gibbs_removal(mri_data, axdir, np)

	### Save output data
	print('    ... saving corrected data')
	try:
		buffer_header = mri_obj.header
		buffer_header.set_data_dtype('float64')   # Make sure we save output as a float64
		out_obj = nib.Nifti1Image(mri_data_corr,mri_obj.affine,buffer_header)
		nib.save(out_obj,args.outfile)
	except:
		raise RuntimeError('ERROR: the output folder does not exist or you do not have permissions to write there.')

	### Done
	print('    ... done')
	print('')

