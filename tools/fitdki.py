### Fitting of Diffusion Kurtosis Imaging on DWI data; wrapper of DiPy
#
#
# Author: Francesco Grussu, UCL Queen Square Institute of Neurology
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

# Load 	useful modules
import multiprocessing
import argparse, sys, os
import numpy as np
import nibabel as nib
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.dki as dki

# Run the module as a script
if __name__ == "__main__":

	# Print help and parse arguments
	parser = argparse.ArgumentParser(description='Run DKI fitting using DiPy. Command-line wrapper for the following DiPy tutorial: http://dipy.org/documentation/1.0.0./examples_built/reconst_dki . Author: Francesco Grussu, <f.grussu@ucl.ac.uk>.')
	parser.add_argument('dwi', help='DWI data set (4D Nifti)')
	parser.add_argument('out', help='Output base name. Output files will append to the base name "_fa.nii" (fractional anisotropy), "_ad.nii" (axial diffusivity, in um^2/ms), "_rd.nii" (radial diffusivity, in um^2/ms), "_md.nii" (mean diffusivity, in um^2/ms), "_ak.nii" (axial kurtosis excess), "_rk.nii" (radial kurtosis excess), "_mk.nii" (mean kurtosis excess).')
	parser.add_argument('bvals', help='b-value files in FSL format (1 x no. images), expressed in [s mm^-2]')
	parser.add_argument('bvecs', help='gradient directions in FSL format (3 x no. images)')	
	parser.add_argument('--mask', metavar='<file>', help='mask (Nifti; 1/0 in voxels to include/exclude)')
	parser.add_argument('--kmin', metavar='<value>', default=0, help='minimum value allowed for kurtosis excess values (default 0)')
	parser.add_argument('--kmax', metavar='<value>', default=3, help='maximum value allowed for kurtosis excess values (default 3)')
	args = parser.parse_args()

	# Get input parameters
	dwifile = args.dwi
	outbase = args.out
	bvalfile = args.bvals
	bvecfile = args.bvecs
	maskid = args.mask
	kmin = float(args.kmin)
	kmax = float(args.kmax)

	if(kmax<=kmin):
		print('')
		print('ERROR: maximum allowed kurtosis must be higher than minimum allowed kurtosis!')
		print('')
		sys.exit(1)

	# Print
	print('')
	print('*************************************************')
	print('    Fitting of DKI model (wrapper of DiPy)       ')
	print('*************************************************')
	print('')
	print('Called on 4D DWI Nifti file: {}'.format(dwifile))
	print('b-values: {}'.format(bvalfile))
	print('gradient directions: {}'.format(bvecfile))
	print('kurtosis excess range: [{}; {}]'.format(kmin,kmax))
	print('Output base name: {}'.format(outbase))
	print('')

	# Load DWI
	dwiobj = nib.load(dwifile)    # Input data
	dwiimg = dwiobj.get_fdata()
	imgsize = dwiimg.shape

	# Load b-values and gradient directions
	bvals, bvecs = read_bvals_bvecs(bvalfile, bvecfile)
	gtab = gradient_table(bvals, bvecs)

	# Load mask if necessary
	if isinstance(maskid, str)==1:

		# The user provided a mask: load it
		maskobj = nib.load(maskid)
		maskvals = np.array(maskobj.get_fdata(),dtype=bool)
	else:

		# No mask provided: flag that you want to analyse all voxels
		maskvals = np.array(np.ones(imgsize[0:3]),dtype=bool)

	# Starting DKI model fitting
	print('')
	print('... starting model fitting. Please wait...')
	print('')

	dkimodel = dki.DiffusionKurtosisModel(gtab)
	my_dki_fit = dkimodel.fit(dwiimg, mask=maskvals)

	# Get output metrics
	FA = np.array(my_dki_fit.fa)
	MD = 1000.000*np.array(my_dki_fit.md)    # Convert to um^2/ms
	AD = 1000.000*np.array(my_dki_fit.ad)    # Convert to um^2/ms
	RD = 1000.000*np.array(my_dki_fit.rd)    # Convert to um^2/ms
	MK = np.array(my_dki_fit.mk(kmin,kmax))
	AK = np.array(my_dki_fit.ak(kmin,kmax))
	RK = np.array(my_dki_fit.rk(kmin,kmax))
	dirs = np.array(my_dki_fit.directions)
	outvec = np.zeros((dirs.shape[0],dirs.shape[1],dirs.shape[2],3))
	outvec[:,:,:,0] = dirs[:,:,:,0,0]
	outvec[:,:,:,1] = dirs[:,:,:,0,1]
	outvec[:,:,:,2] = dirs[:,:,:,0,2]

	# Save output files
	print('')
	print('... saving output files...')
	print('')

	buffer_header = dwiobj.header
	buffer_header.set_data_dtype('float64')   # Make sure we save quantitative maps as float64, even if input header indicates a different data type

	buffer_string=''
	seq_string = (outbase,'_fa.nii')
	fa_outfile = buffer_string.join(seq_string)
		
	buffer_string=''
	seq_string = (outbase,'_ad.nii')
	ad_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (outbase,'_rd.nii')
	rd_outfile = buffer_string.join(seq_string) 

	buffer_string=''
	seq_string = (outbase,'_md.nii')
	md_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (outbase,'_ak.nii')
	ak_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (outbase,'_rk.nii')
	rk_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (outbase,'_mk.nii')
	mk_outfile = buffer_string.join(seq_string)

	buffer_string=''
	seq_string = (outbase,'_dir.nii')
	dir_outfile = buffer_string.join(seq_string)

	fa_obj = nib.Nifti1Image(FA,dwiobj.affine,buffer_header)
	nib.save(fa_obj, fa_outfile)

	ad_obj = nib.Nifti1Image(AD,dwiobj.affine,buffer_header)
	nib.save(ad_obj, ad_outfile)

	rd_obj = nib.Nifti1Image(RD,dwiobj.affine,buffer_header)
	nib.save(rd_obj, rd_outfile)

	md_obj = nib.Nifti1Image(MD,dwiobj.affine,buffer_header)
	nib.save(md_obj, md_outfile)

	ak_obj = nib.Nifti1Image(AK,dwiobj.affine,buffer_header)
	nib.save(ak_obj, ak_outfile)

	rk_obj = nib.Nifti1Image(RK,dwiobj.affine,buffer_header)
	nib.save(rk_obj, rk_outfile)

	mk_obj = nib.Nifti1Image(MK,dwiobj.affine,buffer_header)
	nib.save(mk_obj, mk_outfile)

	dir_obj = nib.Nifti1Image(outvec,dwiobj.affine,buffer_header)
	nib.save(dir_obj, dir_outfile)


	# Done
	print('')
	print('... done. Enjoy your diffusion maps -- hope they are not too noisy!')
	print('')

