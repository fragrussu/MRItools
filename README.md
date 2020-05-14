# MRItools
MRItools is a collection of command line utilities written in Python 3 that some may find useful to conduct MRI research. MRItools may be useful to you if you were planning to:
* fit ADC or DKI signal representations to diffusion MRI data;
* convert real/imaginary MR images to magnitude/phase or vice versa;
* convert diffusion gradient directions from scanner to image space;
* perform repeatability analyses on quantitative MRI metrics;
* rephase complex-valued MR images to get rid of the Rician noise floor and work with real-valued data and Gaussian noise (coming soon).

# Dependencies
To use MRItools you need a Python 3 distribution such as [Anaconda](http://www.anaconda.com/distribution). Additionally, you need the following third party modules/packages:
* [DiPy](http://dipy.org/)
* [NumPy](http://numpy.org)
* [Nibabel](http://nipy.org/nibabel)
* [SciPy](http://www.scipy.org)
* [Nipype](http://nipype.readthedocs.io/en/latest)


# Download 
Getting MRItools is extremely easy: cloning this repository is all you need to do. The tools would be ready for you to run.

If you use Linux or MacOS:

1. Open a terminal;
2. Navigate to your destination folder;
3. Clone MRItools:
```
git clone https://github.com/fragrussu/MRItools.git 
```
4. MRItools is ready for you in `./MRItools` with the tools available here: 
```
./MRItools/tools
```
5. You should now be able to use the code. Try to print the manual of a script, for instance of `fitdki.py`, to make sure this is really the case:
```
python ./MRItools/tools/fitdki.py --help
```

# MRItools available
The following command line tools are available.
* [`fitdki.py`](http://github.com/fragrussu/MRItools/blob/master/tools/fitdki.py): to fit the [Diffusion Kurtosis Tensor](http://doi.org/10.1002/mrm.20508) signal representation to multi-shell diffusion MRI data in NIFTI format. It is essentially a command-line wrapper of [this tutorial](http://dipy.org/documentation/1.0.0./examples_built/reconst_dki) made available within the [DiPy](http://dipy.org/) project;
* [`getADCDKI.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getADCDKI.py): to fit a mono-dimensional diffusion MRI decay (i.e. b-value dependence only, no directional dependence) signal representation (choose between apparent diffusion coefficient or diffusion kurtosis);
* [`getphase.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getphase.py): to calculate MR phase from two NIFTIs storing the real and imaginary signals;
* [`getPM.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getPM.py): to calculate MR phase and magnitude from two NIFTIs storing the real and imaginary signals;
* [`getRI.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getRI.py): to calculate real and imaginary signals from two NIFTIs storing MR magnitude and phase;
* [`rotgrad2imgspace.py`](http://github.com/fragrussu/MRItools/blob/master/tools/rotgrad2imgspace.py): to rotate diffusion MRI gradient directions from scanner space to image space. The code replicates the conversion of gradient directions from scanner to image image space performed by the freely available [`dcm2niix`](http://github.com/rordenlab/dcm2niix) tool;
* [`iccvoxelwise.py`](http://github.com/fragrussu/MRItools/blob/master/tools/iccvoxelwise.py): to evaluate voxel-wise intraclass correlation coefficient (ICC) maps according to the formalism of [Shrout and Fleiss, Psychological Bulletin 1979](http://doi.org/10.1037/0033-2909.86.2.420) (indices ICC31 and ICC11; for ICC11, estimates of variances are obtained "by hand" according to Appendix B of [Grussu et al, NeuroImage 2015](http://doi.org/10.1016/j.neuroimage.2015.01.045)). 



You can run the tools from command line, for instance using a Bash or C shell. Importantly, each tool has a manual: to print it, simply type in your terminal
```
python </PATH/TO/TOOL> --help
```
(for example, `python ./MRItools/tools/getADCDKI.py --help`).

# If you use MRItools
If you use MRItools in your research, please remember to cite our paper:

"Multi-parametric quantitative in vivo spinal cord MRI with unified signal readout and image denoising". Grussu F, Battiston M, Veraart J, Schneider T, Cohen-Adad J, Shepherd TM, Alexander DC, Fieremans E, Novikov DS, Gandini Wheeler-Kingshott CAM; [NeuroImage 2020, 217: 116884](http://doi.org/10.1016/j.neuroimage.2020.116884) (DOI: 10.1016/j.neuroimage.2020.116884).



If you use [`iccvoxelwise.py`](http://github.com/fragrussu/MRItools/blob/master/tools/iccvoxelwise.py), please also cite:

"Neurite orientation dispersion and density imaging of the healthy cervical spinal cord in vivo". Grussu F, Schneider T, Zhang H, Alexander DC, Wheeler–Kingshott CAM; [NeuroImage 2015, 111: 590-601](http://doi.org/10.1016/j.neuroimage.2015.01.045) (DOI: 10.1016/j.neuroimage.2015.01.045).



# License
MRItools is distributed under the BSD 2-Clause License, Copyright (c) 2019, 2020 University College London. All rights reserved.
Link to license [here](http://github.com/fragrussu/MRItools/blob/master/LICENSE).

# Acknowledgements
Funding from the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 634541) and from the United Kingdom Engineering and Physical Sciences Research Council (EPSRC R006032/1 and M020533/1) is acknowledged.


