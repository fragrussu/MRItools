# MRItools
MRItools is a collection of utilities written almost entirely in Python 3 (and in some cases in Matlab) that you may find useful to conduct MRI research. They could help you do things like:
* fit ADC or DKI signal representations to diffusion MRI data;
* convert real/imaginary MR images to magnitude/phase or vice versa;
* convert diffusion gradient directions from scanner to image space;
* perform repeatability analyses on quantitative MRI metrics;
* rephase complex-valued MR images to get rid of the noise floor and work with real-valued data and Gaussian noise;
* perform noise floor mitigation following MP-PCA MR image denoising.

# Dependencies
To use tools written in python you need a Python 3 distribution such as [Anaconda](http://www.anaconda.com/distribution). Additionally, you need the following third party modules/packages:
* [DiPy](http://dipy.org)
* [NumPy](http://numpy.org)
* [Nibabel](http://nipy.org/nibabel)
* [SciPy](http://www.scipy.org)
* [Nipype](http://nipype.readthedocs.io/en/latest)
* [Scikit-image](http://scikit-image.org)
* [tqdm](https://github.com/tqdm/tqdm)
* [joblib](https://joblib.readthedocs.io)
* [MP-PCA](https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py)

To use tools written in Matlab, you need the following toolboxes/functions:
* [MP-PCA](https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/MPdenoising.m).

# Download 
Getting MRItools is extremely easy: cloning this repository is all you need to do. The tools would be ready for you to run.

If you use Linux or MacOS:

1. Open a terminal;
2. Navigate to your destination folder;
3. Clone MRItools:
```
git clone https://github.com/fragrussu/MRItools.git 
```
4. The tools written is python are now ready for you here: 
```
./MRItools/tools
```
You should now be able to use them. Try to print the manual of a script, for instance of `fitdki.py`, to make sure this is really the case:
```
python ./MRItools/tools/fitdki.py --help
```
5. The tools written in Matlab are instead ready for you here:
```
./MRItools/matlabtools
```

# MRItools: python
The following tools written in python are available:
* [`fitdki.py`](http://github.com/fragrussu/MRItools/blob/master/tools/fitdki.py): to fit the [Diffusion Kurtosis Tensor](http://doi.org/10.1002/mrm.20508) signal representation to multi-shell diffusion MRI data in NIFTI format. It is essentially a command-line wrapper of [this tutorial](http://dipy.org/documentation/1.0.0./examples_built/reconst_dki) made available within the [DiPy](http://dipy.org/) project;
* [`getADCDKI.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getADCDKI.py): to fit a mono-dimensional diffusion MRI decay (i.e. b-value dependence only, no directional dependence) signal representation (choose between apparent diffusion coefficient or diffusion kurtosis);
* [`getphase.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getphase.py): to calculate MR phase from two NIFTIs storing the real and imaginary signals;
* [`getPM.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getPM.py): to calculate MR phase and magnitude from two NIFTIs storing the real and imaginary signals;
* [`getRI.py`](http://github.com/fragrussu/MRItools/blob/master/tools/getRI.py): to calculate real and imaginary signals from two NIFTIs storing MR magnitude and phase;
* [`rotgrad2imgspace.py`](http://github.com/fragrussu/MRItools/blob/master/tools/rotgrad2imgspace.py): to rotate diffusion MRI gradient directions from scanner space to image space. The code replicates the conversion of gradient directions from scanner to image image space performed by the freely available [`dcm2niix`](http://github.com/rordenlab/dcm2niix) tool;
* [`iccvoxelwise.py`](http://github.com/fragrussu/MRItools/blob/master/tools/iccvoxelwise.py): to evaluate voxel-wise intraclass correlation coefficient (ICC) maps according to the formalism of [Shrout and Fleiss, Psychological Bulletin 1979](http://doi.org/10.1037/0033-2909.86.2.420) (indices ICC31 and ICC11; for ICC11, estimates of variances are obtained "by hand" according to Appendix B of [Grussu et al, NeuroImage 2015, 111: 590-601](http://doi.org/10.1016/j.neuroimage.2015.01.045));
* [`imgrephaseDC.py`](http://github.com/fragrussu/MRItools/blob/master/tools/imgrephaseDC.py): to rephase complex-valued multi-slice EPI scans data according to the *decorrelated phase filtering method* of [Sprenger et al, Magnetic Resonance in Medicine 2017, 77:559-570](http://doi.org/10.1002/mrm.26138), obtaining real-valued images corrupted by Gaussian noise;
* [`imgrephaseTV.py`](http://github.com/fragrussu/MRItools/blob/master/tools/imgrephaseTV.py): to rephase complex-valued multi-slice EPI scans according to the *total variation real-valued imaging method* of [Eichner et al, NeuroImage 2015, 122:373-384](http://doi.org/10.1016/j.neuroimage.2015.07.074), obtaining real-valued images corrupted by Gaussian noise;
* [`runMPPCA.py`](https://github.com/fragrussu/MRItools/blob/master/tools/runMPPCA.py): to perform [MP-PCA](http:/doi.org/10.1016/j.neuroimage.2016.08.016) denoising (Veraart J et al, NeuroImage 2016, 142: 394-406) and subsequently mitigate the residual noise floor with a custom implementation of the [*method of moments*](http://doi.org/10.1016/j.jmr.2006.01.016) (Koay and Basser, J Magn Reson 2006, 179(2):317-22) - it requires a freely available MP-PCA Python [implementation](https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py);
* [`unring.py`](https://github.com/fragrussu/MRItools/blob/master/tools/unring.py): to perform [Gibbs unringing](https://doi.org/10.1002/mrm.26054) (Kellner et al. Gibbs‐ringing artifact removal based on local subvoxel‐shifts (2016), Magn Reson Med 76:1574–1581) with a freely available Python [implementation](https://github.com/RafaelNH/gibbs-removal/blob/master/gibbs_removal.py).

You can run the tools from command line, for instance using a Bash or C shell. Importantly, each tool has a manual: to print it, simply type in your terminal
```
python </PATH/TO/TOOL> --help
```
(for example, `python ./MRItools/tools/getADCDKI.py --help`).



# MRItools: Matlab
The following tools written in Matlab are available:
* [`MPio_moments.m`](https://github.com/fragrussu/MRItools/blob/master/matlabtools/MPio_moments.m): to perform [MP-PCA](http:/doi.org/10.1016/j.neuroimage.2016.08.016) denoising (Veraart J et al, NeuroImage 2016, 142: 394-406) and subsequently mitigate the residual noise floor with a custom implementation of the [*method of moments*](http://doi.org/10.1016/j.jmr.2006.01.016) (Koay and Basser, J Magn Reson 2006, 179(2):317-22) - it requires a freely available MP-PCA Matlab [implementation](https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/MPdenoising.m). `MPio_moments.m` is essentially equivalent to [`runMPPCA.py`](https://github.com/fragrussu/MRItools/blob/master/tools/runMPPCA.py).
* [`ReadTrajFilesCamino.m`](https://github.com/fragrussu/MRItools/blob/master/matlabtools/ReadTrajFilesCamino.m): to read into Matlab random walker trajectories stored as binary files generated by the [`datasynth`](http://camino.cs.ucl.ac.uk/index.php?n=Man.Datasynth) command in [Camino](http://camino.cs.ucl.ac.uk/index.php) (free open-source software for diffusion MRI data analysis and simulations).


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


