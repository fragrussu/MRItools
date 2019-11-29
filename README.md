# MRItools
MRItools is a collection of command line utilities written in Python 3 that some may find useful to conduct MRI research. MRItools may be useful to you if you were planning to:
* fit ADC or DKI signal representations to diffusion MRI data;
* convert real/imaginary MR images to magnitude/phase or vice versa;
* rephase complex-valued MR images to get rid of the Rician noise floor and work with real-valued data and Gaussian noise;
* convert diffusion gradient directions from scanner to image space;
* perform repeatability analysis (coming soon).

# Dependencies
To use MRItools you need a Python 3 distribution such as [Anaconda](http://www.anaconda.com/distribution). Additionally, you need the following third party modules/packages:
* [DiPy](http://dipy.org/)
* [NumPy](http://numpy.org)
* [Nibabel](http://nipy.org/nibabel)


# Installation (Linux and MacOS)
Gettins MRItools is extremely easy.

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

# Description of tools
The following command line tools are available.
* `getADCDKI.py`: to fit a mono-dimensional diffusion MRI decay (i.e. b-value dependence only, no directional dependence) signal representation (choose between apparent diffusion coefficient or diffusion kurtosis);
* `fitdki.py`: to fit the [Diffusion Kurtosis Tensor](http://doi.org/10.1002/mrm.20508) signal representation to multi-shell diffusion MRI data in NIFTI format. It is essentially a command-line wrapper of [this tutorial](http://dipy.org/documentation/1.0.0./examples_built/reconst_dki) made available within the [DiPy](http://dipy.org/) project;
* `getphase.py`: to calculate MR phase from two NIFTIs storing the real and imaginary signals;
* `getPM.py`: to calculate MR phase and magnitude from two NIFTIs storing the real and imaginary signals;
* `getRI.py`: to calculate real and imaginary signals from two NIFTIs storing MR magnitude and phase.



You can run the tools from command line, for instance using a Bash or C shell. Importantly, each tool has a manual: to print it, simply type in your terminal
```
python </PATH/TO/TOOL> --help
```
(for example, `python ./MRItools/tools/getADCDKI.py --help`).

# If you use MRItools
If you use MRItools in your research, please remember to cite our paper: (just wait 48 hours to find out!)

# License
MRItools is distributed under the BSD 2-Clause License, Copyright (c) 2019, University College London. All rights reserved.
Link to license [here](http://github.com/fragrussu/MRItools/blob/master/LICENSE).

# Acknowledgements
Funding from the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 634541) and from the United Kingdom Engineering and Physical Sciences Research Council (EPSRC R006032/1 and M020533/1) is acknowledged.


