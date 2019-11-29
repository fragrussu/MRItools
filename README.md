# MRItools
MRItools is a collection of command line utilities written in Python 3 that some may find useful to process MRI images.

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
* `fitdki.py`: to fit the [Diffusion Kurtosis Imaging](http://doi.org/10.1002/mrm.20508) (DKI) signal representation to multi-shell diffusion MRI data in NIFTI format. It is essentially a command-line wrapper of [this tutorial](http://dipy.org/documentation/1.0.0./examples_built/reconst_dki) made available within the [DiPy](http://dipy.org/) project. 



Each tool has a manual. To print it, simply type in your terminal
```
python </PATH/TO/TOOL> --help
```
(for example, `python ./MRItools/tools/fitdki.py --help`).

# If you use MRItools
If you use MRItools in your research, please remember to cite our paper: (just wait 48 hours to find out!)

# License
MRItools is distributed under the BSD 2-Clause License, Copyright (c) 2019, University College London. All rights reserved.
Link to license [here](http://github.com/fragrussu/MRItools/blob/master/LICENSE).

# Acknowledgements
Funding from the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 634541) and from the United Kingdom Engineering and Physical Sciences Research Council (EPSRC R006032/1 and M020533/1) is acknowledged.


