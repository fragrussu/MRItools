#!/bin/bash
### Tool for motion correction of 4D multi-contrast images, useful in diffusion-relaxation MRI
#
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
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.

isaffine=0
isdebug=0
ln=3
lp=3
inorder=3
changexfm=0
changetmp=0
PARAMS=""
export FSLOUTPUTTYPE=NIFTI
while (( "$#" )); do
  case "$1" in
    -h|--help)
      echo ""
      echo "Usage: moco4d.sh [OPTION] INPUT OUTPUT"
      echo ""
      echo "Correct motion in 4D INPUT nifti file using NiftyReg" 
      echo "linear registration tool (reg_aladin) and store"
      echo "motion correction results to OUTPUT nifti file."
      echo "Motion correction reference is calculated as the"
      echo "mean of all volumes in INPUT"
      echo ""
      echo "Mandatory arguments to long options are mandatory" 
      echo "for short options too."
      echo ""
      echo ""
      echo "  -d, --debug                do not delete temporary"
      echo "                             folder storing"
      echo "                             intermediate results"
      echo ""
      echo "  -a, --aff                  use full affine registration"
      echo "                             (default rigid registration," 
      echo "                              i.e. reg_aladin with option"
      echo "                              -rigOnly)"
      echo ""
      echo "  -t FOLDER, --tmp FOLDER    use FOLDER to store"
      echo "                             intermediate results"
      echo "                             (default OUTPUT.tmp)"
      echo ""
      echo "  -x FOLDER, --xfm FOLDER    use FOLDER to store"
      echo "                             registration" 
      echo "                             transformations"
      echo "                             (default OUTPUT.xfm)"
      echo "                             Registration transformations"
      echo "                             are stored like vol0.txt,"
      echo "                             vol1.txt, ..., volN.txt"
      echo ""
      echo "  --ln NUM                   use reg_aladin with option" 
      echo "                             -ln NUM (default 3)"
      echo ""
      echo "  --lp NUM                   use reg_aladin with -lp NUM"
      echo "                             -lp NUM (default 3)"
      echo ""
      echo "  -i NUM, --interp NUM       use interpolation order NUM" 
      echo "                             (0, 1, 3, 4; default 3;"
      echo "                             (0=NN, 1=LIN; 3=CUB, 4=SINC)"
      echo ""
      echo "  -h, --help                 print this help manual"
      echo ""
      echo "Dependencies:"
      echo " - FSL https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/"
      echo " - NiftyReg http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftyReg"
      echo ""
      echo "Author: Francesco Grussu <francegrussu@gmail.com>"
      echo ""
      exit	
      shift
      ;;
    -a|--aff)
      isaffine=1
      shift
      ;;
    -d|--debug)
      isdebug=1
      shift
      ;;
    -t|--tmp)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        changetmp=1
        newtmp=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -x|--xfm)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        changexfm=1
        newxfm=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    --ln)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        ln=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    --lp)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        lp=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -i|--interp)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        inorder=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -*|--*=) 
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) 
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

# Get positional arguments
eval set -- "$PARAMS"
niftiin=$1
niftiout=$2

# Get output folder names
if [ $changetmp -eq 1 ]
then
	tmpfolder=$newtmp
else
	tmpfolder=$niftiout".tmp"
fi

if [ $changexfm -eq 1 ]
then
	xfmfolder=$newxfm
else
	xfmfolder=$niftiout".xfm"
fi

# Print feedback to the user
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                        runmoco.sh                         "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo ""
echo "** 4D input Nifti file          : "$niftiin
echo "** 4D output Nifti file         : "$niftiout
echo "** output transformation folder : "$xfmfolder
echo "** temporary working folder     : "$tmpfolder
echo "** debug mode                   : "$isdebug
echo "** affine registration          : "$isaffine
echo "** reg_aladin --ln parameter    : "$ln
echo "** reg_aladin --lp parameter    : "$lp
echo "** interpolation order          : "$inorder
echo ""
echo ""
echo ""

# Create temporary and transformation folder
mkdir $xfmfolder
mkdir $tmpfolder

# Get number of 4D volumes
nvol=`fslval $niftiin dim4`
let nvol_minus_1=nvol-1

# Create reference
echo "    ... creating registration target"
echo ""
fslmaths $niftiin -Tmean $tmpfolder"/ref.nii"

# Extract 4D volumes and register them to reference
echo "    ... correcting volumes:"
for vv in `seq 0 $nvol_minus_1`
do
	echo "                           volume "$vv
	fslroi $niftiin $tmpfolder"/vol"$vv".nii" $vv 1
	if [ $isaffine -eq 1 ]
	then
		reg_aladin -ref $tmpfolder"/ref.nii" -flo $tmpfolder"/vol"$vv".nii" -res $tmpfolder"/res_vol"$vv".nii" -ln $ln -lp $lp -aff $xfmfolder"/vol"$vv".txt" -interp $inorder  
	else
		reg_aladin -ref $tmpfolder"/ref.nii" -flo $tmpfolder"/vol"$vv".nii" -res $tmpfolder"/res_vol"$vv".nii" -ln $ln -lp $lp -aff $xfmfolder"/vol"$vv".txt" -interp $inorder -rigOnly
	fi

done 
echo ""

# Merge volumes
echo "    ... merging corrected volumes:"
echo ""
fulllist=""
for vv in `seq 0 $nvol_minus_1`
do
	fulllist=$fulllist" "$tmpfolder"/res_vol"$vv".nii"
done
fslmerge -t $tmpfolder"/res_all.nii" $fulllist


# Store the final output motion-corrected time series and remove temporary folder if required
echo "    ... cleaning up"
cp $tmpfolder"/res_all.nii" $niftiout
if [ $isdebug -eq 0 ]
then
	rm -r -f $tmpfolder
fi
echo ""

# Done
echo "    ... done"
echo ""



