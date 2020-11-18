git clone https://github.com/NYU-DiffusionMRI/mppca_denoise
mv -v mppca_denoise/mpdenoise.py tools
rm -r -f -v mppca_denoise

git clone https://github.com/RafaelNH/gibbs-removal
mv -v gibbs-removal/gibbs_removal.py tools
rm -r -f -v gibbs-removal
