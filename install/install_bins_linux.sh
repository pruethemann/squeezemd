#!/bin/bash

# Script needs to be executed once to install bins for MAC:
# - poSco
# - foldX

# Download archive containing all bins
echo "Download software in current folder"
curl -L -o squeezemd-bin.tar.gz "https://www.dropbox.com/scl/fi/4v5f8z42m38yypj48biwb/squeezemd-bin.tar.gz?rlkey=04m2g789pry8pqiribd7u0uii&dl=1"

# Unpack archiv in ~/tools/
INSTALLDIR=~/tools/
mkdir -p $INSTALLDIR
tar -xvf squeezemd-bin.tar.gz -C $INSTALLDIR

# Save paths in bashrc (extended bash on current Mac)
echo "# foldX
export PATH=\$PATH:~/tools/foldx/linux64" >> ~/.bashrc

# PosCo
echo "# foldX
export PATH=\$PATH:~/tools/interaction-analyzer/linux" >> ~/.bashrc

# source
source ~/.bashrc

rm squeezemd-bin.tar.gz
