#!/bin/bash

# Script needs to be executed once to install posco and foldX

# Download archive containing all bins
echo "Download software in current folder"
curl -L -o squeezemd-bin.tar.gz "https://www.dropbox.com/scl/fi/vrsujfeh4ka11xdfgt6b4/squeezemd-bin.tar.gz?rlkey=ibj0f6pxez6hcfooowqwvzdsr&st=ff0w04s5&dl=1"

# Unpack archiv in ~/tools/
INSTALLDIR=~/tools/
mkdir -p $INSTALLDIR
tar -xvf squeezemd-bin.tar.gz -C $INSTALLDIR

# Save paths in bashrc (extended bash on current Mac)
echo "# foldX
export PATH=\$PATH:~/tools/foldX/foldx_Linux" >> ~/.bashrc

# PosCo
echo "# Po-Sco
export PATH=\$PATH:~/tools/po-sco" >> ~/.bashrc

# source
source ~/.bashrc

rm squeezemd-bin.tar.gz
