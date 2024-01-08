#!/bin/bash

# Script needs to be executed once to download interaction-analyzer

# Download archive
wget https://drive.switch.ch/index.php/s/zCl4uiFDNbtdhOS/download  -O squeezemd-bin.tar.gz

# Unpack archiv in ~/tools/
INSTALLDIR=~/toolz/
mkdir $INSTALLDIR
tar -xvf squeezemd-bin.tar.gz -C $INSTALLDIR

# Save paths in bashrc
echo "# foldX
export PATH=\$PATH:/home/pixelline/tools/foldx/

# Interaction Analyzer of Martin
export PATH=\$PATH:/home/pixelline/tools/interaction-analyzer/" >> ~/.bashrc

# source
source ~/.bashrc
