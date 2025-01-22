#!/bin/bash

# Script needs to be executed once to install bins for MAC:
# - poSco
# - foldX

# Download archive containing all bins
echo "Download software in current folder"
curl -L -o squeezemd-bin.tar.gz "https://www.dropbox.com/scl/fi/ngo4nhh0mcj09goob10wk/squeezemd-bin.tar.gz?rlkey=1kw69tmcrrgmqzbx096ksd7l1&st=ouhq22rf&dl=1"

# Unpack archiv in ~/tools/
INSTALLDIR=~/tools/
mkdir -p $INSTALLDIR
tar -xvf squeezemd-bin.tar.gz -C $INSTALLDIR

# Save paths in bashrc (extended bash on current Mac)
echo "# foldX
export PATH=\$PATH:~/tools/foldX/foldx5_1Linux64_0" >> ~/.bashrc

# PosCo
echo "# Po-Sco
export PATH=\$PATH:~/tools/po-sco" >> ~/.bashrc

# source
source ~/.bashrc

rm squeezemd-bin.tar.gz
