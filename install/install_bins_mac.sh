#!/bin/zsh

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

# Save paths in zshrc (extended bash on current Mac)
echo "# foldX
export PATH=\$PATH:~/tools/foldx/foldxMacC11_0" >> ~/.zshrc

# PosCo
echo "# foldX
export PATH=\$PATH:~/tools/interaction-analyzer/mac" >> ~/.zshrc

# source
source ~/.zshrc

rm squeezemd-bin.tar.gz
