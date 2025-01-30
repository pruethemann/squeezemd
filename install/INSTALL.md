# INSTALL
***

# Description

The *squeezeMD* code can be installed using a conda environment. This installs
all required modules and the squeezemd code. Additional dependencies can be easly installed
following this script. I recommand to use mamba or in particular micromamba to install
everything. This allows a clean and fast install

## Linux

1. Install micromamba
> "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
2. Install the conda environment
> micromamba create -f environment.yml  \
> micromamba activate squeeze-dev
3. Install additional binaries (Posco and foldX)
chmod +x install_bins_linux.sh \
./install_bins_linux.sh
4. source .bashrc
5. Upgrade the current squeezeMD code
cd squeezemd  \
python3 setup.py sdist && pip3 install --upgrade . 
5. Test binaries by executing:
> po-sco
> foldx_20251231


## Test install
1. Test Cuda and OpenMM install
> python3 -m openmm.testInstallation

You should expect sth like this:
>There are 3 Platforms available:

>1 Reference - Successfully computed forces
>2 CPU - Successfully computed forces
>3 CUDA - Successfully computed forces

if it fails it's probably a cuda dependency problem:
Downgrade cuda:
> mamba install -c conda-forge cudatoolkit=11.4

if libtinfo.so.5 is missing
> sudo apt-get install libtinfo5

2. Run the demo workflow in the folder demo
> cd demo  \
> micromamba activate squeeze  \
> squeeze -j4 -n  \
> squeeze -j4

Make sure algoriths completes successfully and uses GPU during Molecular dynamics step



