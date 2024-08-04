# INSTALL
***

# Description

The *squeezeMD* code can be installed using a conda environment. This installs
all required modules and the squeezemd code. Additional dependencies can be easly installed
following this script. I recommand to use mamba or in particular micromamba to install
everything. This allows a clean and fast install

## Micromamba install

1. This install micromamba on Linux 64 in the folder bin # TODO create tools/micromamba_bin
> curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj ~/tools/micromamba_bin \
> curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
> ~/tools/micromamba_bin/micromamba shell init -s bash -p ~/tools/micromamba \
> source ~/.bashrc
2. Execute the following command the root of this github folder
> micromamba create -f squeeze_env.yml
3. Activate the environment after a successful install
> micromamba activate squeeze
4. Install additional bins (foldX, posCo)
chmod +x install_bins.sh
./install_bins.sh

Version 1.5.8
1. "${SHELL}" <(curl -L micro.mamba.pm/install.sh)



## Test install
1. Test Cuda and OpenMM install
> python3 -m openmm.testInstallation
if it fails it's probably a cuda dependency problem:
Downgrad cuda:
> mamba install -c conda-forge cudatoolkit=11.4

if libtinfo.so.5 is missing
> sudo apt-get install libtinfo5



