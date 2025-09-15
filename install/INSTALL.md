# INSTALL

## Description

The *squeezeMD* workflow can be installed on **Linux (Ubuntu recommended)** using a conda environment.
This installs all required Python modules, *squeezeMD*, and additional third-party tools via the provided script.

---

## Requirements

1. Install the following system dependencies first, if not already installed:

```sh
sudo apt update
sudo apt install -y git wget
```

2. Install nvidia drivers [NVIDIA drivers (latest recommended)](https://documentation.ubuntu.com/server/how-to/graphics/install-nvidia-drivers/)

---

## Install Miniconda

*(Skip if conda is already installed)*
[Official instructions](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2)

```sh
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

---

## Install squeezeMD and third-party tools

```sh
# Clone repository
git clone https://github.com/pruethemann/squeezemd.git
cd squeezemd

# Create and activate conda environment
conda env create -f install/environment.yml
conda activate squeeze

# Install FoldX (v5.1) and PoSCo (v1.24)
chmod +x install/install_bins_linux.sh
./install/install_bins_linux.sh
```

Verify installation:

```sh
foldx_20251231 --version
po-sco --version
```

---

## Upgrade squeezeMD

If a newer development version is available on GitHub. Execute this command in the squeezemd git folder.

```sh
pip install --upgrade .
```

---

## Test Installation

```sh
python3 -m openmm.testInstallation
```

✅ Make sure the run completes successfully on GPU.

---

## Troubleshooting

* Ensure **cudatoolkit** matches your NVIDIA driver. Pinning cudatoolkit may be required.
  [CUDA compatibility guide](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)
* Installation currently supports **Linux/Ubuntu only**.

---

⚠️ **Known limitations**

* Only tested on **Ubuntu Linux** (no support for Windows/macOS).
* Requires a **recent NVIDIA driver**; older drivers may not work.
* No full end-to-end validation script provided (only binary version checks + OpenMM test).
