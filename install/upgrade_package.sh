#!/bin/bash

curent_path=pwd

# Make sure repo is there or at least with a symbolic link
cd ~/Dropbox/code/SqueezeMD/squeezemd
python3 setup.py sdist && pip3 install --upgrade .

cd $current_path
