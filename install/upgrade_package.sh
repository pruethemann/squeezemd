#!/bin/bash

curent_path=pwd

# Make sure repo is there or at least with a symbolic link
cd ~/Dropbox/code/squeezemd
pip3 install --upgrade .
# python3 setup.py sdist && 

cd $current_path
