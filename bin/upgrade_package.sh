#!/bin/bash

curent_path=pwd


cd /home/pixelline/ownCloud/Institution/code/SqueezeMD/squeezemd
python3 setup.py sdist && pip3 install --upgrade .

cd $current_path
