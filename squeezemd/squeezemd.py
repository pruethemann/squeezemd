"""Main module."""

def add(a,b):
    return a + b

def subtract(a,b):
    return a - b

import math

def circle_area(r):
    return math.pi * r**2

import subprocess
import os

def call_my_binary(arg1, arg2, arg3):
    binary_path = '/home/pixelline/ownCloud/Institution/code/SqueezeMD/squeezemd/bin/7_interaction-analyzer.x'
    result = subprocess.run([binary_path, str(arg1), str(arg2), str(arg3)], capture_output=True)
    return result.stdout

