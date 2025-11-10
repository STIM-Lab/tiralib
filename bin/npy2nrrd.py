#!/usr/bin/env python3

import sys
import numpy as np
import nrrd
from collections import OrderedDict

# get the number of command-line arguments
nargs = len(sys.argv)

# the user has to provide at least an input and output directory
if(nargs < 2):
    print("Usage: stack2nrrd.py <input_directory> <output_file> <x_spacing> <y_spacing> <z_spacing>")
    exit()
    
input_file = sys.argv[1]


# store the spacing based on user input (default is isotropic)
spacing = [1.0, 1.0, 1.0]
if(nargs > 3):
    spacing[2] = float(sys.argv[3])
if(nargs > 4):
    spacing[1] = float(sys.argv[4])
if(nargs > 5):
    spacing[0] = float(sys.argv[5])
    


if(len(sys.argv) == 2):
	output_file = input_file + ".nrrd"
else:
	output_file = sys.argv[2]


S = np.load(input_file)

# generate the header file
output_header = OrderedDict()
output_header['dimension'] = 3
output_header['sizes'] = S.shape
output_header['type'] = 'uint8'
output_header['encoding'] = 'gzip'
output_header['spacings'] = spacing

nrrd.write(output_file, S, header=output_header)