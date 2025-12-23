#!/usr/bin/env python3

import skimage as ski
import sys
import os
import numpy as np
from sys import exit
from tqdm import tqdm
import nrrd
from collections import OrderedDict


# get the number of command-line arguments
nargs = len(sys.argv)

# the user has to provide at least an input and output directory
if(nargs < 3):
    print("Usage: stack2nrrd.py <input_directory> <output_file> <x_spacing> <y_spacing> <z_spacing>")
    exit()
    
input_path = sys.argv[1]
output_file = sys.argv[2]

# store the spacing based on user input (default is isotropic)
spacing = [1.0, 1.0, 1.0]
if(nargs > 3):
    spacing[2] = float(sys.argv[3])
if(nargs > 4):
    spacing[1] = float(sys.argv[4])
if(nargs > 5):
    spacing[0] = float(sys.argv[5])

# load the file list and provide the user with information about the stack
fnames = os.listdir(input_path)



# create an array large enough to store the resulting files
nfiles = len(fnames)
proto_image_path = os.path.join(input_path, fnames[0])
proto_image = ski.io.imread(proto_image_path)

# calculate the stack size
stack_size = [nfiles, proto_image.shape[0], proto_image.shape[1]]
data_type = proto_image.dtype

print(str(nfiles) + " images detected for processing into a " + str(stack_size[0]) + "x" + str(stack_size[1]) + "x" + str(stack_size[2]) + " nrrd file")
print("output file: " + output_file)
print("voxel spacing: " + "x = " + str(spacing[2]) + ", y = " + str(spacing[1]) + ", z = " + str(spacing[0]))
confirm_yn = input("continue [Y/n] ")

if(confirm_yn != 'Y'):
    exit()


# create a volume to fit the stack
S = np.zeros(stack_size, dtype=data_type)

for fi in tqdm(range(nfiles)):
    f = os.path.join(input_path, fnames[fi])
    image = ski.io.imread(f)
    S[fi, :, :] = np.transpose(image)
    
# generate the header file
output_header = OrderedDict()
output_header['dimension'] = 3
output_header['sizes'] = S.shape
output_header['type'] = 'uint8'
output_header['encoding'] = 'gzip'
output_header['spacings'] = spacing

nrrd.write(output_file, S, header=output_header)