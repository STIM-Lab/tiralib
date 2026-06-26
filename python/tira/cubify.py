import numpy as np
import skimage as ski
import os
import re
from tqdm import tqdm
from pathlib import Path


## @package cubify
#  This package contains functions for converting between image stacks
#  and sets of cubes that can be processed in parallel.

# performs a ceiling division
def ceildiv(a, b):
    return -(a // -b)

def cubestr(xi, yi, zi):
    return str(xi) + "_" + str(yi) + "_" + str(zi)

# build the cubify directory structure
def build_directories(out_directory, cube_dim):
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)
        
    for xi in range(cube_dim[0]):
        for yi in range(cube_dim[1]):
            for zi in range(cube_dim[2]):
                subdir_name = cubestr(xi, yi, zi) + "/"
                if not os.path.exists(out_directory + subdir_name):
                    os.makedirs(out_directory + subdir_name)
        
        
# split up the images into separate stacks
def split_stacks(in_directory, out_directory, cube_size):
    
    if in_directory[-1] != "/":
        in_directory = in_directory + "/"
    if out_directory[-1] != "/":
        out_directory = out_directory + "/"
    
    # get the list of files from the input directory
    files = [f for f in os.listdir(in_directory) if os.path.isfile(os.path.join(in_directory, f))]
    files.sort()
    
    # load the first file to determine the X and Y dimensions of the stack
    F0 = ski.io.imread(in_directory + files[0])
    
    # calculate the total number of cubes in each dimension
    cube_dim = (ceildiv(F0.shape[0], cube_size), ceildiv(F0.shape[1], cube_size), ceildiv(len(files), cube_size))
    
    # build the directory structure
    build_directories(out_directory, cube_dim)
    
    # for each file in the large image stack
    print("Splitting Images Into Sub-Directories...")
    for fi in tqdm(range(len(files))):
        
        zi = int(fi / cube_size)
        # load the file
        F = ski.io.imread(in_directory + files[fi])
        for xi in range(cube_dim[0]):
            xi_min = xi * cube_size
            xi_max = (xi+1) * cube_size
            if xi == cube_dim[0] - 1:
                xi_max = F0.shape[0];
                
            for yi in range(cube_dim[1]):
                yi_min = yi * cube_size
                yi_max = (yi+1) * cube_size
                if yi == cube_dim[1] - 1:
                    yi_max = F0.shape[1];
                    
                P = F[xi_min:xi_max, yi_min:yi_max]
                ski.io.imsave(out_directory + cubestr(xi, yi, zi) + "/" + files[fi], P, check_contrast=False)
    return cube_dim
                
# convert the image stacks to numpy files
def stacks_to_npy(out_directory, cube_dim):
    print("Converting Image Stacks into .npy Arrays")
    cubes = cube_dim[0] * cube_dim[1] * cube_dim[2]
    pbar = tqdm(total=cubes)
    
    for xi in range(cube_dim[0]):
        for yi in range(cube_dim[1]):
            for zi in range(cube_dim[2]):
                
                # get the name of the subdirectory
                stack_name = cubestr(xi, yi, zi)
                
                # load all of the files in the associated subdirectory
                stack_directory = out_directory + stack_name + "/"
                files = [f for f in os.listdir(stack_directory) if os.path.isfile(os.path.join(stack_directory, f))]
                files.sort()
                
                # allocate an array to store all of the files
                F0 = ski.io.imread(stack_directory + files[0])
                Stack = np.zeros((F0.shape[0], F0.shape[1], len(files)), dtype=np.float32)
                
                for fi in range(len(files)):
                    Stack[:, :, fi] = ski.io.imread(stack_directory + files[fi]).astype(np.float32)
                stack_filename = out_directory + stack_name + ".npy"
                np.save(stack_filename, Stack)
                pbar.update(1)
    pbar.close()
        
    
## cubify script splits a large image stack into manageable cubes of size cube_size
def cubify(in_directory, out_directory, cube_size):
    
    if in_directory[-1] != "/":
        in_directory = in_directory + "/"
    if out_directory[-1] != "/":
        out_directory = out_directory + "/"
    
    # build substacks
    cube_dimension = split_stacks(in_directory, out_directory, cube_size)
    
    # convert each stack to a NumPy file
    stacks_to_npy(out_directory, cube_dimension)

def decubify(in_directory, out_directory, cube_size=200):

    in_directory = 'filename'
    out_directory = 'filename'
    os.makedirs(out_directory, exist_ok=True)
    
    print(f"Reading 3D cubes from: {in_directory}")
    print(f"Saving 2D slices to: {out_directory}\n")
    
    # file proccesing
    for filename in os.listdir(in_directory):
        if filename.endswith('.npy'):
            print(f"Decubifying file: {filename}")
            
            # Pull the grid coordinates out of the filename
            numbers = re.findall(r'\d+', filename)
            if len(numbers) >= 3:
                cube_x = int(numbers[0])
                cube_y = int(numbers[1])
                cube_z = int(numbers[2])
            else:
                print(f"Skipped {filename}")
                continue
    
            # load the 3D numpy volume data [X, Y, Z]
            full_file_path = os.path.join(in_directory, filename)
            cube_3d = np.load(full_file_path, allow_pickle=True)
    
            # finding the midpoint along the z axis
            mid_z = cube_3d.shape[2] // 2
            
            # Slicing syntax: [all_x, all_y, exact_z_layer]
            # This collapses the 3D array down into a flat 2D matrix
            tile_2d = cube_3d[:, :, mid_z]
    
            # saving the 2d data and keeping the coordinates
            output_filename = f"cube_2d_X{cube_x}_Y{cube_y}_Z{cube_z}.npy"
            full_output_path = os.path.join(out_directory, output_filename)
            
            np.save(full_output_path, tile_2d)
            print(f"saved 2D slice data to: {output_filename}")

    
