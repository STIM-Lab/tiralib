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
    print("Splitting images into subdirectories...")
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
    print("Converting image stacks into .npy arrays")
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
#
#  @param in_directory is the directory containing volume slices that will be turned into cubes
#  @param out_directory is the destination directory for the cubes
#  @param cube_size is the maximum cube size: all of the "internal" cubes will be this size, but cubes at
#         the edges may be smaller
def cubify(in_directory, out_directory, cube_size):
    
    if in_directory[-1] != "/":
        in_directory = in_directory + "/"
    if out_directory[-1] != "/":
        out_directory = out_directory + "/"
    
    # build substacks
    cube_dimension = split_stacks(in_directory, out_directory, cube_size)
    
    # convert each stack to a NumPy file
    stacks_to_npy(out_directory, cube_dimension)
    
## Retrieve the cube tile indices from a filename
def tilenum(filename):
    # Pull the grid coordinates out of the filename
    numbers = re.findall(r'\d+', filename)
    if len(numbers) < 3:
        raise Exception("A filename (" + filename + ") has an unexpected format: at least three numbers are expected in the filename")
    cube_x = int(numbers[0])
    cube_y = int(numbers[1])
    cube_z = int(numbers[2])
    
    return (cube_x, cube_y, cube_z)
    
## Generate a grid of tiles, where each grid position provides the filename of the associated tile
def namegrid(filenames):
    
    #find the maximum index for each dimension
    max_x = 0
    max_y = 0
    max_z = 0
    
    for f in filenames:
        indices = tilenum(f.name)
        if max_x < indices[0]:
            max_x = indices[0]
        if max_y < indices[1]:
            max_y = indices[1]
        if max_z < indices[2]:
            max_z = indices[2]
    
    NameGrid = np.full((max_x + 1, max_y + 1, max_z + 1), "", dtype=object)
    
    for f in filenames:
        indices = tilenum(f.name)
        NameGrid[indices[0], indices[1], indices[2]] = f
        
    return NameGrid

## Calculates the size of the entire volume from a NameGrid. The main challenge
#  here is that the last tile along any dimension may be smaller than the other tiles inside
#  the grid. This is because the entire volume may not be evenly divisible by the cube size.
#  Since the cube at the starting point of the grid (0, 0, 0) will be the largest, we start
#  there and then update the size based on the size of the last cube along that dimension.
def volshape(namegrid):
    
    # get the number of cubes along each dimension
    num_cubes = namegrid.shape
    
    # load the size of the corner cube at 0,0,0
    filename = namegrid[0, 0, 0]
    CornerCube = np.load(filename)
    corner_shape = CornerCube.shape
    
    # get the size of the slice EXCLUDING the last cube (which may be a different size)
    x_size = corner_shape[0] * (num_cubes[0] - 1)
    y_size = corner_shape[1] * (num_cubes[1] - 1)
    z_size = corner_shape[2] * (num_cubes[2] - 1)
    
    # add the size of the last cube
    LastCube = np.load(namegrid[-1, -1, -1])
    x_size = x_size + LastCube.shape[0]
    y_size = y_size + LastCube.shape[1]
    z_size = z_size + LastCube.shape[2]
    
    return (x_size, y_size, z_size)

## Returns the data type of the cubes in a tiled grid
def dtype(namegrid):
    filename = namegrid[0, 0, 0]
    CornerCube = np.load(filename)
    return CornerCube.dtype

## Returns the number of slices in the current Z plank (based on the number of slices in the
#  corner cube)
def cubeshape(namegrid, zi):
    filename = namegrid[0, 0, 0]
    CornerCube = np.load(filename)
    return CornerCube.shape

## Combines a set of equally-sized cubes (generated using the cubify function) into
#  a set of 2D slices across the entire volume.

def decubify(in_directory, out_directory):

    # create the output directory (if it doesn't already exist)
    os.makedirs(out_directory, exist_ok=True)
    
    # get all of the NPY files in the input directory (we assume these are all cubes)
    npy_files = list(Path(in_directory).glob("*.npy"))

    # create a tiled "filename grid" of all of the NPY files    
    NameGrid = namegrid(npy_files)
    
    # calculate the size of a 2D slice from the FileGrid
    VolumeShape = volshape(NameGrid)
    
    # allocate a 2D slice using the same data type as the NPY cubes
    data_type = dtype(NameGrid)
    Slice = np.zeros((VolumeShape[0], VolumeShape[1]), dtype=data_type)
    
    slice_counter = 0
    print("Saving volume slices")
    pbar = tqdm(total=VolumeShape[2])
    # for each "plank" (layer of cubes) along the Z axis
    for zi in range(NameGrid.shape[2]):
        
        # get the number of slices in the plank
        cube_shape = cubeshape(NameGrid, zi)
        
        # iterate through each slice in the plank
        for si in range(cube_shape[2]):
            
            # iterate through each XY tile in the slice
            for xi in range(NameGrid.shape[0]):
                for yi in range(NameGrid.shape[1]):
                    
                    # load the corresponding cube
                    Cube = np.load(NameGrid[xi, yi, zi])
                    
                    # copy the slice from the cube to the volume slice
                    x_start = xi * cube_shape[0]
                    x_end = x_start + Cube.shape[0]
                    y_start = yi * cube_shape[1]
                    y_end = y_start + Cube.shape[1]
                    Slice[x_start:x_end, y_start:y_end] = Cube[:, :, si]
                    
            # save the entire slice to disk
            slice_filename = out_directory + "/" + str(slice_counter) + ".npy"
            np.save(slice_filename, Slice)
            slice_counter = slice_counter + 1
            pbar.update(1)
    pbar.close()

    
