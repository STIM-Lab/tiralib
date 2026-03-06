import numpy as np
import skimage as ski
import os

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
def split_stacks(in_directory, out_directory, cube_dim):
    # for each file in the large image stack
    for fi in range(len(files)):
        
        zi = int(fi / S)
        # load the file
        F = ski.io.imread(in_directory + files[fi])
        for xi in range(cube_dim[0]):
            xi_min = xi * 200
            xi_max = (xi+1) * 200
            if xi == cube_dim[0] - 1:
                xi_max = F0.shape[0];
                
            for yi in range(cube_dim[1]):
                yi_min = yi * 200
                yi_max = (yi+1) * 200
                if yi == cube_dim[0] - 1:
                    yi_max = F0.shape[0];
                    
                P = F[xi_min:xi_max, yi_min:yi_max]
                ski.io.imsave(out_directory + cubestr(xi, yi, zi) + "/" + files[fi], P)
                
# convert the image stacks to numpy files
def stacks_to_npy(out_directory, cube_dim):
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
                Stack = np.zeros((F0.shape[0], F0.shape[1], len(files)))
                
                for fi in range(len(files)):
                    Stack[:, :, fi] = ski.io.imread(stack_directory + files[fi])
                stack_filename = out_directory + stack_name + ".npy"
                np.save(stack_filename, Stack)
        
    

in_directory = "/home/david/tissue/kesm_small_0/images/"
out_directory = "/home/david/test/"

S = 200

# get the list of files from the input directory
files = [f for f in os.listdir(in_directory) if os.path.isfile(os.path.join(in_directory, f))]
files.sort()

# load the first file to determine the X and Y dimensions of the stack
F0 = ski.io.imread(in_directory + files[0])

# calculate the total number of cubes in each dimension
cube_dim = (ceildiv(F0.shape[0], S), ceildiv(F0.shape[1], S), ceildiv(len(files), S))

# build the directory structure
build_directories(out_directory, cube_dim)

# build substacks
split_stacks(in_directory, out_directory, cube_dim)

# convert each stack to a NumPy file
stacks_to_npy(out_directory, cube_dim)




