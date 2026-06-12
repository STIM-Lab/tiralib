import numpy as np
import skimage as ski
import os
from tqdm import tqdm
from pathlib import Path


## @package cubify
#  This package contains functions for converting a large image stack
#  into a set of cubes that can be processed in parallel.

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

    
def otsu_batch(input_npy_dir, output_mask_dir):

    input_npy_dir = Path(input_npy_dir)
    output_mask_dir = Path(output_mask_dir)

    # make output folder
    output_mask_dir.mkdir(parents=True, exist_ok=True)
    
    sorted_npy = sorted(input_npy_dir.glob("*.npy"))

    # process all npy cubes
    for ni in tqdm(range(len(sorted(input_npy_dir.glob("*.npy"))))):
        
        npy_file = sorted_npy[ni]

        #print(f"processing : {npy_file.name}")

        # load cube
        vol = np.load(npy_file).astype(np.float32)

        # compute global 3d otsu threshold
        thresh = ski.filters.threshold_otsu(vol)

        #print(f"otsu threshold : {thresh}")
        
        inside = vol <= thresh

        # binary segmentation
        mask = np.ones(vol.shape, dtype=np.float32) * 255.0
        mask[inside] = 0.0

        # convert boolean to uint8
        #mask = mask.astype(np.float32) * 255.0

        # save binary cube
        out_file = output_mask_dir / f"{npy_file.stem}_otsu3d.npy"
        np.save(out_file, mask)

## This function executes the RSF GPU code for each volume in a specified directory.
#  The volumes are assumed consist of (1) a raw image volume and (2) an initial
#  binary segmentation.
#
#  @param exec_dir is the directory of the rsf_gpu executable
#  @param bin_dir is the directory containing the initial binary segmentation
#  @param raw_dir is the directory containing the raw image volumes corresponding to the binary segmentations
#  @param sigma_l is the standard deviation of the blur kernel used in the localization function
#  @param sigma_k is the standard deviation of the blur kernel used in the fitting function
#  @param T is the number of time steps to simulate
#  @param dt is the size of each time step
#  @param wr is the weight for the regularization term
#  @param ws is the weight for the smoothing term
#  @param wf is the weight for the fitting term
#  @param cuda is a true/false value defining whether or not the GPU is used
#  @param dp is a triple defining the size of the volume along each spatial dimension
def rsf_batch(exec_dir, bin_dir, raw_dir, sigma_l=3.0, sigma_k=3.0, T=50, dt=0.1, dp=(1.0, 1.0, 1.0), wr=0.1, ws=0.1, wf=0.1, cuda=True):

    exec_dir = Path(exec_dir)
    bin_dir = Path(bin_dir)
    raw_dir = Path(raw_dir)

    # create output directory
    output_dir = raw_dir.parent / "rsf_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # locate rsf executable
    exe_path = exec_dir / "rsf_gpu.exe"

    if not exe_path.exists():
        raise FileNotFoundError("could not find rsf_gpu executable")

    # get all raw volumes
    raw_files = sorted(raw_dir.glob("*.npy"))

    print("Running RSF Batch Processing...")
    pbar = tqdm(total=len(raw_files))

    for ri in range(len(raw_files)):

        raw_file = raw_files[ri]

        # locate corresponding binary volume
        bin_file = bin_dir / f"{raw_file.stem}_otsu3d.npy"

        if not bin_file.exists():
            print("missing binary file for " + raw_file.name)
            pbar.update(1)
            continue

        # generate output filename
        out_file = output_dir / f"{raw_file.stem}_rsf.npy"

        # select cpu or gpu execution
        if cuda:
            cuda_device = 0
        else:
            cuda_device = -1

        # build command line arguments
        cmd = [
            str(exe_path),
            "--binary", str(bin_file),
            "--image", str(raw_file),
            "--output", str(out_file),
            "--sigma", str(sigma_l),
            "--sigmaf", str(sigma_k),
            "--t", str(T),
            "--dt", str(dt),
            "--wr", str(wr),
            "--ws", str(ws),
            "--wf", str(wf),
            "--dx", str(dp[0]),
            "--dy", str(dp[1]),
            "--dz", str(dp[2]),
            "--cuda", str(cuda_device),
            
        ]

        # execute rsf
        result = subprocess.run(
                 cmd,
                 cwd=str(exec_dir),
                 capture_output=True,
                 text=True
         )
        
        if result.returncode != 0:
           print("rsf failed on " + raw_file.name)
           print(result.stdout)
           print(result.stderr)


        pbar.update(1)

    pbar.close()




