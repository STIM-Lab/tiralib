import numpy as np
import skimage as ski
from pathlib import Path
from tqdm import tqdm
from skimage import measure
import os
import subprocess
import trimesh
import re
import matplotlib.pyplot as plt

# This function sums all tiff images in a directory 
def sum(directory):
    sum_imgs = 0  
    files = sorted (os.listdir(directory))
    for f in files:
        fpath = os.path.join(directory, f)
        img = ski.io.imread(fpath)
        sum_imgs = sum_imgs + img
    return sum_imgs

# This function computes the maximum intensity projection of all tiff images in a directory
def mip(directory):
    files = sorted (os.listdir(directory))
    images = []
       
    for f in files:
        img = ski.io.imread(os.path.join(directory, f))
        images.append(img)
            
    mip = np.max(images, axis = 0)
    return mip


def otsu(input_npy_dir, output_mask_dir):

    input_path = Path(input_npy_dir)
    output_path = Path(output_mask_dir)

    # make output folder
    output_path.mkdir(parents=True, exist_ok=True)
    
    sorted_npy = sorted(input_path.glob("*.npy"))
    
    print("Applying Otsu's method to " + str(len(sorted_npy)) + " .npy arrays")

    # process all npy cubes
    for ni in tqdm(range(len(sorted(input_path.glob("*.npy"))))):
        
        npy_file = sorted_npy[ni]
#from skimage import measure
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
        out_file = output_path / f"{npy_file.stem}_otsu3d.npy"
        np.save(out_file, mask)
        
## This function executes the RSF GPU code for each volume in a specified directory.
#  The volumes are assumed consist of (1) a raw image volume and (2) an initial
#  binary segmentation.
#
#  @param exec_dir is the directory of the rsf_gpu executable
#  @param bin_dir is the directory containing the initial binary segmentation
#  @param raw_dir is the directory contanpy filesining the raw image volumes corresponding to the binary segmentations
#  @param sigma_l is the standard deviation of the blur kernel used in the localization function
#  @param sigma_k is the standard deviation of the blur kernel used in the fitting function
#  @param T is the number of time steps to simulate
#  @param dt is the size of each time step
#  @param wr is the weight for the regularization term
#  @param ws is the weight for the smoothing term
#  @param wf is the weight for the fitting termin_directory, out_directory
#  @param cuda is a true/false value defining whether or not the GPU is used
#  @param dp is a triple defining the size of the volume along each spatial dimension
def rsf(exec_dir, bin_dir, raw_dir, sigma_l=3.0, sigma_k=3.0, T=50, dt=0.1, dp=(1.0, 1.0, 1.0), wr=0.1, ws=0.1, wf=0.1, cuda=True):

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
    
def npy2obj(input_npy_dir, output_dir, level):
    input_npy_dir = os.path.join('foldername')

    filenames = os.listdir(input_npy_dir)
    for fid in range(len(filenames)):
        filename = filenames[fid]
        if filename.endswith('.npy'):
            print(f"processing: {filename}...")
            full_file_path = os.path.join(input_npy_dir, filename)
            sdf_data = np.load(full_file_path)


            verts, faces, normals, values, = measure.marching_cubes(sdf_data, level=0.0)

            output_file = output_dir + filename.replace('.npy', '_surface.obj')

            mesh = trimesh.Trimesh(vertices=verts, faces=faces)
            mesh.export(output_file)
            
def cmap2d(in_directory, out_directory):
   
    in_directory = ('foldername')
    out_directory = ('foldername')
    os.makedirs(out_directory, exist_ok=True)

    CUBE_SIZE = 200

    print(f"Reading slices from: {in_directory}\n")
    print(f"saving images to: {out_directory}\n")

    def check_file_exists(path):
        try:
            os.stat(path)
        except (OSError, ValueError):
            return False
        return True 
    file_count = 0


    for filename in os.listdir(in_directory):   
        if filename.endswith('.npy'):
            print(f"processing files : {filename}")
            numbers = re.findall(r'\d+', filename)
            z_layer = numbers[0] if numbers else "Unknown"
             
            full_file_path = os.path.join(in_directory, filename)
            slice_2d = np.load(full_file_path, allow_pickle=True)
             
            plt.figure(figsize=(8, 8))
            plt.imshow(slice_2d, cmap='RdYlBu')
            plt.title(f"Master Volume - Global Z Slice: {z_layer}")
            plt.axis('off')
             
            output_filename = f"plot_Z{z_layer}.png"
            full_output_path = os.path.join(out_directory, output_filename)
            
     
        #  Save it in_directory = ('foldername')
        out_directory = ('foldername')
        os.makedirs(out_directory, exist_ok=True)

        CUBE_SIZE = 200

        print(f"Reading slices from: {in_directory}\n")
        print(f"saving images to: {out_directory}\n")

        def check_file_exists(path):
            try:
                os.stat(path)
            except (OSError, ValueError):
                return False
            return True 
        file_count = 0

        for filename in os.listdir(in_directory):   
            if filename.endswith('.npy'):
                print(f"processing files : {filename}")
                numbers = re.findall(r'\d+', filename)
                z_layer = numbers[0] if numbers else "Unknown"
                 
                full_file_path = os.path.join(in_directory, filename)
                slice_2d = np.load(full_file_path, allow_pickle=True)
                 
                plt.figure(figsize=(8, 8))
                plt.imshow(slice_2d, cmap='RdYlBu')
                plt.title(f"Master Volume - Global Z Slice: {z_layer}")
                plt.axis('off')
                 
                output_filename = f"plot_Z{z_layer}.png"
                full_output_path = os.path.join(out_directory, output_filename)
                
         
            #  Save it 
                plt.savefig(full_output_path, bbox_inches='tight')
                plt.close()
            plt.savefig(full_output_path, bbox_inches='tight')
            plt.close()

    print(f"Successfully saved slice to {full_output_path}")

def cmap3d(in_directory, out_directory):
    in_directory = ('foldername')
    out_directory = ('foldername')
    os.makedirs(out_directory, exist_ok=True)

    CUBE_SIZE = 200

    print(f"Reading slices from: {in_directory}\n")
    print(f"saving images to: {out_directory}\n")

    def check_file_exists(path):
        try:
            os.stat(path)
        except (OSError, ValueError):
            return False
        return True 
    file_count = 0

    for filename in os.listdir(in_directory):   
        if filename.endswith('.npy'):
            print(f"processing files : {filename}")
            numbers = re.findall(r'\d+', filename)
            if len(numbers) >= 3:
                cube_x, cube_y, cube_z = int(numbers[0]), int(numbers[1]), int(numbers[2])
            
                global_x = cube_x * CUBE_SIZE
                global_y = cube_y * CUBE_SIZE
                global_z = cube_z * CUBE_SIZE
            
            else:
                print("skipped")
                continue
            
            full_file_path = os.path.join(in_directory, filename)
            sdf_data = np.load(full_file_path, allow_pickle=True)
            
            mid_z = sdf_data.shape[2] // 2
            slice_2d = sdf_data[:, :, mid_z]
            
            # Plot and save the slice as a viewable picture
            plt.figure(figsize=(6, 6))
            plt.imshow(slice_2d, cmap='RdYlBu')
            plt.title(f"Origin: X={global_x}, Y={global_y}, Z={global_z}")
            plt.axis('off')
            
            output_filename = f"slice_X{global_x}_Y{global_y}_Z{global_z}.png"
            full_output_path = os.path.join(out_directory, output_filename)
            print(f"Successfully saved slice to {full_output_path}")
