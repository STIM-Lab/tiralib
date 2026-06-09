import numpy as np
import skimage as ski

def npy2obj(input_file, output_name, level=0.0):
    """
    Convert a numpy array to an OBJ file.

    Args:
        input_file (str or np.ndarray): Input numpy array or file path.
        output_name (str): Output OBJ file path.
        level (float): Level for marching cubes.
    """
    if isinstance(input_file, str):
        volume = np.load(input_file)
    elif isinstance(input_file, np.ndarray):
        volume = input_file
    else:
        raise ValueError("Input must be a numpy array or file path.")
        
    # marching cubes
    verts, faces, _, _ = ski.measure.marching_cubes(volume, level=level)
    
    with open(output_name, 'w') as f:
        for v in verts:
            f.write(f'v {v[0]} {v[1]} {v[2]}\n')
        for face in faces:
            f.write(f'f {face[0]+1} {face[1]+1} {face[2]+1}\n')
