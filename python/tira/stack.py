import os
import numpy as np
import skimage as ski


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
