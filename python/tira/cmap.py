import numpy as np

def normalize_intensity(image, min_val=None, max_val=None):
    fp_image = image.astype(np.float32)
    if min_val == None:
        min_val = np.min(fp_image)
    if max_val == None:
        max_val = np.max(fp_image)
    
    #  normalize the image and threshold the outliers
    norm_image = np.clip((fp_image - min_val)/(max_val - min_val), 0, 1)
    return norm_image

def green(image, min_val=None, max_val=None):
    fp_image = image.astype(np.float32)
    if min_val == None:
        min_val = np.min(fp_image)
    if max_val == None:
        max_val = np.max(fp_image)
      
    norm_image = normalize_intensity(fp_image, min_val, max_val)
    image_size = norm_image.shape
    color_image_size = image_size + (3,)
    
    color_image = np.zeros(color_image_size)
    color_image[:, :, 1] = norm_image
    return color_image