
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['figure.dpi'] = 300

import time, os, sys, math, argparse
from types import SimpleNamespace
from pprint import pprint
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from cellpose import utils, models, plot, io

import skimage.io
import cv2
from skimage import data, img_as_int, img_as_uint, img_as_float, filters, feature, morphology, color, exposure, measure
from scipy import ndimage as nd

import glob
import czifile
 
    
#########################################################
                   ## Main Function ##
#########################################################
    
def main(params):
    
    # Make output path if it does not exist
    if not os.path.isdir(params.output_path):
        os.mkdir(params.output_path)
        
    # parse directory
    file_list = parse_input_path(params)
    params.file_list = file_list

    # segment nuclei
    segment_nucleus_3D(params, preprocessing=True, save=True)
              

#########################################################
                ## FUNCTION DEFINITIONS ##
#########################################################

def parse_arguments(parser):

    # required arguments
    parser.add_argument('parent_path', type=str,
                        help='Full path to folder that contains subfolders of experiments with data')
    parser.add_argument('output_path', type=str,
                        help='Full path to folder where you want the output to be stored. The folder will be made if it does not exist')
    parser.add_argument('channel', type=str,
                        help='Unique ID for channel that corresponds to nuclei')                        
    parser.add_argument('diameter', type=int,
                        help='Approximate diameter of cells in the image in pixels.')   
    
    # optional arguments
    parser.add_argument("--min_area", type=float, default=10000,
                        help='Optional minimum area for nuclei in pixels^2. Any nuclei below this will be discarded. Default is 10000')

    params = parser.parse_args()

    return params


def parse_input_path(params):
    
    '''
     Inputs
    ---------
    parent_path: This function will search this entire directoy (including subdirectories) to find TIFs that contain the
    channel name
    
    channel: This is a string corresponding to a unique identifier for this channel (e.g. "488 GFP", "561 CY3", "642 CY5")
    
     Outputs
    ---------
    file_list:  is a sorted list of file paths specifying the images used for nuclei detection. This will be the direct input
    for cellpose.
    '''
    
    parent_path = params.parent_path
    channel = params.channel
    
    file_list = []
    sample_names = []
    file_extension = get_microscope_details(parent_path,params)
    
    # For RPI/Andor
    if file_extension.lower() == ".tif":
        for root, dirs, files in os.walk(parent_path):
            for file in files:
                if any([file_extension in file]) and channel in file:
                    file_list.append(os.path.join(root, file))
    # For Airyscan
    elif file_extension.lower() == ".czi":
        for root, dirs, files in os.walk(parent_path):
            for file in files:
                if file_extension in file:
                    file_list.append(os.path.join(root, file))
    
    file_list.sort()
    
    return file_list


def load_images_czi(file): # Order of img shape is (B, V, C, T, Z, Y, X, 0)
    
    with czifile.CziFile(file) as czi:
        img = czi.asarray()
        num_of_channels = img.shape[2]
        print(num_of_channels)
      
        
        img = np.squeeze(img)  # gets rid of all singleton dimensions
        print(img.shape)
        
        # load in images
        channel_images = []
        for i in range(num_of_channels):
            if num_of_channels > 1:
                channel_images.append(img_as_float(img[i, :, :, :]))
            else:
                channel_images.append(img_as_float(img[:, :, :]))
        metadata = czi.metadata()

    root = ET.fromstring(metadata)

    channel_names = []
    excitations = []
    for channels in root.findall("./Metadata/Information/Image//Channels/"):
            channel_names.append("ch" + channels.attrib['Name'].replace(" ", "_"))
            excitation = channels.find("ExcitationWavelength").text
            excitations.append(excitation)
            
    
    print("Excitations: ", excitations)
    print("Found the following channel names", channel_names)
    
    for i in range(0,len(excitations)):
        if "633" in excitations[i]: #change based on excitation corresponding to nuclear image
            nuclear_img = channel_images[i]
        else:
            nuclear_img = channel_images[0]
    return nuclear_img

def segment_nucleus_3D(params, preprocessing=False, save=False):
    
    '''
      Inputs
    ------------
    
    file_list: is the output from parse_input_path. A list of file paths to the TIF images of nuclei
    
    diameter: is empirically-determined average diameter of the cells. Find this using ImageJ for now.
    In the future, we will automagically determine this.
    
    preprocessing: is a flag that when activated, uses some filters to make the nuclear calling easier
    
    save: is a flag that determines whether PNGs and python files are saved to output
    
    NOTE: The prior way of getting sample_names was wrong, because they were being sorted, which put them out of order
    with "file_list" because that list contains the full path, whereas sample_names was just the basename
    file_list = params.file_list
    
      Outputs
    ------------
    
    _seg.npy: file with nuclear segmentaiton for each sample outputted to the output path
    
    _Nuclei.png: image of called nuclei for each sample outputted to the output path
    
      Returns
    ------------
    
    NA
    
    '''
    
    file_list = params.file_list
    parent_path = params.parent_path
    delimiter = params.delimiter
    
    sample_names = [os.path.basename(f)[:os.path.basename(f).find(delimiter)].replace(".", "dot") for f in file_list]
    params.sample_names = sample_names
    
    diameter = params.diameter
    

    if params.file_extension == ".TIF" or params.file_extension == ".tif":
        imgs = [skimage.io.imread(f) for f in file_list]
    elif params.file_extension == ".czi":
        imgs = [load_images_czi(f) for f in file_list]     
    
    # Preprocess images with max and median filter
    if preprocessing:
        imgs = preprocess_imgs(imgs)
    
    # List of 2D images
    max_imgs = []
    for img in imgs:
        
        # 3D max projection
        if len(img.shape) == 3:  # if z-stack
            max_img = max_project(img)
            max_imgs.append(max_img)
        
        # 2D option
        else:
            max_imgs.append(img)
    
    # Generate cellpose model
    model = models.Cellpose(gpu=False, model_type='nuclei')
    channels = [0,0]
    
    # Run segmentation on list of 2D images
    masks, flows, styles, diams = model.eval(max_imgs, diameter=diameter, channels=channels)
    
    filtered_masks = []
    for m in masks:
        temp_mask = filter_nuclei_below_size_threshold(params, m)
        temp_mask = fix_mask_labels(temp_mask)
        
        filtered_masks.append(temp_mask)
        
    output_names = [os.path.join(params.output_path, s) for s in sample_names]
 
    io.masks_flows_to_seg(max_imgs, filtered_masks, flows, diams, output_names)
    
    # plotting
    nimgs = len(imgs)
    
    for n in range(nimgs):
        fig = plt.figure(figsize=(20, 10))
        masksi = filtered_masks[n]
        flowsi = flows[n][0]
        
        plot.show_segmentation(fig, max_imgs[n], masksi, flowsi, channels=channels)
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(params.output_path, sample_names[n] + '_Nuclei.png'), dpi=300)
        else:
            plt.show()
            
        plt.close()


# Old preprocess step
'''
def preprocess_imgs(imgs):
     Inputs
    --------
    
    imgs: is a list of the numpy arrays of nuclear imgs
    
    Outputs
    --------
    
    output: is a list of numpy arrays of the nuclear imgs with the maximum and median filter applied
    
    output = []
    for img in imgs:
        
        # 3D filtering
        if len(img.shape) == 3:
            output_img = nd.filters.maximum_filter(img, size=10) # this may throw errors if the z-dimension is smaller than 10
            output_img = nd.filters.median_filter(output_img, size=(1, 20, 20))
        
        # 2D filtering
        elif len(img.shape) == 2:
            output_img = nd.filters.maximum_filter(img, size=10)
            output_img = nd.filters.median_filter(output_img, size=10)
            
        output.append(output_img)

    return output
'''
    
def preprocess_imgs(imgs):
   '''                                                                                                     
   Inputs                                                                                                 
   --------                                                                                                
   imgs: is a list of the numpy arrays of nuclear imgs                                                        
                                                                                                        
   Outputs                                                                                                 
   --------                                                                                               
   output: is a list of numpy arrays of the nuclear imgs with the maximum and median filter applied
   '''

   # imgs is a list of the numpy arrays of nuclear imgs

   output = []
   for img in imgs:

       if len(img.shape) == 3:
           print("image shape is 3")

           output_img = nd.filters.median_filter(img, size=(1, 20, 20))
       elif len(img.shape) == 2:
           output_img = img
           threshold = np.mean(output_img) + 2*np.std(output_img)
           output_img[output_img > threshold] = np.mean(output_img)
           output_img = nd.filters.median_filter(output_img, size=10)
           print(output_img.shape)

       output.append(output_img)

   return output


def filter_nuclei_below_size_threshold(params, mask):
    
    # returns a labeled mask where nuclei that are smaller than params.area_threshold have their labels
    # changed to 0.
    
    img = mask.copy()
    num_of_nuclei = np.max(img)

    for n in np.arange(1, num_of_nuclei + 1):  # because the 0 label is background
        nuc_area = np.sum(img == n)
        if nuc_area < params.min_area:
            img[img == n] = 0

    return img


def fix_mask_labels(mask):
    
    # This function will take a label mask and renumber the labels so that they are contiguous
    
    unique_labels = np.unique(mask)
    x = len(unique_labels) - 1  # The number of labels should equal the max label # - 1 because 0 is included

    while np.max(mask) > x:
        for i in np.arange(0, x+1):  # The +1 is to include the last label
            if i > 0:
                label_test = np.any(mask == i)
                if not label_test:
                    mask[np.nonzero(mask > i)] = mask[np.nonzero(mask > i)] - 1

    return mask


def max_project(img):
    projection = np.max(img, axis=0)
    
    return projection


def find_img_channel_name(file_name):
    str_idx = file_name.find('Conf ')  
    channel_name = file_name[str_idx + 5 : str_idx + 8]
    channel_name = 'ch' + channel_name

    return channel_name


def find_region_area(r):
    a = ((r[0].stop - r[0].start)) * ((r[1].stop - r[1].start))
    
    return a
    

def filter_nuc_regions(nuc_mask, nuc_regions, params):
    for idx, region in enumerate(nuc_regions):
        if region.area < params.nuc_area_threshold:
            coords = region.coords
            row = [r[0] for r in coords]
            col = [c[1] for c in coords]
            
            nuc_mask[row, col] = False
        
    return nuc_mask


def set_region_to_zero(img, coords):
    output_img = img
    r_coords = [r[1] for r in coords]
    c_coords = [c[2] for c in coords]
    
    output_img[r_coords, c_coords] = 0
    
    return output_img

def get_microscope_details(parent_path, params):
    
    # This function will extract file extension, and delimiter used to get the sample name based on the file type
    
    if glob.glob(parent_path+"/*.TIF"):
        print("Microscope type: RPI or Andor")
        file_extension = ".TIF"
        delimiter = "_w"
    elif glob.glob(parent_path+"/*.tif"):
        print("Microscope type: RPI or Andor")
        file_extension = ".tif"
        delimiter = "_w"
    elif glob.glob(parent_path+"/*.czi"):
        print("Microscope type: Airyscan")
        file_extension = ".czi"
        delimiter = "_Airyscan" 
    else:
        print("Do not recognize file types")
        sys.exit(0)
    
    params.delimiter = delimiter
    params.file_extension = file_extension
    
    return file_extension


def clear_axis_ticks(ax):
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
        
                           
if __name__ == "__main__":

    # parse input
    parser = argparse.ArgumentParser()
    params = parse_arguments(parser)
    params.parent_path = params.parent_path.replace("Volumes","lab")
    params.output_path = params.output_path.replace("Volumes","lab")
    
    main(params)
        
    print('--------------------------------------')
    print('Completed at: ', datetime.now())
