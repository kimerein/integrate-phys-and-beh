"""
Created on Mon Jan 24 10:48:54 2022

@author: sabatini
"""

import os
import numpy as np
import cv2
from tqdm import tqdm
from scipy.io import savemat
from deeplabcut.utils import auxiliaryfunctions
from pathlib import Path

def twoD_to_threeD(
):
    """
    2D DeeplabCut to 3D coordinates
    
    Before running this, need to check that 3D config.yaml file has the correct bodyparts and correct camera names.
    Change bodyparts in 2D config.yaml files.
    Change the project path in the config.yaml file for 2D if copied cameras.
    """

    bodyparts = ['currentPellet','middle_knuckle','middle_tip','palm','pinky_knuckle','pinky_tip','pointer_knuckle','pointer_tip','ring_knuckle','ring_tip','thumb_knuckle','thumb_tip','wrist']
    df = pd.read_hdf(r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\March_C_2019-05-30\camera-1_2019-05-30-152018-0001DLC_resnet50_Testing2DJan4shuffle1_500000.h5')

    # Rename columns
    df_new = df
    for b in bodyparts:
        df_new = df_new.rename(columns={b+'Side':b})
    
    for b in df_new.columns.get_level_values(1):
        if b not in bodyparts:
            df_new = df_new.drop(columns=b, level=1)

    # To check that it worked
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    df_new.head()

    # Write new file
    df_new.to_hdf(r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\March_C_2019-05-30\sidecam_2019-05-30-152018-0001DLC_resnet50_Testing2DJan4shuffle1_500000.h5', key='df_with_missing', format='table', mode='w')

    # Repeat for second camera
    df = pd.read_hdf(r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\March_C_2019-05-30\camera-1_2019-05-30-152018-0001DLC_resnet50_Testing2DJan4shuffle1_500000.h5')

    # Rename columns
    df_new = df
    for b in bodyparts:
        df_new = df_new.rename(columns={b+'Side':b})
    
    for b in df_new.columns.get_level_values(1):
        if b not in bodyparts:
            df_new = df_new.drop(columns=b, level=1)

    # To check that it worked
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    df_new.head()

    # Write new file
    df_new.to_hdf(r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\March_C_2019-05-30\sidecam_2019-05-30-152018-0001DLC_resnet50_Testing2DJan4shuffle1_500000.h5', key='df_with_missing', format='table', mode='w')

twoD_to_threeD()