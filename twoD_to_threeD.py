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
    
    # For each .h5 file, rename the columns to match the bodyparts in the 3D config.yaml file
    # Read in all the .h5 files from directory
    directory=r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\March_C_2019-05-30'
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".h5"):
            print(['Now processing ' os.path.join(directory, filename)])
            df = pd.read_hdf(os.path.join(directory, filename))

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
            currfname=os.path.join(directory, filename)
            # Change camera-1 to sidecam in filename
            newfname=currfname.replace('camera-1','sidecam')
            df_new.to_hdf(newfname, key='df_with_missing', format='table', mode='w')

            # Repeat for second camera
            # Rename columns
            df_new = df
            for b in bodyparts:
                df_new = df_new.rename(columns={b+'Bottom':b})
            
            for b in df_new.columns.get_level_values(1):
                if b not in bodyparts:
                    df_new = df_new.drop(columns=b, level=1)

            # To check that it worked
            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_rows', None)
            df_new.head()

            # Write new file
            # Change camera-1 to sidecam in filename
            newfname=currfname.replace('camera-1','undercam')
            df_new.to_hdf(newfname, key='df_with_missing', format='table', mode='w')
        else:
            continue
    
twoD_to_threeD()