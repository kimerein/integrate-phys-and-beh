# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 10:11:24 2022

@author: Kim

Script stereo calibrates a pair of mirrors (or views) imaged by a single camera, with assumptions:
    Camera distortion is nearly zero (i.e., imaged field of view is small enough such that no fisheye/pinhole camera distortion)
    Mirrors introduce essentially no distortion
    Relative positions of mirrors and camera are similar to Kim's rig (this pertains only to the final reference frame rotations)
    Fixed camera focal length around 13-16 mm, imaged fov around 4x4x4 cm, and/or changing camera focal length creates minimal distortion (small enough fov)
    Known relationship between real paw length in world and number of pixels in camera image
    Known approximate prinicipal points for the two mirror views
    Essentially, these assumptions are required to construct the camera matrix and distortion coefficients for each view
"""

import glob
import os
import pickle
from pathlib import Path

import cv2
import scipy
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes._axes import _log as matplotlib_axes_logger

from deeplabcut.utils import auxiliaryfunctions
from deeplabcut.utils import auxiliaryfunctions_3d

matplotlib_axes_logger.setLevel("ERROR")

config=r'C:\Users\sabatini\Documents\DLC_3D\Test2022-Kim-2022-01-04-3d\config.yaml'

# this code is only for one camera/mirror/view pair, called for example 'sidecam' and 'undercam'
# can rerun this calibration for a different camera pair
views=['side', 'under']
camera_names=[views[0]+'cam', views[1]+'cam']

# define files to import from Matlab
objectpointsfile=r'C:\Users\sabatini\Documents\objpoints.mat'
sideimgpointsfile=os.path.join(r'C:\User\sabatini\Documents', ''.join([views[0], 'points.mat']))
underimgpointsfile=os.path.join(r'C:\User\sabatini\Documents', ''.join([views[1], 'points.mat']))
mirrorcentersfile=r'C:\Users\sabatini\Documents\mirrorcenters.mat'
imgsizefile=r'C:\Users\sabatini\Documents\imgsize.mat'

# known facts about camera / rig setup
paw_length=8 # in mm
paw_length_in_pixels=70
camera_focal_length_pixels=2*paw_length_in_pixels
f_x=camera_focal_length_pixels # assume same focal lengths in x and y
f_y=camera_focal_length_pixels

# read image size
imagesize={}
imagesize.setdefault('dim1')
imagesize.setdefault('dim2')
ims=scipy.io.loadmat(imgsizefile)
imagesize['dim1']=ims['imagesize'][0][0] # x dimension
imagesize['dim2']=ims['imagesize'][1][0] # y dimension

# Read the config file
cfg_3d = auxiliaryfunctions.read_config(config)
(
    img_path,
    path_corners,
    path_camera_matrix,
    path_undistort,
    path_removed_images,
) = auxiliaryfunctions_3d.Foldernames3Dproject(cfg_3d)

# get object and image points both views, mirror centers
objpoints={}
imgpoints={}
mirrorcenters={}
dist_pickle={}
for cam in camera_names:
    objpoints.setdefault(cam,[])
    imgpoints.setdefault(cam,[])
    mirrorcenters.setdefault(cam,[])
    dist_pickle.setdefault(cam,[])

# object points
mat_objpoints=scipy.io.loadmat(objectpointsfile)
objp=np.zeros((len(mat_objpoints['objpoints']),3),np.float32)
temp=mat_objpoints['objpoints']
for i in range(0,len(temp)):
    a=temp[i]
    temp2=a
    objp[i]=temp2
for cam in camera_names:
    objpoints[cam]=objp

# image points view 1
mat_imgpoints=scipy.io.loadmat(sideimgpointsfile)
objp=np.zeros((len(mat_imgpoints[views[0]+'points']),3),np.float32)
temp=mat_imgpoints[views[0]+'points']
for i in range(0,len(temp)):
    a=temp[i]
    temp2=a
    objp[i]=temp2
imgpoints[camera_names[0]]=objp

# image points view 2
mat_imgpoints=scipy.io.loadmat(underimgpointsfile)
objp=np.zeros((len(mat_imgpoints[views[1]+'points']),3),np.float32)
temp=mat_imgpoints[views[1]+'points']
for i in range(0,len(temp)):
    a=temp[i]
    temp2=a
    objp[i]=temp2
imgpoints[camera_names[1]]=objp

# mirror centers are the prinicipal points
mat_mirrorcenters=scipy.io.loadmat(mirrorcentersfile)
mirrorcenters[camera_names[0]] = {
        "p_x": mat_mirrorcenters['mirror_vs'][0][0],
        "p_y": mat_mirrorcenters['mirror_vs'][0][1] 
    }
mirrorcenters[camera_names[1]] = {
        "p_x": mat_mirrorcenters['mirror_vs'][0][2],
        "p_y": mat_mirrorcenters['mirror_vs'][0][3]  
    }

for cam in camera_names:
    # camera/mirror/view calibrations
    camera_matrix=[[f_x, 0, mirrorcenters[cam]['p_x']],[0, f_y, mirrorcenters[cam]['p_y']],[0, 0, 1]];
    camera_matrix=np.array(camera_matrix)
    distCoeffs=np.array([0,0,0,0,0]) # assume zero distortion
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints[cam],imgpoints[cam],[imagesize['dim1'],imagesize['dim2']],camera_matrix,distCoeffs,None,None,flags=cv2.CALIB_USE_INTRINSIC_GUESS)
    # save for pickle 
    dist_pickle[cam] = {
                "mtx": mtx,
                "dist": dist,
                "rvecs": rvecs,
                "tvecs": tvecs,
                "objpoints": objpoints[cam],
                "imgpoints": imgpoints[cam],
            }
    # to 3D deeplabcut
    pickle.dump(
        dist_pickle,
        open(
            os.path.join(path_camera_matrix, cam + "_intrinsic_params.pickle"),
            "wb",
        ),
    )
    print(
        "Saving intrinsic camera calibration matrices for %s as a pickle file in %s"
        % (cam, os.path.join(path_camera_matrix))
    )

# get the reprojected points consistent with a single self-consistent homography 

