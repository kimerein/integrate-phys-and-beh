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
    Essentially, these assumptions are required to construct the camera matrix, distortion coefficients and principal points for each view
"""

import os
import pickle

import cv2
import scipy
import numpy as np

from deeplabcut.utils import auxiliaryfunctions
from deeplabcut.utils import auxiliaryfunctions_3d

def main():
    config=r'C:\Users\sabatini\Documents\DLC_3D\Test3D-Whitney-2022-06-21-3d\config.yaml'
    testimg=r'C:\Users\sabatini\Documents\test_vid_july3.jpg'
     
    # this code is only for one camera/mirror/view pair, called for example 'sidecam' and 'undercam'
    # can rerun this calibration for a different camera pair
    views=['side', 'under']
    camera_names=[views[0]+'cam', views[1]+'cam']
    
    # define files to import from Matlab
    objectpointsfile=r'C:\Users\sabatini\Documents\objpoints.mat'
    sideimgpointsfile=os.path.join(r'C:\Users\sabatini\Documents', ''.join([views[0], 'points.mat']))
    underimgpointsfile=os.path.join(r'C:\Users\sabatini\Documents', ''.join([views[1], 'points.mat']))
    mirrorcentersfile=r'C:\Users\sabatini\Documents\mirrorcenters.mat'
    imgsizefile=r'C:\Users\sabatini\Documents\imgsize.mat'
    
    # known facts about camera / rig setup
    # paw_length=8 # in mm
    paw_length_in_pixels=70
    camera_focal_length_pixels=2*paw_length_in_pixels
    f_x=camera_focal_length_pixels # assume same focal lengths in x and y
    f_y=camera_focal_length_pixels
    obj_x_offset=0 #200 #70 # in mm
    obj_y_offset=0 #100 #66 # in mm
    obj_z_offset=0 #60 #20 # in mm
    
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
    stereo_from_uncalib={}
    stereo_params={}
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
    #objp=objp-np.min(objp)
    for i in range(len(temp)):
        objp[i][0]=objp[i][0]+obj_x_offset # shift object world coordinates in x
        objp[i][1]=objp[i][1]+obj_y_offset # shift object world coordinates in y
        objp[i][2]=objp[i][2]+obj_z_offset # shift object world coordinates in z
    print('object points')
    print(objp)
    for cam in camera_names:
        objpoints[cam].append(objp)
    
    # image points view 1
    mat_imgpoints=scipy.io.loadmat(sideimgpointsfile)
    objp=np.zeros((len(mat_imgpoints[views[0]+'points']),2),np.float32)
    temp=mat_imgpoints[views[0]+'points']
    for i in range(0,len(temp)):
        a=temp[i]
        temp2=a
        objp[i]=temp2
    imgpoints[camera_names[0]].append(objp)
    
    # image points view 2
    mat_imgpoints=scipy.io.loadmat(underimgpointsfile)
    objp=np.zeros((len(mat_imgpoints[views[1]+'points']),2),np.float32)
    temp=mat_imgpoints[views[1]+'points']
    for i in range(0,len(temp)):
        a=temp[i]
        temp2=a
        objp[i]=temp2
    imgpoints[camera_names[1]].append(objp)
    
    # mirror centers are estimates of the prinicipal points
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
        ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints[cam],
                                                           imgpoints[cam],
                                                           [imagesize['dim1'],imagesize['dim2']],
                                                           camera_matrix,distCoeffs,
                                                           None,
                                                           None,
                                                           flags=cv2.CALIB_USE_INTRINSIC_GUESS+
                                                           cv2.CALIB_FIX_PRINCIPAL_POINT+
                                                           cv2.CALIB_ZERO_TANGENT_DIST+
                                                           cv2.CALIB_FIX_FOCAL_LENGTH+
                                                           cv2.CALIB_FIX_K1+
                                                           cv2.CALIB_FIX_K2+
                                                           cv2.CALIB_FIX_K3+
                                                           cv2.CALIB_FIX_K4+
                                                           cv2.CALIB_FIX_K5+
                                                           cv2.CALIB_FIX_K6)
        print("camera matrix for " + cam)
        print(mtx)
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
    
    # get the reprojected points consistent with one homography 
    for cam in camera_names:
        testimgpoints, jcb = cv2.projectPoints(objpoints[cam][0],dist_pickle[cam]["rvecs"][0],dist_pickle[cam]["tvecs"][0],dist_pickle[cam]["mtx"],dist_pickle[cam]["dist"])
        dist_pickle[cam].update({"reproj_points": testimgpoints})
        
    # get the fundamental matrix for stereo calibration
    # F, mask = cv2.findFundamentalMat(dist_pickle[camera_names[0]]["reproj_points"], dist_pickle[camera_names[1]]["reproj_points"], cv2.FM_8POINT) # use nearly all the points
    F, mask = cv2.findFundamentalMat(dist_pickle[camera_names[0]]["reproj_points"], dist_pickle[camera_names[1]]["reproj_points"], cv2.FM_RANSAC, 0.1, 0.99) # should use all if using reproj points
    stereo_from_uncalib["F"] = F
    # get and save homography for reference although don't use it
    H = cv2.stereoRectifyUncalibrated(dist_pickle[camera_names[0]]["reproj_points"],dist_pickle[camera_names[1]]["reproj_points"],stereo_from_uncalib["F"],(imagesize['dim1'], imagesize['dim2']),0)
    stereo_from_uncalib["H1"] = H[1]
    stereo_from_uncalib["H2"] = H[2]
    # to 3D deeplabcut
    pickle.dump(
        stereo_from_uncalib,
        open(
            os.path.join(path_camera_matrix, "stereo_from_uncalib.pickle"),
            "wb",
        ),
    )
    print(
        "Saving fundamental matrix and homographies for stereo calibration"
    )
    
    # get essential matrix
    focal_length=f_x
    E = np.matmul(np.matmul(np.matrix.transpose(dist_pickle[camera_names[1]]["mtx"]),stereo_from_uncalib["F"]*focal_length),dist_pickle[camera_names[0]]["mtx"])
    trydecomp=cv2.decomposeEssentialMat(E)
    try1_R=trydecomp[0]
    try2_R=trydecomp[1]
    t=trydecomp[2]
    # produces four possible poses: [try1_R,t], [try1_R,-t], [try2_R,t], [try2_R,-t]
    
    img1=cv2.imread(testimg)
    cv2.imwrite(os.path.join(path_camera_matrix, 'input_image_for_test_undistort.jpg'),img1)
    pair = [camera_names[0], camera_names[1]]
    print(
        "Saving the stereo parameters for pair of cameras as a pickle file in %s"
        % str(os.path.join(path_camera_matrix))
    )
    for i in range(4):
        # stereorectify
        if i==0:
            tryrot=try1_R
            tryt=t
        if i==1:
            tryrot=try1_R
            tryt=-t
        if i==2:
            tryrot=try2_R
            tryt=t
        if i==3:
            tryrot=try2_R
            tryt=-t
        R1, R2, P1, P2, Q, roi1, roi2 = cv2.stereoRectify(dist_pickle[camera_names[0]]["mtx"], 
                                                          dist_pickle[camera_names[0]]["dist"], 
                                                          dist_pickle[camera_names[1]]["mtx"], 
                                                          dist_pickle[camera_names[1]]["dist"], 
                                                          (imagesize['dim1'], imagesize['dim2']), 
                                                          tryrot, tryt, alpha=0.9)
        map1_x, map1_y = cv2.initUndistortRectifyMap(dist_pickle[camera_names[0]]["mtx"],
                                                     dist_pickle[camera_names[0]]["dist"],
                                                     R1,
                                                     P1,
                                                     (imagesize['dim1'], imagesize['dim2']),
                                                     cv2.CV_16SC2)
        map2_x, map2_y = cv2.initUndistortRectifyMap(dist_pickle[camera_names[1]]["mtx"],
                                                     dist_pickle[camera_names[1]]["dist"],
                                                     R2,
                                                     P2,
                                                     (imagesize['dim1'], imagesize['dim2']),
                                                     cv2.CV_16SC2)
        im_remapped = cv2.remap(img1, map1_x, map1_y, cv2.INTER_LINEAR, cv2.BORDER_TRANSPARENT)
        # count up how many valid mapped points
        #validpts = 0
        #for j in range(len(im_remapped)):
        #    for k in range(len(im_remapped[j])):
        #        validpts = validpts + np.any(im_remapped[j][k])
        #print('valid mapped points for view 1')
        #print(validpts)
        cv2.imwrite(os.path.join(path_camera_matrix, 'test_undistort_pose' + str(i) + 'A.jpg'),im_remapped)
        im_remapped = cv2.remap(img1, map2_x, map2_y, cv2.INTER_LINEAR, cv2.BORDER_TRANSPARENT)
        # count up how many valid mapped points
        #validpts = 0
        #for j in range(len(im_remapped)):
        #    for k in range(len(im_remapped[j])):
        #        validpts = validpts + np.any(im_remapped[j][k])
        #print('valid mapped points for view 2')
        #print(validpts)
        cv2.imwrite(os.path.join(path_camera_matrix, 'test_undistort_pose' + str(i) + 'B.jpg'),im_remapped)
        
        # save for pickle
        stereo_params[pair[0] + "-" + pair[1]] = {"cameraMatrix1": dist_pickle[camera_names[0]]["mtx"],
                                                  "cameraMatrix2": dist_pickle[camera_names[1]]["mtx"],
                                                  "distCoeffs1": dist_pickle[camera_names[0]]["dist"],
                                                  "distCoeffs2": dist_pickle[camera_names[1]]["dist"],
                                                  "R": tryrot,
                                                  "T": tryt,
                                                  "E": E,
                                                  "F": stereo_from_uncalib["F"],
                                                  "R1": R1,
                                                  "R2": R2,
                                                  "P1": P1,
                                                  "P2": P2,
                                                  "roi1": roi1,
                                                  "roi2": roi2,
                                                  "Q": Q,
                                                  "image_shape":[(imagesize['dim1'], imagesize['dim2']), (imagesize['dim1'], imagesize['dim2'])]}
        # write stereo params for all possible poses
        # user selects the correct pose and changes one of these files to stereo_params.pickle for 3d deeplabcut
        auxiliaryfunctions.write_pickle(
            os.path.join(path_camera_matrix, "stereo_params" + str(i) + ".pickle"), stereo_params
        )
    
        # get error of this stereo calibration
        # convert all to world coordinate points in view 1's rectified world coordinate reference frame
    
        # triangulate
        # note that output will be in view 1's rectified world coordinates, since P1, R1, P2, R2 are from stereorectify
        points1=dist_pickle[camera_names[0]]['reproj_points']
        points2=dist_pickle[camera_names[1]]['reproj_points']
        testpoints1 = cv2.undistortPoints(src=points1.astype(np.float32),cameraMatrix=dist_pickle[camera_names[0]]["mtx"], distCoeffs=dist_pickle[camera_names[0]]["dist"],P=stereo_params[pair[0] + "-" + pair[1]]["P1"],R=stereo_params[pair[0] + "-" + pair[1]]["R1"])
        testpoints2 = cv2.undistortPoints(src=points2.astype(np.float32),cameraMatrix=dist_pickle[camera_names[1]]["mtx"], distCoeffs=dist_pickle[camera_names[1]]["dist"],P=stereo_params[pair[0] + "-" + pair[1]]["P2"],R=stereo_params[pair[0] + "-" + pair[1]]["R2"])
        X = cv2.triangulatePoints(stereo_params[pair[0] + "-" + pair[1]]["P1"][:3],stereo_params[pair[0] + "-" + pair[1]]["P2"][:3], testpoints1, testpoints2)
        X = X / X[3]
        # rotate X
        # this step specific to exact rig orientation of mirrors and cameras
        rotateX=np.matmul([[1, 0, 0],[0, -1, 0],[0, 0, 1]], np.matmul([[1, 0, 0],[0, 1, 0],[0, 0, 1]],X[:3]))
        savedata={}
        savedata["triangulationOutput"]=rotateX
        
        # use Rodrigues to get rotations from the world coordinates to each camera's rectified world coordinates
        outOfRod = cv2.Rodrigues(dist_pickle[camera_names[0]]["rvecs"][0])
        Rotate = outOfRod[0]
        outOfRod = cv2.Rodrigues(dist_pickle[camera_names[1]]["rvecs"][0])
        Rotate2 = outOfRod[0]
    
        # input from view 2 to view 1's rectified world coordinates
        # plus a final rotation specific to exact rig orientation of mirrors and cameras
        view2to1rect=np.matmul([[1, 0, 0],[0, -1, 0],[0, 0, 1]], 
                  np.matmul(stereo_params[pair[0] + "-" + pair[1]]["R1"], 
                            np.matmul(np.linalg.inv(stereo_params[pair[0] + "-" + pair[1]]["R"]), 
                                      np.matmul(Rotate2,
                                                np.matrix.transpose(dist_pickle[camera_names[1]]['objpoints'][0])))))
        savedata["view2to1rect"]=view2to1rect
        
        # input from view 1 to view 1's rectified world coordinates 
        # plus a final rotation specific to exact rig orientation of mirrors and cameras
        view1to1rect=np.matmul(stereo_params[pair[0] + "-" + pair[1]]["R1"],
                  np.matmul(Rotate,np.matrix.transpose(dist_pickle[camera_names[0]]['objpoints'][0])))
        savedata["view1to1rect"]=view1to1rect
        scipy.io.savemat(os.path.join(path_camera_matrix, 'compareStereos' + str(i) + '.mat'),savedata)
        
        
main()
