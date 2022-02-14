# -*- coding: utf-8 -*-
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

def getRigEvents(
):
    """
    Kim added this code for detecting cue, pellet approaching, and distractor events.

    Output: 
        
    """
    
    print("Starting getRigEvents")
    
    # codedir = r'C:\Users\sabatini\Documents\GitHub\integrate-phys-and-beh'
    codedir = os.getcwd()
    videos = [r'Z:\Kim\for_orchestra\deeplabcut test vids\dLight2_2021-02-10']
    # videos: list
    #    A list of strings containing the full paths to videos for analysis or a path to the directory, where all the videos with same extension are stored.
    getDifferenceEventsBeforeSaving = True
    cueThresh = 100
    rawThresh = 1
    distractorThresh = 100
    
    videotype = ".avi"
    wheelZone = [139, 223, 1, 91] # [x_start, x_end, y_start, y_end]
    cueZone = [1, 71, 295, 420] # [x_start, x_end, y_start, y_end]
    distractorZone = [1, 30, 452, 478] # [x_start, x_end, y_start, y_end]
    
    ##################################################
    # Looping over videos
    ##################################################
    print("Finding videos")
    Videos = auxiliaryfunctions.Getlistofvideos(videos, videotype)
    Videos = sorted(Videos)
    
    if len(Videos) > 0:
        rawDiffs = []
        cueDiffs = []
        distractorDiffs = []
        wheelDiffs = []
        for video in Videos:
            ##################################################
            # Loading the video
            ##################################################
            print("Starting to analyze % ", video)
            destfolder = None
            if destfolder is None:
                destfolder = str(Path(video).parents[0])
            auxiliaryfunctions.attempttomakefolder(destfolder)
            print("Loading ", video)
            cap = cv2.VideoCapture(video)
            if not cap.isOpened():
                raise IOError(
                    "Video could not be opened. Please check that the the file integrity."
                )
            fps = cap.get(
                5
            )  # https://docs.opencv.org/2.4/modules/highgui/doc/reading_and_writing_images_and_video.html#videocapture-get
            nframes = int(cap.get(7))
            duration = nframes * 1.0 / fps
            size = (int(cap.get(4)), int(cap.get(3)))
    
            ny, nx = size
            print(
                "Duration of video [s]: ",
                round(duration, 2),
                ", recorded with ",
                round(fps, 2),
                "fps!",
            )
            print(
                "Overall # of frames: ",
                nframes,
                " found with (before cropping) frame dimensions: ",
                nx,
                ny,
            )
            
            pbar = tqdm(total=nframes)
            counter = 0
            step = max(10, int(nframes / 100))
            prevFrame = None
            while cap.isOpened():
                if counter % step == 0:
                    pbar.update(step)

                ret, frame = cap.read()
                if ret:
                    ny, nx, nc = np.shape(frame)
                    frame = rgb2gray(frame)
                    ny, nx = np.shape(frame)
                    
                    # get the difference in raw pixels from frame to frame
                    # IN THE WHEEL ZONE
                    # a big jump suggests a fixed time before the cue, because acquisition was triggered 
                    # on pellet moving into position
                    # OR OVER THE WHOLE FRAME, NOT THE WHEEL ZONE 
                    # a large intensity change across the image suggests distractor
                    # OR HAVE USER SELECT A DISTRACTOR ZONE
                    # OR IN THE CUE ZONE
                    # a large intensity change in the cue zone suggests cue
                    if prevFrame is None:
                        prevFrame = frame
                    else:
                        rawDiffs.append(np.mean(frame - prevFrame))
                        cueDiffs.append(np.mean(frame[cueZone[0]:cueZone[1], cueZone[2]:cueZone[3]]  - prevFrame[cueZone[0]:cueZone[1], cueZone[2]:cueZone[3]]))
                        distractorDiffs.append(np.mean(frame[distractorZone[0]:distractorZone[1], distractorZone[2]:distractorZone[3]]  - prevFrame[distractorZone[0]:distractorZone[1], distractorZone[2]:distractorZone[3]]))
                        wheelDiffs.append(np.mean(frame[wheelZone[0]:wheelZone[1], wheelZone[2]:wheelZone[3]]  - prevFrame[wheelZone[0]:wheelZone[1], wheelZone[2]:wheelZone[3]]))
                        prevFrame = frame
                elif counter >= nframes:
                    break
                counter += 1
            pbar.close()
    # once have tested a video so know the right thresholds, can set getDifferenceEventsBeforeSaving to True
    # in order to save smaller files for Matlab
    if getDifferenceEventsBeforeSaving:
        def condition(x, thresh): return x > thresh
        def condminus(x, thresh): return x < -thresh
        rawDiffEvs_plus = [idx for idx, element in enumerate(rawDiffs) if condition(element, rawThresh)]
        cueDiffEvs_plus = [idx for idx, element in enumerate(cueDiffs) if condition(element, cueThresh)]
        distractorDiffEvs_plus = [idx for idx, element in enumerate(distractorDiffs) if condition(element, distractorThresh)]
        rawDiffEvs_minus = [idx for idx, element in enumerate(rawDiffs) if condminus(element, rawThresh)]
        cueDiffEvs_minus = [idx for idx, element in enumerate(cueDiffs) if condminus(element, cueThresh)]
        distractorDiffEvs_minus = [idx for idx, element in enumerate(distractorDiffs) if condminus(element, distractorThresh)]
        howmanyframes = len(cueDiffs)
        # save event indices
        outp = {"rawDiffEvs_plus": rawDiffEvs_plus, "cueDiffEvs_plus": cueDiffEvs_plus, "distractorDiffEvs_plus": distractorDiffEvs_plus,
                "rawDiffEvs_minus": rawDiffEvs_minus, "cueDiffEvs_minus": cueDiffEvs_minus, "distractorDiffEvs_minus": distractorDiffEvs_minus, "howmanyframes": howmanyframes}
    else:
        # save results
        outp = {"rawDiffs": rawDiffs, "cueDiffs": cueDiffs, "distractorDiffs": distractorDiffs, "wheelDiffs": wheelDiffs}
    savemat("rig_events.mat", outp)
    os.chdir(codedir)
    return nframes, rawDiffs, cueDiffs, distractorDiffs, wheelDiffs
            

def rgb2gray(rgb):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
    return gray

getRigEvents()