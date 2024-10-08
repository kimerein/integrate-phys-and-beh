# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:48:54 2022

@author: sabatini
"""

import os
import numpy as np
import cv2
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.io import savemat
#from deeplabcut.utils import auxiliaryfunctions
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
    #videos = r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20191127\April_short\DLC vids'
    videos = r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190510\3F_white\DLC output\test'
    # videos: list
    #    A list of strings containing the full paths to videos for analysis or a path to the directory, where all the videos with same extension are stored.
    getDifferenceEventsBeforeSaving = True
    
    videotype = ".avi"
    # dLight2_2021-02-10_00001.avi
    #wheelZone = [139, 223, 1, 91] # [x_start, x_end, y_start, y_end]
    #cueZone = [1, 71, 295, 420] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 30, 452, 478] # [x_start, x_end, y_start, y_end]
    #cueThresh = 100
    #rawThresh = 3
    #wheelThresh = 20
    #distractorThresh = 100
    # March_C first day
    #wheelZone = [200, 300, 1, 125] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 30, 340, 375] # [x_start, x_end, y_start, y_end]
    #reachZone = [106, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 10
    #wheelThresh = 2
    #rawThresh = 3
    #distractorThresh = 5
    #reachThresh = 3
    # March_C second day
    #wheelZone = [200, 300, 1, 185] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 30, 340, 375] # [x_start, x_end, y_start, y_end]
    #reachZone = [106, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 10
    #wheelThresh = 2
    #rawThresh = 3
    #distractorThresh = 5
    #reachThresh = 3
    # March_C 20190621
    #wheelZone = [200, 300, 1, 195] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 30, 340, 375] # [x_start, x_end, y_start, y_end]
    #reachZone = [106, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 10
    #wheelThresh = 2
    #rawThresh = 3
    #distractorThresh = 5
    #reachThresh = 3
    # March_C 20190715
    #wheelZone = [200, 300, 1, 195] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [6, 13, 405, 411] # [x_start, x_end, y_start, y_end]
    #reachZone = [106, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 2
    #rawThresh = 3
    #distractorThresh = 7
    #reachThresh = 3
    # March_A 20190624
    #wheelZone = [200, 300, 1, 195] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 21, 351, 381] # [x_start, x_end, y_start, y_end]
    #reachZone = [90, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1.5
    #rawThresh = 2
    #distractorThresh = 15
    #reachThresh = 3
    # March_A 20190725
    #wheelZone = [200, 300, 1, 195] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 125, 450, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 21, 351, 381] # [x_start, x_end, y_start, y_end]
    #reachZone = [90, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 15
    #reachThresh = 3
    # March_C 20190604
    #wheelZone = [200, 300, 1, 195] # [x_start, x_end, y_start, y_end]
    #cueZone = [75, 86, 492, 525] # [x_start, x_end, y_start, y_end]
    #distractorZone = [12, 23, 414, 421] # [x_start, x_end, y_start, y_end]
    #reachZone = [90, 148, 98, 144] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 15
    #reachThresh = 3
    # March_A 20190815
    #wheelZone = [191, 300, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [80, 109, 468, 543] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 36, 342, 389] # [x_start, x_end, y_start, y_end]
    #reachZone = [104, 150, 111, 198] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 15
    #reachThresh = 1.5
    # March_D 20190628
    #wheelZone = [207, 290, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [80, 134, 477, 559] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 28, 413, 453] # [x_start, x_end, y_start, y_end]
    #reachZone = [104, 150, 156, 208] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 3
    #reachThresh = 1.5
    # March_D 20190715
    #wheelZone = [207, 290, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [71, 111, 451, 559] # [x_start, x_end, y_start, y_end]
    #distractorZone = [19, 44, 528, 558] # [x_start, x_end, y_start, y_end]
    #reachZone = [88, 150, 102, 171] # [x_start, x_end, y_start, y_end]
    #cueThresh = 6
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 3
    #reachThresh = 1.5
    # March_D 20190715
    #wheelZone = [207, 262, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [0, 110, 501, 559] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 60, 338, 403] # [x_start, x_end, y_start, y_end]
    #reachZone = [110, 195, 115, 171] # [x_start, x_end, y_start, y_end]
    #cueThresh = 1
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 3
    #reachThresh = 1.5
    # April_short 20190927
    #wheelZone = [207, 262, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [63, 103, 477, 559] # [x_start, x_end, y_start, y_end]
    #distractorZone = [1, 49, 348, 395] # [x_start, x_end, y_start, y_end]
    #reachZone = [119, 195, 108, 155] # [x_start, x_end, y_start, y_end]
    #cueThresh = 3
    #wheelThresh = 1
    #rawThresh = 1
    #distractorThresh = 3
    #reachThresh = 1.5
    # April_short 20191011
    #wheelZone = [207, 262, 1, 251] # [x_start, x_end, y_start, y_end]
    #cueZone = [69, 140, 474, 559] # [x_start, x_end, y_start, y_end]
    #distractorZone = [38, 59, 347, 394] # [x_start, x_end, y_start, y_end]
    #reachZone = [107, 189, 108, 155] # [x_start, x_end, y_start, y_end]
    #cueThresh = 10
    #wheelThresh = 1
    #rawThresh = 2
    #distractorThresh = 2
    #reachThresh = 1.5
    # April_short 20191127
    #wheelZone = [150, 250, 1, 235] # [x_start, x_end, y_start, y_end]
    #cueZone = [59, 82, 446, 512] # [x_start, x_end, y_start, y_end]
    #distractorZone = [38, 59, 347, 394] # [x_start, x_end, y_start, y_end]
    #reachZone = [70, 142, 75, 120] # [x_start, x_end, y_start, y_end]
    #cueThresh = 5
    #wheelThresh = 1
    #rawThresh = 2
    #distractorThresh = 2
    #reachThresh = 1.5
    # Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190510\3F_white\DLC output\test
    wheelZone = [195, 254, 139, 190] # [x_start, x_end, y_start, y_end]
    cueZone = [1, 28, 432, 502] # [x_start, x_end, y_start, y_end]
    distractorZone = [1, 15, 341, 358] # [x_start, x_end, y_start, y_end]
    reachZone = [74, 124, 109, 171] # [x_start, x_end, y_start, y_end]
    optoZone = [113, 143, 412, 463] # [x_start, x_end, y_start, y_end]
    cueThresh = 5
    wheelThresh = 1.5
    rawThresh = 1
    distractorThresh = 30
    reachThresh = 7
    optoThresh = 20

    ##################################################
    # Looping over videos
    ##################################################
    print("Finding videos")
    Videos = []
    for filename in os.listdir(videos):
        if filename.endswith(videotype):
            Videos.append(filename)
    #Videos = auxiliaryfunctions.Getlistofvideos(videos, videotype)
    Videos = sorted(Videos)
    
    if len(Videos) > 0:
        rawDiffs = []
        cueDiffs = []
        distractorDiffs = []
        wheelDiffs = []
        vidStarts = []
        vidEnds = []
        reachDiffs = []
        optoDiffs = []
        # Go to location of videos
        os.chdir(videos)
        for video in Videos:
            ##################################################
            # Show zones on example frame
            ##################################################
            # If this is the first video
            if len(rawDiffs) == 0:
                plotZonesOnExampleFrame(video, wheelZone, cueZone, distractorZone, reachZone, optoZone)

            ##################################################
            # Loading the video
            ##################################################
            print("Starting to analyze ", video)
            #auxiliaryfunctions.attempttomakefolder(destfolder)
            print("Loading ", video)
            vidStarts.append(len(rawDiffs))
            cap = cv2.VideoCapture(video)
            if not cap.isOpened():
                raise IOError(
                    "Video could not be opened. Please check the file integrity."
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
                        reachDiffs.append(np.mean(frame[reachZone[0]:reachZone[1], reachZone[2]:reachZone[3]]  - prevFrame[reachZone[0]:reachZone[1], reachZone[2]:reachZone[3]]))
                        optoDiffs.append(np.mean(frame[optoZone[0]:optoZone[1], optoZone[2]:optoZone[3]]  - prevFrame[optoZone[0]:optoZone[1], optoZone[2]:optoZone[3]]))
                        prevFrame = frame
                elif counter >= nframes:
                    vidEnds.append(len(rawDiffs))
                    break
                counter += 1
            pbar.close()
    # Plot rawDiffs
    plt.figure()
    plt.plot(rawDiffs)
    # Plot rawThresh as line
    plt.plot([0, len(rawDiffs)], [rawThresh, rawThresh], 'r-')
    plt.plot([0, len(rawDiffs)], [-rawThresh, -rawThresh], 'r-')
    plt.title('rawDiffs')
    plt.show()
    # Plot wheelDiffs
    plt.figure()
    plt.plot(wheelDiffs)
    # Plot rawThresh as line
    plt.plot([0, len(wheelDiffs)], [wheelThresh, wheelThresh], 'r-')
    plt.plot([0, len(wheelDiffs)], [-wheelThresh, -wheelThresh], 'r-')
    plt.title('wheelDiffs')
    plt.show()
    # Plot cueDiffs
    plt.figure()
    plt.plot(cueDiffs)
    # Plot cueThresh as line
    plt.plot([0, len(cueDiffs)], [cueThresh, cueThresh], 'r-')
    plt.plot([0, len(cueDiffs)], [-cueThresh, -cueThresh], 'r-')
    plt.title('cueDiffs')
    plt.show()
    # Plot distractorDiffs
    plt.figure()
    plt.plot(distractorDiffs)
    # Plot distractorThresh as line
    plt.plot([0, len(distractorDiffs)], [distractorThresh, distractorThresh], 'r-')
    plt.plot([0, len(distractorDiffs)], [-distractorThresh, -distractorThresh], 'r-')
    plt.title('distractorDiffs')
    plt.show()
    # Plot reachDiffs
    plt.figure()
    plt.plot(reachDiffs)
    # Plot rawThresh as line
    plt.plot([0, len(reachDiffs)], [reachThresh, reachThresh], 'r-')
    plt.plot([0, len(reachDiffs)], [-reachThresh, -reachThresh], 'r-')
    plt.title('reachThresh')
    plt.show()
    # Plot optoDiffs
    plt.figure()
    plt.plot(optoDiffs)
    # Plot optoThresh as line
    plt.plot([0, len(optoDiffs)], [optoThresh, optoThresh], 'r-')
    plt.plot([0, len(optoDiffs)], [-optoThresh, -optoThresh], 'r-')
    plt.title('optoDiffs')
    plt.show()
    # once have tested a video so know the right thresholds, can set getDifferenceEventsBeforeSaving to True
    # in order to save smaller files for Matlab
    if getDifferenceEventsBeforeSaving:
        def condition(x, thresh): return x > thresh
        def condminus(x, thresh): return x < -thresh
        rawDiffEvs_plus = [idx for idx, element in enumerate(rawDiffs) if condition(element, rawThresh)]
        wheelDiffEvs_plus = [idx for idx, element in enumerate(wheelDiffs) if condition(element, wheelThresh)]
        cueDiffEvs_plus = [idx for idx, element in enumerate(cueDiffs) if condition(element, cueThresh)]
        distractorDiffEvs_plus = [idx for idx, element in enumerate(distractorDiffs) if condition(element, distractorThresh)]
        reachDiffEvs_plus = [idx for idx, element in enumerate(reachDiffs) if condition(element, reachThresh)]
        optoDiffEvs_plus = [idx for idx, element in enumerate(optoDiffs) if condition(element, optoThresh)]
        rawDiffEvs_minus = [idx for idx, element in enumerate(rawDiffs) if condminus(element, rawThresh)]
        wheelDiffEvs_minus = [idx for idx, element in enumerate(wheelDiffs) if condminus(element, wheelThresh)]
        cueDiffEvs_minus = [idx for idx, element in enumerate(cueDiffs) if condminus(element, cueThresh)]
        distractorDiffEvs_minus = [idx for idx, element in enumerate(distractorDiffs) if condminus(element, distractorThresh)]
        reachDiffEvs_minus = [idx for idx, element in enumerate(reachDiffs) if condminus(element, reachThresh)]
        optoDiffEvs_minus = [idx for idx, element in enumerate(optoDiffs) if condminus(element, optoThresh)]
        howmanyframes = len(cueDiffs)
        # save event indices
        outp = {"rawDiffEvs_plus": rawDiffEvs_plus, "cueDiffEvs_plus": cueDiffEvs_plus, "distractorDiffEvs_plus": distractorDiffEvs_plus, "wheelDiffEvs_plus": wheelDiffEvs_plus, "reachDiffEvs_plus": reachDiffEvs_plus, "optoDiffEvs_plus": optoDiffEvs_plus,
                "rawDiffEvs_minus": rawDiffEvs_minus, "cueDiffEvs_minus": cueDiffEvs_minus, "distractorDiffEvs_minus": distractorDiffEvs_minus, "wheelDiffEvs_minus": wheelDiffEvs_minus, "reachDiffEvs_minus": reachDiffEvs_minus, "optoDiffEvs_minus": optoDiffEvs_minus, "howmanyframes": howmanyframes,
                "vidStarts": vidStarts, "vidEnds": vidEnds}
    else:
        # save results
        outp = {"rawDiffs": rawDiffs, "cueDiffs": cueDiffs, "distractorDiffs": distractorDiffs, "wheelDiffs": wheelDiffs, "reachDiffs": reachDiffs, "vidStarts": vidStarts, "vidEnds": vidEnds}
    savemat("rig_events.mat", outp)
    # Print location of output file
    print("Saved rig_events.mat to ", os.getcwd())
    os.chdir(codedir)
    return nframes, rawDiffs, cueDiffs, distractorDiffs, wheelDiffs
            

def plotZonesOnExampleFrame(video, wheelZone, cueZone, distractorZone, reachZone, optoZone):
    # Load example frame
    cap = cv2.VideoCapture(video)
    ret, frame = cap.read()
    # Plot wheel zone
    plt.figure()
    plt.imshow(frame)
    plt.plot([wheelZone[2], wheelZone[2]], [wheelZone[0], wheelZone[1]], 'r-')
    plt.plot([wheelZone[3], wheelZone[3]], [wheelZone[0], wheelZone[1]], 'r-')
    plt.plot([wheelZone[2], wheelZone[3]], [wheelZone[0], wheelZone[0]], 'r-')
    plt.plot([wheelZone[2], wheelZone[3]], [wheelZone[1], wheelZone[1]], 'r-')
    plt.title('wheelZone')
    plt.show()
    # Plot cue zone
    plt.figure()
    plt.imshow(frame)
    plt.plot([cueZone[2], cueZone[2]], [cueZone[0], cueZone[1]], 'r-')
    plt.plot([cueZone[3], cueZone[3]], [cueZone[0], cueZone[1]], 'r-')
    plt.plot([cueZone[2], cueZone[3]], [cueZone[0], cueZone[0]], 'r-')
    plt.plot([cueZone[2], cueZone[3]], [cueZone[1], cueZone[1]], 'r-')
    plt.title('cueZone')
    plt.show()
    # Plot distractor zone
    plt.figure()
    plt.imshow(frame)
    plt.plot([distractorZone[2], distractorZone[2]], [distractorZone[0], distractorZone[1]], 'r-')
    plt.plot([distractorZone[3], distractorZone[3]], [distractorZone[0], distractorZone[1]], 'r-')
    plt.plot([distractorZone[2], distractorZone[3]], [distractorZone[0], distractorZone[0]], 'r-')
    plt.plot([distractorZone[2], distractorZone[3]], [distractorZone[1], distractorZone[1]], 'r-')
    plt.title('distractorZone')
    plt.show()
    # Plot reach zone
    plt.figure()
    plt.imshow(frame)
    plt.plot([reachZone[2], reachZone[2]], [reachZone[0], reachZone[1]], 'r-')
    plt.plot([reachZone[3], reachZone[3]], [reachZone[0], reachZone[1]], 'r-')
    plt.plot([reachZone[2], reachZone[3]], [reachZone[0], reachZone[0]], 'r-')
    plt.plot([reachZone[2], reachZone[3]], [reachZone[1], reachZone[1]], 'r-')
    plt.title('reachZone')
    plt.show()
    # Plot opto zone
    plt.figure()
    plt.imshow(frame)
    plt.plot([optoZone[2], optoZone[2]], [optoZone[0], optoZone[1]], 'r-')
    plt.plot([optoZone[3], optoZone[3]], [optoZone[0], optoZone[1]], 'r-')
    plt.plot([optoZone[2], optoZone[3]], [optoZone[0], optoZone[0]], 'r-')
    plt.plot([optoZone[2], optoZone[3]], [optoZone[1], optoZone[1]], 'r-')
    plt.title('optoZone')
    plt.show()
    # Close video
    cap.release()


def rgb2gray(rgb):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
    return gray


getRigEvents()