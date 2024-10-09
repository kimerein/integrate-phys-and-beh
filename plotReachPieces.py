import os
import numpy as np
import cv2
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy.io 
#from deeplabcut.utils import auxiliaryfunctions
from pathlib import Path
import glob

def plotReachPieces(
):
    """
    Analyze epochs of reach during opto or no opto.
    Output: Plots the average of X, Y, and Z values around positive deflections in X separately when optoOn is True or False.
    """

    #videos = r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20191127\April_short\DLC vids'
    videos = r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20200302\Oct_0\DLC output'
    videotype = ".avi"
    # User set these values
    # Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190510\3F_white\DLC output
    optoZone = [217, 234, 315, 336] # [x_start, x_end, y_start, y_end]
    optoThresh = 150  # Threshold for optoOn
    window_size = 200  # Window around the peak (e.g., 40 frames before and after)
    Rise = 20  # Minimum rise in X to be considered a deflection
    RiseTimeInFrames = 80  # Maximum time in frames for the rise to occur
    showThisManyFiles=5 # Show the first this many files for debugging

    ##################################################
    # Looping over videos
    ##################################################
    print("Finding videos")
    Videos = []
    for filename in os.listdir(videos):
        if filename.endswith(videotype):
            Videos.append(filename)
    Videos = sorted(Videos)
    print("Found ", len(Videos), " videos")

    # Initialize arrays to store the average X, Y, and Z values for optoOn True and False with NaN
    avg_X_true = np.full(2 * window_size, np.nan)
    avg_Y_true = np.full(2 * window_size, np.nan)
    avg_Z_true = np.full(2 * window_size, np.nan)
    avg_tip_true = np.full(2 * window_size, np.nan)
    avg_X_false = np.full(2 * window_size, np.nan)
    avg_Y_false = np.full(2 * window_size, np.nan)
    avg_Z_false = np.full(2 * window_size, np.nan)
    avg_tip_false = np.full(2 * window_size, np.nan)
    num_deflections_true = 0
    num_deflections_false = 0
    all_X_true = []
    all_Y_true = []
    all_Z_true = []
    all_tip_true = []
    all_X_false = []
    all_Y_false = []
    all_Z_false = []
    all_tip_false = []
    
    if len(Videos) > 0:
        # Go to location of videos
        os.chdir(videos)
        for idx, video in enumerate(Videos):
            ##################################################
            # Show zones on example frame
            ##################################################
            # If this is the first video
            if idx == 0:
                plotZonesOnExampleFrame(video, optoZone)

            ##################################################
            # Loading the video
            ##################################################
            print("Starting to analyze ", video)
            print("Loading ", video)
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
            optoOn = []
            optoVals = []
            while cap.isOpened():
                if counter % step == 0:
                    pbar.update(step)

                ret, frame = cap.read()
                if ret:
                    ny, nx, nc = np.shape(frame)
                    frame = rgb2gray(frame)
                    ny, nx = np.shape(frame)

                    # Append to optoOn whether the value in optoZone is above threshold optoThresh
                    optoOn.append(np.mean(frame[optoZone[0]:optoZone[1], optoZone[2]:optoZone[3]]) > optoThresh)
                    optoVals.append(np.mean(frame[optoZone[0]:optoZone[1], optoZone[2]:optoZone[3]]))
                elif counter >= nframes:
                    break
                counter += 1
            pbar.close()
            
            # Plot optoVals and threshold of example video if this is first video
            # if video is one of the first 10 videos
            if idx < showThisManyFiles:
                plt.figure()
                plt.plot(optoVals, label='optoVals')
                plt.axhline(optoThresh, color='r', linestyle='--', label='Threshold')
                plt.title("optoVals")
                plt.show()

             # Use glob to find a .mat file that matches the wildcard pattern
            mat_pattern = video.split(".")[0] + "*" + "3D.mat"
            mat_files = glob.glob(mat_pattern)
            tip_pattern = video.split(".")[0] + "*" + "tip_distances.mat"
            tip_files = glob.glob(tip_pattern)
            print(f"Found {len(mat_files)} .mat files matching the pattern {mat_pattern}")
            print(mat_files)

            if len(mat_files) > 0:
                # Load the first .mat file that matches the pattern
                mat_filename = mat_files[0]
                print(f"Loading .mat file: {mat_filename}")
                mat_data = scipy.io.loadmat(mat_filename)
                tip_distances = scipy.io.loadmat(tip_files[0])['tip_distances']

                # Assuming 'X' is the key in the .mat file containing the vector of values for each frame
                X = mat_data['X'].flatten()  # Ensure X is a 1D vector
                Y = mat_data['Y'].flatten()  # Ensure Y is a 1D vector
                Z = mat_data['Z'].flatten()  # Ensure Z is a 1D vector
                tips = tip_distances.flatten()  # Ensure tips is a 1D vector
                optoOn_array = np.array(optoOn)  # Convert optoOn to numpy array for easier indexing

                # Plot optoOn_array 
                if idx < showThisManyFiles:
                    plt.figure()
                    plt.plot(optoOn_array, label='optoOn')
                    plt.title("optoOn")
                    plt.show()

                # Plot X
                #plt.figure()
                #plt.plot(X, label='X')
                #plt.title("X")
                #plt.show()

                 # Find positive-going deflections in X when optoOn is True and False
                optoOn_true_deflections = []
                optoOn_false_deflections = []
                copyX = X.copy()
                for i in range(len(X) - RiseTimeInFrames):
                    # Check for a positive change greater than Rise in less than RiseTimeInFrames frames
                    if copyX[i + RiseTimeInFrames] - copyX[i] > Rise:
                        #print("Found deflection at frame", i + np.argmax(X[i:i + RiseTimeInFrames]))
                        peak_idx = i + np.argmax(copyX[i:i + RiseTimeInFrames])  # Find peak in the deflection window
                        #print("optoOn: ", optoOn_array[i + RiseTimeInFrames])
                        # If optoOn_array was true for any time in the deflection window, add this deflection to the list
                        if np.any(optoOn_array[i:i + RiseTimeInFrames]):
                            optoOn_true_deflections.append(peak_idx)
                        else:
                            optoOn_false_deflections.append(peak_idx)
                        # Set all these timepoints in X to NaN to avoid double counting
                        copyX[i:i + RiseTimeInFrames] = np.nan
                        if idx < showThisManyFiles:
                            plt.figure()
                            plt.plot(X[peak_idx - window_size:peak_idx + window_size], label='X')
                            plt.axvline(window_size, color='r', linestyle='--', label='Peak')
                            plt.title("Deflection")
                            plt.show()

                # Loop through the detected deflections for optoOn == True
                for peak_i in optoOn_true_deflections:
                    avg_X_true, avg_Y_true, avg_Z_true, avg_tip_true, num_deflections_true = update_sum(peak_i, X, Y, Z, tips, avg_X_true, avg_Y_true, avg_Z_true, avg_tip_true, window_size, num_deflections_true)
                    all_X_true, all_Y_true, all_Z_true, all_tip_true = update_append(peak_i, X, Y, Z, tips, all_X_true, all_Y_true, all_Z_true, all_tip_true, window_size)

                # Loop through the detected deflections for optoOn == False
                for peak_i in optoOn_false_deflections:
                    avg_X_false, avg_Y_false, avg_Z_false, avg_tip_false, num_deflections_false = update_sum(peak_i, X, Y, Z, tips, avg_X_false, avg_Y_false, avg_Z_false, avg_tip_false, window_size, num_deflections_false)
                    all_X_false, all_Y_false, all_Z_false, all_tip_false = update_append(peak_i, X, Y, Z, tips, all_X_false, all_Y_false, all_Z_false, all_tip_false, window_size)
                    # Plot current avg_X_false
                    #plt.figure()
                    #plt.plot(avg_X_false, label='X (optoOn=False)')
                    #plt.axvline(window_size, color='r', linestyle='--', label='Peak (optoOn=False)')
                    #plt.legend()
                    #plt.title(f"Averaged X around positive deflections (optoOn=False, N={num_deflections_false})")
                    #plt.xlabel("Frames (relative to peak)")
                    #plt.ylabel("Position")
                    #plt.show()
                
    print("num_deflections_true: ", num_deflections_true)
    print("num_deflections_false: ", num_deflections_false)

    # Divide through by num_deflections to get the average
    avg_X_true /= num_deflections_true
    avg_Y_true /= num_deflections_true
    avg_Z_true /= num_deflections_true
    avg_tip_true /= num_deflections_true
    avg_X_false /= num_deflections_false
    avg_Y_false /= num_deflections_false
    avg_Z_false /= num_deflections_false
    avg_tip_false /= num_deflections_false
    
    # Plot the averaged X, Y, and Z for optoOn == True
    plotAverageAroundDeflections(avg_X_true, True, window_size, num_deflections_true, "X")
    plotAverageAroundDeflections(avg_Y_true, True, window_size, num_deflections_true, "Y")
    plotAverageAroundDeflections(avg_Z_true, True, window_size, num_deflections_true, "Z")
    plotAverageAroundDeflections(avg_tip_true, True, window_size, num_deflections_true, "tip")

    # Plot the averaged X, Y, and Z for optoOn == False
    plotAverageAroundDeflections(avg_X_false, False, window_size, num_deflections_false, "X")
    plotAverageAroundDeflections(avg_Y_false, False, window_size, num_deflections_false, "Y")
    plotAverageAroundDeflections(avg_Z_false, False, window_size, num_deflections_false, "Z")
    plotAverageAroundDeflections(avg_tip_false, False, window_size, num_deflections_false, "tip")

    # Save to .mat files
    print("Saving to .mat files")
    scipy.io.savemat("avg_X_true.mat", {"avg_X_true": avg_X_true})
    scipy.io.savemat("avg_Y_true.mat", {"avg_Y_true": avg_Y_true})
    scipy.io.savemat("avg_Z_true.mat", {"avg_Z_true": avg_Z_true})
    scipy.io.savemat("avg_tip_true.mat", {"avg_tip_true": avg_tip_true})
    scipy.io.savemat("avg_X_false.mat", {"avg_X_false": avg_X_false})
    scipy.io.savemat("avg_Y_false.mat", {"avg_Y_false": avg_Y_false})
    scipy.io.savemat("avg_Z_false.mat", {"avg_Z_false": avg_Z_false})
    scipy.io.savemat("avg_tip_false.mat", {"avg_tip_false": avg_tip_false})
    scipy.io.savemat("all_X_true.mat", {"all_X_true": all_X_true})
    scipy.io.savemat("all_Y_true.mat", {"all_Y_true": all_Y_true})
    scipy.io.savemat("all_Z_true.mat", {"all_Z_true": all_Z_true})
    scipy.io.savemat("all_tip_true.mat", {"all_tip_true": all_tip_true})
    scipy.io.savemat("all_X_false.mat", {"all_X_false": all_X_false})
    scipy.io.savemat("all_Y_false.mat", {"all_Y_false": all_Y_false})
    scipy.io.savemat("all_Z_false.mat", {"all_Z_false": all_Z_false})
    scipy.io.savemat("all_tip_false.mat", {"all_tip_false": all_tip_false})
    print("Done!")



def update_sum(peak_idx, X, Y, Z, tip, avg_X, avg_Y, avg_Z, avg_tip, window_size, num_deflections):
    
    if peak_idx - window_size >= 0 and peak_idx + window_size < len(X):
        x_segment = X[peak_idx - window_size:peak_idx + window_size]
        y_segment = Y[peak_idx - window_size:peak_idx + window_size]
        z_segment = Z[peak_idx - window_size:peak_idx + window_size]
        tip_segment = tip[peak_idx - window_size:peak_idx + window_size]
                        
        # Average ignoring NaNs
        avg_X = np.nansum([avg_X, x_segment], axis=0)
        avg_Y = np.nansum([avg_Y, y_segment], axis=0)
        avg_Z = np.nansum([avg_Z, z_segment], axis=0)
        avg_tip = np.nansum([avg_tip, tip_segment], axis=0)

        num_deflections += 1
    return avg_X, avg_Y, avg_Z, avg_tip, num_deflections



def update_append(peak_idx, X, Y, Z, tip, all_X, all_Y, all_Z, all_tip, window_size):
    """
    Appends each x_segment, y_segment, z_segment, and tip_segment to growing arrays.
    
    Parameters:
    - peak_idx: The index of the detected peak.
    - X, Y, Z, tip: The arrays containing X, Y, Z, and tip data.
    - all_X, all_Y, all_Z, all_tip: Arrays where segments will be appended.
    - window_size: Number of frames before and after the peak to include.
    
    Returns:
    - all_X, all_Y, all_Z, all_tip: Arrays with appended segments.
    """
    if peak_idx - window_size >= 0 and peak_idx + window_size < len(X):
        x_segment = X[peak_idx - window_size:peak_idx + window_size]
        y_segment = Y[peak_idx - window_size:peak_idx + window_size]
        z_segment = Z[peak_idx - window_size:peak_idx + window_size]
        tip_segment = tip[peak_idx - window_size:peak_idx + window_size]
        
        # Append the segments to the growing lists
        all_X.append(x_segment)
        all_Y.append(y_segment)
        all_Z.append(z_segment)
        all_tip.append(tip_segment)

    return all_X, all_Y, all_Z, all_tip



def plotAverageAroundDeflections(av_X, optoOn_trueOrFalse, window_size, N, name):
    plt.figure()
    plt.plot(av_X, label=f"{name} (optoOn={optoOn_trueOrFalse})")
    plt.axvline(window_size, color='r', linestyle='--', label=f"Peak")
    plt.legend()
    plt.title(f"Averaged {name} around positive deflections (optoOn={optoOn_trueOrFalse}, N={N})")
    plt.xlabel("Frames (relative to peak)")
    plt.ylabel("Position")
    plt.show()
    

def plotZonesOnExampleFrame(video, optoZone):
    # Load example frame
    cap = cv2.VideoCapture(video)
    ret, frame = cap.read()
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
    
    
plotReachPieces()