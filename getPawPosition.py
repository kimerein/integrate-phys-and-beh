import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt 
import cv2
import statsmodels.api as sm

def getPawPosition(
):
    """
    Read output of DLC
    """

    #bodyparts = ['currentPellet','middle_knuckle','middle_tip','palm','pinky_knuckle','pinky_tip','pointer_knuckle','pointer_tip','ring_knuckle','ring_tip','thumb_knuckle','thumb_tip','wrist']
    bodyparts = ['palm','thumb_tip','wrist']
    # bodyparts = ['palm']
    # Add Bottom to all bodyparts
    # bodyparts_under = [b+'Bottom' for b in bodyparts]
    bodyparts_under = bodyparts
    # Add Side to all bodyparts
    # bodyparts_side = [b+'Side' for b in bodyparts]
    bodyparts_side = bodyparts
    p_cutoff=0.04

    # For each .h5 file, get horizontal and vertical position of each bodypart
    directory=r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190530\March_C\DLC vids'
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".h5"):
            # If filename contains side, then use bodyparts_side
            if 'side' in filename:
                pp = getPawPosData(directory,filename,bodyparts_side,p_cutoff)
                pawPos_x = pp[0]
                pawPos_y = pp[1]
            if 'under' in filename:
                pp = getPawPosData(directory,filename,bodyparts_under,p_cutoff)
                pawPos_x = pp[0]
                pawPos_y = pp[1]
            # If neither under nor side, then use bodyparts, first side then under
            if ~('side' in filename) and ~('under' in filename):
                # Concatenate "Side" and "Bottom" to bodyparts
                bodyparts_side = [b+'Side' for b in bodyparts]
                bodyparts_under = [b+'Bottom' for b in bodyparts]
                pp = getPawPosData(directory,filename,bodyparts_side,p_cutoff)
                pawPos_x = pp[0]
                pawPos_y = pp[1]
                # Save to .mat file
                savemat(os.path.join(directory, filename[:-4]+'side'+'.mat'),{'pawPos_x':pawPos_x,'pawPos_y':pawPos_y})
                pp = getPawPosData(directory,filename,bodyparts_under,p_cutoff)
                pawPos_x_under = pp[0]
                pawPos_y_under = pp[1]
                savemat(os.path.join(directory, filename[:-4]+'under'+'.mat'),{'pawPos_x':pawPos_x_under,'pawPos_y':pawPos_y_under})
                X, Y, Z, X_from_under = twoDtothreeD(pawPos_x,pawPos_y,pawPos_x_under,pawPos_y_under)
                #X, Y, Z, X_from_under = filt3D(X,Y,Z,X_from_under)
                savemat(os.path.join(directory, filename[:-4]+'3D'+'.mat'),{'X':X,'Y':Y,'Z':Z,'X_from_under':X_from_under})
        elif filename.endswith(".avi"):
            # If there is no corresponding .h5 file with the same name, then create a .mat file with nan
            if os.path.exists(os.path.join(directory, filename[:-4]+'DLC_resnet50_Testing2DJan4shuffle1_500000.h5')):
                continue
            # Count how many frames in this .avi file
            cap = cv2.VideoCapture(os.path.join(directory, filename))
            length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
            cap.release()
            # Fill in pawPos_x and pawPos_y with nan
            pawPos_x = np.full((length,1),np.nan)
            pawPos_y = np.full((length,1),np.nan)
            # Save to .mat file
            savemat(os.path.join(directory, filename[:-4]+'side'+'.mat'),{'pawPos_x':pawPos_x,'pawPos_y':pawPos_y})
            savemat(os.path.join(directory, filename[:-4]+'under'+'.mat'),{'pawPos_x':pawPos_x,'pawPos_y':pawPos_y})
        else:
            continue


def filt3D(X,Y,Z,X_from_under):
    # Filter X, Y and Z
    X_filt = arimafilt(X)
    Y_filt = arimafilt(Y)
    Z_filt = arimafilt(Z)
    X_from_under_filt = arimafilt(X_from_under)
    return X_filt, Y_filt, Z_filt, X_from_under_filt


def arimafilt(X):
    model = sm.tsa.ARIMA(X, order=(3,1,10))
    results = model.fit()
    filtered_data = results.fittedvalues
    return filtered_data


def twoDtothreeD(side_pawPos_x,side_pawPos_y,under_pawPos_x,under_pawPos_y):
    # This uses my definitions of X, Y and Z
    Z = side_pawPos_x # minus is more dorsal, plus is more ventral
    # X is average of side_pawPos_y and under_pawPos_y
    X = side_pawPos_y # minus is more toward pellet, plus is more toward perch/body
    # Rotate under view to line up perch axis with Y axis and perch-to-pellet axis with X axis
    #cosine = np.cos(-np.pi/4)
    #sine = np.sin(-np.pi/4)
    #rotationMatrix = np.array([[cosine, -sine], [sine, cosine]]) # rotate counterclockwise by 45 degrees
    # For each element of under_pawPos_x and under_pawPos_y, rotate by 45 degrees
    #under_pawPos_x_rotated = np.empty_like(under_pawPos_x) * np.nan
    #under_pawPos_y_rotated = np.empty_like(under_pawPos_y) * np.nan
    #for i in range(len(under_pawPos_x)):
    #    under_pawPos_x_rotated[i],under_pawPos_y_rotated[i] = np.dot(rotationMatrix, np.array([under_pawPos_x[i],under_pawPos_y[i]]))
    #under_pawPos_x_rotated = - under_pawPos_x_rotated # minus is more toward pellet, plus is more toward perch/body
    #X_from_under = under_pawPos_x_rotated
    #Y = under_pawPos_y_rotated
    Y = under_pawPos_y
    X_from_under = under_pawPos_x
    return X,Y,Z,X_from_under


def fillinwithnearest(data):
    # Step through data and fill in nan values with nearest value
    # Find first value that is not nan in data
    # If all is nan, return
    if np.all(np.isnan(data)):
        return data
    first = np.where(~np.isnan(data))[0][0]
    # data at first
    defaultfill = data[first] 
    data[1] = defaultfill
    for i in range(len(data)):
        if np.isnan(data[i]):
            data[i] = data[i-1]
    return data


def getPawPosData(directory,filename,bodyparts,p_cutoff):
    # Get paw position data for these bodyparts from this file
    print('Now processing ')
    print(os.path.join(directory, filename))
    df = pd.read_hdf(os.path.join(directory, filename))

    print(df.head())

    # For each bodypart
    for b in bodyparts:
        # Find column in level 1 of df that matches bodypart
        for c in df.columns.get_level_values(1):
            if c == b:
                #print(df.loc[:,(slice(None),c,[d == 'x' for d in df.columns.get_level_values(2)])])
                temp=df.loc[:,(slice(None),c,[d == 'x' for d in df.columns.get_level_values(2)])].xs('x', axis=1, level=2)
                cf=df.loc[:,(slice(None),c,[d == 'likelihood' for d in df.columns.get_level_values(2)])].xs('likelihood', axis=1, level=2)
                # Convert temp to np array
                temp = temp.to_numpy()
                # Make temp nan if likelihood is less than p_cutoff
                temp[cf<p_cutoff]=np.nan
                temp = fillinwithnearest(temp)
                if b == bodyparts[0]:
                    pawPos_x = temp
                    # Change in neighboring elements of pawPos_x
                    pawPos_x_change = np.diff(temp,1,0)
                else:
                    # Sum pawPos_x with temp ignoring nan values
                    pawPos_x = np.nansum(np.dstack((pawPos_x,temp)),axis=2)
                    pawPos_x_change = np.nansum(np.dstack((pawPos_x_change,np.diff(temp,1,0))),axis=2)

                temp=df.loc[:,(slice(None),c,[d == 'y' for d in df.columns.get_level_values(2)])].xs('y', axis=1, level=2)
                # Convert temp to np array
                temp = temp.to_numpy()
                # Make temp nan if likelihood is less than p_cutoff
                temp[cf<p_cutoff]=np.nan
                # Fill in nan elements of temp with closest value
                temp = fillinwithnearest(temp)
                if b == bodyparts[0]:
                    pawPos_y = temp
                    # Change in neighboring elements of pawPos_y
                    pawPos_y_change = np.diff(temp,1,0)
                else:
                    pawPos_y = np.nansum(np.dstack((pawPos_y,temp)),axis=2)
                    pawPos_y_change = np.nansum(np.dstack((pawPos_y_change,np.diff(temp,1,0))),axis=2)

    pawPos_x = pawPos_x/len(bodyparts)
    pawPos_y = pawPos_y/len(bodyparts)
    pawPos_x_change = pawPos_x_change/len(bodyparts)
    pawPos_y_change = pawPos_y_change/len(bodyparts)

    # Integrate pawPos_x_change and pawPos_y_change to get pawPos_x and pawPos_y
    # pawPos_x = np.cumsum(pawPos_x_change,0)
    # pawPos_y = np.cumsum(pawPos_y_change,0)
    
    # # Plot paw position x and y with color changing as function of index
    # plt.figure()
    # # plt.scatter(pawPos_x_change,pawPos_y_change,c=np.arange(len(pawPos_x_change)))
    # plt.scatter(pawPos_x,pawPos_y,c=np.arange(len(pawPos_x)))
    # # Drop nan from pawPos_x and pawPos_y
    # # pawPos_x = pawPos_x[~np.isnan(pawPos_x)]
    # # pawPos_y = pawPos_y[~np.isnan(pawPos_y)]

    pawPos_x_change = pawPos_x_change[~np.isnan(pawPos_x_change)]
    pawPos_y_change = pawPos_y_change[~np.isnan(pawPos_y_change)]
    
    # # Drop pawPos_x and pawPos_y if they are both less than 100
    # # Plot line connecting points, but skip points where either x or y is less than 100
    # # plt.plot(pawPos_x[pawPos_x>100],pawPos_y[pawPos_x>100])
    # plt.plot(pawPos_x,pawPos_y)
    # # plt.plot(pawPos_x_change,pawPos_y_change)
    # # Set x limits
    # plt.xlim(0, 350)
    # plt.ylim(0, 350)

    # plt.xlim(0, 1500)
    # plt.ylim(0, 1500)

    # # plt.xlim(-500, 500)
    # # plt.ylim(-500, 500)
    # plt.colorbar()
    # plt.show()
    return pawPos_x, pawPos_y, pawPos_x_change, pawPos_y_change


getPawPosition()