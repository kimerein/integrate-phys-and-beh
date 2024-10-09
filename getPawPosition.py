import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt 
import cv2
import statsmodels.api as sm
from itertools import combinations

def getPawPosition(
):
    """
    Read output of DLC
    """

    #bodyparts = ['currentPellet','middle_knuckle','middle_tip','palm','pinky_knuckle','pinky_tip','pointer_knuckle','pointer_tip','ring_knuckle','ring_tip','thumb_knuckle','thumb_tip','wrist']
    #bodyparts = ['palm','thumb_tip','wrist'] # this is what I used for paw position only
    bodyparts =['palm','wrist','middle_knuckle','middle_tip','pinky_knuckle','pinky_tip','pointer_knuckle','pointer_tip','ring_knuckle','ring_tip','thumb_knuckle','thumb_tip']
    # bodyparts = ['palm']
    # Add Bottom to all bodyparts
    # bodyparts_under = [b+'Bottom' for b in bodyparts]
    bodyparts_under = bodyparts
    # Add Side to all bodyparts
    # bodyparts_side = [b+'Side' for b in bodyparts]
    bodyparts_side = bodyparts
    p_cutoff=0.1 #0.04
    p_cutoff_fingers=0.01

    # For each .h5 file, get horizontal and vertical position of each bodypart
    #directory=r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20191127\April_short\DLC vids'
    directory=r'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20200302\Oct_0\DLC output'
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
                # Get 3D positions of other bodyparts
                positions_side = getBodyPartPositions(directory, filename, bodyparts_side, p_cutoff_fingers, 'Side')
                positions_under = getBodyPartPositions(directory, filename, bodyparts_under, p_cutoff_fingers, 'Bottom')
                positions3D = get3DBodyPartDataFromTwoViews(positions_side, positions_under)
                # Calculate summed pointwise distance of all pairs of the finger body parts
                knuckle_distances = calculateSummedPointwiseDistance(positions3D, ['middle_knuckle','pinky_knuckle','pointer_knuckle','ring_knuckle','thumb_knuckle'])
                tip_distances = calculateSummedPointwiseDistance(positions3D, ['middle_tip','pinky_tip','pointer_tip','ring_tip','thumb_tip'])
                all_distances = calculateSummedPointwiseDistance(positions3D, ['middle_knuckle','pinky_knuckle','pointer_knuckle','ring_knuckle','thumb_knuckle','middle_tip','pinky_tip','pointer_tip','ring_tip','thumb_tip'])
                # Calculate distance between knuckles and finger tips
                distance_knuckletotip = calculateSummedDistanceKnuckleTip(positions3D, ['middle','pinky','pointer','ring','thumb'])
                # Save these to a .mat file
                savemat(os.path.join(directory, filename[:-4]+'knuckle_distances'+'.mat'),{'knuckle_distances':knuckle_distances})
                savemat(os.path.join(directory, filename[:-4]+'tip_distances'+'.mat'),{'tip_distances':tip_distances})
                savemat(os.path.join(directory, filename[:-4]+'all_distances'+'.mat'),{'all_distances':all_distances})
                savemat(os.path.join(directory, filename[:-4]+'distance_knuckletotip'+'.mat'),{'distance_knuckletotip':distance_knuckletotip})
                
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


def calculateSummedDistanceKnuckleTip(bodypart_positions_3D, include_these_bodyparts):
    """
    Calculate the summed distance between body parts labeled '_knuckle' and '_tip' at each time point
    for all pairs of body parts passed in that have matching names followed by '_knuckle' or '_tip'.
    
    Parameters:
    - bodypart_positions_3D: A dictionary containing the 3D positions of each body part, stored as separate X, Y, Z vectors.
    - include_these_bodyparts: A list of body parts to include in the distance calculation (should have matching '_knuckle' and '_tip' parts).
    
    Returns:
    - A vector (array) representing the summed distance between '_knuckle' and '_tip' labeled body parts for each time point.
    """
    
    # Initialize a vector to accumulate the summed distance for each time point
    n_timepoints = len(next(iter(bodypart_positions_3D.values()))['X'])  # Assuming all body parts have the same number of time points
    total_distance = np.zeros(n_timepoints)

    # Variable to count the number of valid knuckle-tip pairs
    count_pairs = np.zeros(n_timepoints)
    
    # Loop through the provided list to find matching pairs of '_knuckle' and '_tip'
    for bodypart in include_these_bodyparts:
        knuckle_label = bodypart + '_knuckle'
        tip_label = bodypart + '_tip'
        
        # Check if both '_knuckle' and '_tip' labels are present in the dictionary
        if knuckle_label in bodypart_positions_3D and tip_label in bodypart_positions_3D:
            # Extract the X, Y, Z coordinates for both '_knuckle' and '_tip'
            knuckle_positions = np.column_stack((
                bodypart_positions_3D[knuckle_label]['X'],
                bodypart_positions_3D[knuckle_label]['Y'],
                bodypart_positions_3D[knuckle_label]['Z']
            ))
            
            tip_positions = np.column_stack((
                bodypart_positions_3D[tip_label]['X'],
                bodypart_positions_3D[tip_label]['Y'],
                bodypart_positions_3D[tip_label]['Z']
            ))
            
            # Calculate the pointwise distance between '_knuckle' and '_tip'
            distances = np.linalg.norm(knuckle_positions - tip_positions, axis=1)
            
            # Ignore NaN values by using np.nan_to_num (convert NaN to zero)
            distances = np.nan_to_num(distances, nan=0.0)
            
            # Add distances to the total distance vector
            total_distance += distances

            # Increment count only where distances are valid (not NaN)
            count_pairs += ~np.isnan(distances)
        else:
            print(f"Warning: Missing pair for body part '{bodypart}' with '_knuckle' or '_tip' label.")
    
    # Calculate the average distance
    # Avoid division by zero by using np.where to only divide where count_pairs > 0
    average_distance = np.where(count_pairs > 0, total_distance / count_pairs, np.nan)
    
    return average_distance


def calculateSummedPointwiseDistance(bodypart_positions_3D, include_these_bodyparts):
    """
    Calculate the summed pointwise distance of all pairs of specified body parts based on their 3D positions.
    
    Parameters:
    - bodypart_positions_3D: A dictionary containing the 3D positions of each body part.
    - include_these_bodyparts: A list of body parts to include in the distance calculation.
    
    Returns:
    - A scalar value representing the summed pointwise distance of all specified pairs of body parts.
    """
    
    # Filter the body part positions to include only the specified body parts
    filtered_positions = {part: pos for part, pos in bodypart_positions_3D.items() if part in include_these_bodyparts}
    
    # Get the list of filtered body parts
    filtered_bodyparts = list(filtered_positions.keys())
    
    # Initialize a vector to accumulate the summed distance for each time point
    n_timepoints = len(next(iter(filtered_positions.values()))['X'])  # Assuming all body parts have the same number of time points
    total_distance = np.zeros(n_timepoints)
    
    # Generate all unique pairs of filtered body parts using combinations
    for part1, part2 in combinations(filtered_bodyparts, 2):
        # Extract the X, Y, Z coordinates for each pair and combine them into a 3D array
        positions1 = np.column_stack((
            filtered_positions[part1]['X'],
            filtered_positions[part1]['Y'],
            filtered_positions[part1]['Z']
        ))
        
        positions2 = np.column_stack((
            filtered_positions[part2]['X'],
            filtered_positions[part2]['Y'],
            filtered_positions[part2]['Z']
        ))
        
        # Calculate the pointwise distance between the two body parts
        # Assuming positions are stored as arrays with shape (n_timepoints, 3)
        distances = np.linalg.norm(positions1 - positions2, axis=1)
        
        # Ignore NaN values by using np.nan_to_num (convert NaN to zero)
        distances = np.nan_to_num(distances, nan=0.0)
        
        # Add distances to the total distance vector
        total_distance += distances
    
    return total_distance


def get3DBodyPartDataFromTwoViews(bodypart_positions_side, bodypart_positions_under):
    """
    Process the position data for each body part from two different views (side and under).
    
    Parameters:
    - bodypart_positions_side: A dictionary containing position and change data for each body part from the side view.
    - bodypart_positions_under: A dictionary containing position and change data for each body part from the under view.
    
    The function performs operations on corresponding data for each body part in both views.
    """
    
    # Get the set of body parts from both dictionaries
    bodyparts_side = set(bodypart_positions_side.keys())
    bodyparts_under = set(bodypart_positions_under.keys())
    
    # Find common body parts in both views
    common_bodyparts = bodyparts_side.intersection(bodyparts_under)

    # Define dictionary bodypart_positions_3D
    bodypart_positions_3D = {}
    
    # Iterate through each common body part in the dictionaries
    for bodypart in common_bodyparts:
        print(f"Processing data for body part: {bodypart}")
        
        # Access x and y positions for the side view
        x_positions_side = bodypart_positions_side[bodypart]['x']
        y_positions_side = bodypart_positions_side[bodypart]['y']
        
        # Access x and y positions for the under view
        x_positions_under = bodypart_positions_under[bodypart]['x']
        y_positions_under = bodypart_positions_under[bodypart]['y']
        
        # Access changes in x and y positions for both views
        #x_changes_side = bodypart_positions_side[bodypart]['x_change']
        #y_changes_side = bodypart_positions_side[bodypart]['y_change']
        #x_changes_under = bodypart_positions_under[bodypart]['x_change']
        #y_changes_under = bodypart_positions_under[bodypart]['y_change']

        # Convert 2D positions to 3D positions
        X, Y, Z, X_from_under = twoDtothreeD(x_positions_side,y_positions_side,x_positions_under,y_positions_under)
        
        # Return the 3D positions as a dictionary of all the common_bodyparts
        bodypart_positions_3D[bodypart] = {    'X': X, 'Y': Y, 'Z': Z, 'X_from_under': X_from_under    }
        
        print("\n")  # Add a newline for better readability between body parts

    return bodypart_positions_3D


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


def getBodyPartPositions(directory, filename, bodyparts, p_cutoff, substring_to_remove):
    """
    Get position data for multiple body parts from an HDF5 file.
    
    Parameters:
    - directory: The directory containing the HDF5 file.
    - filename: The name of the HDF5 file.
    - bodyparts: A list of body parts to find positions for.
    - p_cutoff: A likelihood cutoff below which positions are considered invalid.
    - substring_to_remove: The substring to be removed from the body part names.
    
    Returns:
    - A dictionary containing the x and y positions, and changes in these positions for each body part.
    """
    
    print('Now processing:')
    print(os.path.join(directory, filename))
    
    df = pd.read_hdf(os.path.join(directory, filename))
    print(df.head())
    
    # Initialize dictionaries to store the results
    bodypart_positions = {}
    
    # For each bodypart
    for b in bodyparts:
        # Find the column in level 1 of df that matches bodypart
        for c in df.columns.get_level_values(1):
            if c == b:
                # Remove the specified substring from the body part name
                bodypart_name_cleaned = c.replace(substring_to_remove, "")

                # Process x-coordinates
                temp_x = df.loc[:, (slice(None), c, [d == 'x' for d in df.columns.get_level_values(2)])].xs('x', axis=1, level=2)
                cf = df.loc[:, (slice(None), c, [d == 'likelihood' for d in df.columns.get_level_values(2)])].xs('likelihood', axis=1, level=2)
                temp_x = temp_x.to_numpy()
                temp_x[cf < p_cutoff] = np.nan
                temp_x = fillinwithnearest(temp_x)
                
                # Process y-coordinates
                temp_y = df.loc[:, (slice(None), c, [d == 'y' for d in df.columns.get_level_values(2)])].xs('y', axis=1, level=2)
                temp_y = temp_y.to_numpy()
                temp_y[cf < p_cutoff] = np.nan
                temp_y = fillinwithnearest(temp_y)
                
                # Calculate changes in positions
                temp_x_change = np.diff(temp_x, 1, 0)
                temp_y_change = np.diff(temp_y, 1, 0)
                
                # Store the results in the dictionary
                bodypart_positions[bodypart_name_cleaned] = {
                    'x': temp_x,
                    'y': temp_y,
                    'x_change': temp_x_change,
                    'y_change': temp_y_change
                }
    
    return bodypart_positions


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