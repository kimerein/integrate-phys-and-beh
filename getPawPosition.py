import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt 

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
    bodyparts_side = [b+'Side' for b in bodyparts]
    p_cutoff=0.04

    # For each .h5 file, get horizontal and vertical position of each bodypart
    directory=r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\curr'
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".h5"):
            print('Now processing ')
            print(os.path.join(directory, filename))
            df = pd.read_hdf(os.path.join(directory, filename))

            print(df.head())

            # For each bodypart
            for b in bodyparts_under:
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
                        if b == bodyparts_under[0]:
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
                        if b == bodyparts_under[0]:
                            pawPos_y = temp
                            # Change in neighboring elements of pawPos_y
                            pawPos_y_change = np.diff(temp,1,0)
                        else:
                            pawPos_y = np.nansum(np.dstack((pawPos_y,temp)),axis=2)
                            pawPos_y_change = np.nansum(np.dstack((pawPos_y_change,np.diff(temp,1,0))),axis=2)

            pawPos_x = pawPos_x/len(bodyparts_under)
            pawPos_y = pawPos_y/len(bodyparts_under)
            pawPos_x_change = pawPos_x_change/len(bodyparts_under)
            pawPos_y_change = pawPos_y_change/len(bodyparts_under)

            # Integrate pawPos_x_change and pawPos_y_change to get pawPos_x and pawPos_y
            # pawPos_x = np.cumsum(pawPos_x_change,0)
            # pawPos_y = np.cumsum(pawPos_y_change,0)
            
            # Plot paw position x and y with color changing as function of index
            plt.figure()
            # plt.scatter(pawPos_x_change,pawPos_y_change,c=np.arange(len(pawPos_x_change)))
            plt.scatter(pawPos_x,pawPos_y,c=np.arange(len(pawPos_x)))
            # Drop nan from pawPos_x and pawPos_y
            # pawPos_x = pawPos_x[~np.isnan(pawPos_x)]
            # pawPos_y = pawPos_y[~np.isnan(pawPos_y)]
            pawPos_x_change = pawPos_x_change[~np.isnan(pawPos_x_change)]
            pawPos_y_change = pawPos_y_change[~np.isnan(pawPos_y_change)]
            # Drop pawPos_x and pawPos_y if they are both less than 100
            # Plot line connecting points, but skip points where either x or y is less than 100
            # plt.plot(pawPos_x[pawPos_x>100],pawPos_y[pawPos_x>100])
            plt.plot(pawPos_x,pawPos_y)
            # plt.plot(pawPos_x_change,pawPos_y_change)
            # Set x limits
            plt.xlim(0, 350)
            plt.ylim(0, 350)
            # plt.xlim(-500, 500)
            # plt.ylim(-500, 500)
            plt.colorbar()
            plt.show()


        else:
            continue
    

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

getPawPosition()