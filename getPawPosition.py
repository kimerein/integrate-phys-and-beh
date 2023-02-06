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
    bodyparts = ['currentPellet','palm','thumb_tip','wrist']
    # Add Bottom to all bodyparts
    bodyparts_under = [b+'Bottom' for b in bodyparts]
    # Add Side to all bodyparts
    bodyparts_side = [b+'Side' for b in bodyparts]
    p_cutoff=0.08

    # For each .h5 file, get horizontal and vertical position of each bodypart
    directory=r'Z:\MICROSCOPE\Kim\for_orchestra\DLC_testvids\curr'
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".h5"):
            print('Now processing ')
            print(os.path.join(directory, filename))
            df = pd.read_hdf(os.path.join(directory, filename))

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
                        if b == bodyparts_under[0]:
                            pawPos_x = temp
                        else:
                            # Sum pawPos_x with temp ignoring nan values
                            pawPos_x = np.nansum(np.dstack((pawPos_x,temp)),axis=2)

                        temp=df.loc[:,(slice(None),c,[d == 'y' for d in df.columns.get_level_values(2)])].xs('y', axis=1, level=2)
                        # Convert temp to np array
                        temp = temp.to_numpy()
                        # Make temp nan if likelihood is less than p_cutoff
                        temp[cf<p_cutoff]=np.nan
                        if b == bodyparts_under[0]:
                            pawPos_y = temp
                        else:
                            pawPos_y = np.nansum(np.dstack((pawPos_y,temp)),axis=2)

            pawPos_x = pawPos_x/len(bodyparts_under)
            pawPos_y = pawPos_y/len(bodyparts_under)
            
            # Plot paw position x and y with color changing as function of index
            plt.figure()
            plt.scatter(pawPos_x,pawPos_y,c=np.arange(len(pawPos_x)))
            plt.colorbar()
            plt.show()


        else:
            continue
    
getPawPosition()