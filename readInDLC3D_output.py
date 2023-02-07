import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt 

def readInDLC3D_output(
):

    """
    Reads in output of DLC after 3D triangulate
    """

    directory=r'C:\Users\sabatini\Documents\DLC_3D\Test2022-Kim-2022-01-04-3d\test pipeline'
    filename=r'2019-05-30-152018-0000_DLC_3D.h5'

    df = pd.read_hdf(os.path.join(directory, filename))
    # Display all columns
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    print(df.head())

    bodyparts_under = ['palm','thumb_tip','wrist']

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

readInDLC3D_output()