import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt

def plotDLC3d(
):
    """
    Plot output of DLC 3D triangulation
    """

    directory=r'C:\Users\sabatini\Documents\DLC_3D\Test2022-Kim-2022-01-04-3d\test pipeline'
    filename='2019-05-30-152018-0001_DLC_resnet50_Testing2DJan4shuffle1_500000.h5'
    df = pd.read_hdf(os.path.join(directory, filename))
    # Print  columns of df
    print(df)

    # Check if any of df are not nan
    print(any(df.notna()))

    # Get rows of df that are not nan
    df = df[df.notna()]
    bodyparts = ['palm','thumb_tip','wrist']
    # For each bodypart
    for b in bodyparts:
        # Find column in level 1 of df that matches bodypart
        for c in df.columns.get_level_values(1):
            if c == b:
                x=df.loc[:,(slice(None),c,[d == 'x' for d in df.columns.get_level_values(2)])].xs('x', axis=1, level=2)
                y=df.loc[:,(slice(None),c,[d == 'y' for d in df.columns.get_level_values(2)])].xs('y', axis=1, level=2)
                z=df.loc[:,(slice(None),c,[d == 'z' for d in df.columns.get_level_values(2)])].xs('z', axis=1, level=2)
                # Get index of x
                coor = x.index
                # Convert coor to np array
                coor = coor.to_numpy()
                # x to numpy
                x = x.to_numpy()
                y = y.to_numpy()
                z = z.to_numpy()
                # Show size of x
                print(x.shape)
                # Which x are not nan
                nonan = ~np.isnan(x)
                # Make nonan 1d
                nonan = nonan.flatten()
                print(any(nonan == True))
                # Get which indices of nonan are True
                x = x[nonan == True]
                y = y[nonan == True]
                z = z[nonan == True]
                # Get only coor for nonan is True
                coor = coor[nonan == True]
                # Plot x
                print(x)
                plt.plot(x)
                plt.show()
                # 3D scatter plot color coded by coor
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(x, y, z, c=coor, marker='o')
                ax.set_xlabel('X Label')
                ax.set_ylabel('Y Label')
                ax.set_zlabel('Z Label')
                # Make figure title b
                plt.title(b)
                plt.show()
                return
                
                
                

plotDLC3d()