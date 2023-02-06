import os
import pandas as pd
from scipy.io import savemat

def plotDLC3d(
):
    """
    Plot output of DLC 3D triangulation
    """

    directory=r'C:\Users\sabatini\Documents\DLC_3D\Test2022-Kim-2022-01-04-3d\test pipeline'
    filename='2019-05-30-152018-0010_DLC_resnet50_Testing2DJan4shuffle1_500000.h5'
    df = pd.read_hdf(os.path.join(directory, filename))
    print('got here')
    # Print  columns of df
    print(df)

plotDLC3d()