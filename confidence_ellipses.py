import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipses():
    print("Plotting confidence ellipses")

    currcolor='red'
    
    # Read in .mat data
    vals_x_axis = scipy.io.loadmat(r'C:\Users\sabatini\Documents\confidence ellipse\vals_x_axis.mat')
    vals_y_axis = scipy.io.loadmat(r'C:\Users\sabatini\Documents\confidence ellipse\vals_y_axis.mat')

    # Discard if nan in vals_x_axis or vals_y_axis
    vals_x_axis = vals_x_axis['vals_x_axis']
    vals_y_axis = vals_y_axis['vals_y_axis']
    tokeep = ~np.isnan(vals_x_axis) & ~np.isnan(vals_y_axis)
    vals_x_axis = vals_x_axis[tokeep]
    vals_y_axis = vals_y_axis[tokeep]

    # Plot scatter plot of vals_x_axis versus vals_y_axis
    plt.scatter(vals_x_axis, vals_y_axis, s=10, c=currcolor, alpha=0.2)
    # Get axes of scatter plot
    ax = plt.gca()
    plt.xlabel('vals_x_axis')
    plt.ylabel('vals_y_axis')
    
    #confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=2.0, edgecolor=currcolor)
    confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=1.5, facecolor=currcolor, zorder=0, alpha=0.1)
    confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=1, facecolor=currcolor, zorder=0, alpha=0.1)
    confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=0.5, facecolor=currcolor, zorder=0, alpha=0.1)
    confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=0.25, facecolor=currcolor, zorder=0, alpha=0.1)
    confidence_ellipse(vals_x_axis, vals_y_axis, ax, n_std=0.1, facecolor=currcolor, zorder=0, alpha=0.1)

    # Show plot
    plt.show()


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


confidence_ellipses()


