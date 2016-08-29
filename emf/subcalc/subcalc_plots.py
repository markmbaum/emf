import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import subcalc_funks

#rcparams for more static global formatting changes
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.figsize'] = (12, 6)
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.color'] = (.2, .2, .2)
mpl.rcParams['axes.labelcolor'] = (.2, .2, .2)
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.borderaxespad'] = 0 #mpl default is None
mpl.rcParams['xtick.color'] = (.2, .2, .2)
mpl.rcParams['ytick.color'] = (.2, .2, .2)

def contour_plot(mod):
    """Generate a contour plot from the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kwargs:

    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        CS - matplotlib QuadContourSet object"""

    #generate the figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #adjust axis width
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
    #plot contours
    CS = ax.contour(mod.grid_x, mod.grid_y, mod.grid_Bmax)
    #plot max location
    idx = np.argmax(mod.Bmax)
    peak_B = mod.Bmax[idx]
    h, = ax.plot(mod.x[idx], mod.y[idx], 'r*')
    #legend
    ax.legend(CS.collections + [h],
            [str(subcalc_funks._check_intable(i)) + ' mG'
                for i in CS.cvalues]
                + ['Peak Magnetic Field\n(%g mG)' % peak_B],
            loc = 'center left', bbox_to_anchor = (1.01, 0.5),
            fontsize = 12, numpoints = 1)
    #text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    return(fig, ax, CS)

def ion():
    plt.ion()

def show():
    plt.show()
