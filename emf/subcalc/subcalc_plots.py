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

#other more specific/dynamic global formatting variables
_subcalc_plots_footprint_alpha = .75
_subcalc_plots_footprint_zorder = -1
_subcalc_plots_field_marker_size = 8

def _label_footprint_group(ax, fps):
    """Put text into ax for the footprint group fps
    args:
        ax - Axes object to plot in
        fps - list of footprints"""
    if(len(fps) > 1):
        x = subcalc_funks._flatten([fp.x for fp in fps])
        y = subcalc_funks._flatten([fp.y for fp in fps])
        x = (min(x) + max(x))/2.0
        y = (min(y) + max(y))/2.0
        ax.text(x, y, fps[0].group, ha = 'center', va = 'center')
    elif(len(fps) == 1):
        fp = fps[0]
        x = fp.x
        y = fp.y
        if(fp.draw_as_loop):
            x = (min(x) + max(x))/2.0
            y = (min(y) + max(y))/2.0
        else:
             x = x[0]
             y = y[0]
        ax.text(x, y, fp.group)

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

    #assess dimensions of grid
    xmax, xmin = np.max(mod.x), np.min(mod.x)
    ymax, ymin = np.max(mod.y), np.min(mod.y)
    aspect_ratio = (ymax - ymin)/(xmax - xmin)
    xmarg, ymarg = (xmax - xmin)*0.005, (ymax - ymin)*0.005
    #generate the figure
    fig = plt.figure(figsize = (12 + 4, 12*aspect_ratio))
    ax = fig.add_subplot(1,1,1)
    #adjust axis dimensions by transforming between figure and display coords
    box = ax.get_position()
    T = fig.transFigure.transform
    I = fig.transFigure.inverted().transform
    width = I((T((0., box.height))[1]/aspect_ratio, 0.))[0]
    ax.set_position([box.x0, box.y0, width, box.height])
    #plot contours
    CS = ax.contour(mod.X, mod.Y, mod.B, lindwidth = 5,
            levels = [.1,.5,1.,5.,10.,25.,50.],
            locator = mpl.ticker.LogLocator(), zorder = -1)
    #plot location of maximum field
    handles, labels = [], []
    peak_B, yidx, xidx = subcalc_funks._2Dmax(mod.B)
    peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
    handles.append(ax.plot(mod.x[xidx], mod.y[yidx], 'ro',
            markersize = _subcalc_plots_field_marker_size)[0])
    labels.append('Maximum Modeled\nMagnetic Field')
    ax.text(mod.x[xidx] + xmarg, mod.y[yidx] + ymarg, peak_B)
    #plot footprints
    for g in mod.footprint_groups:
        #get footprints in the group and the group name
        fps = [mod.footprints[i] for i in g]
        group_name = fps[0].group
        #plot individual footprints
        #check for the power line
        power_line_check = [fp.power_line for fp in fps]
        if(True in power_line_check):
            idx = power_line_check.index(True)
            fp = fps[idx]
            handles.append(ax.plot(fp.x, fp.y, 'k--', linewidth = 2,
                    alpha = _subcalc_plots_footprint_alpha,
                    zorder = _subcalc_plots_footprint_zorder)[0])
            labels.append(fp.name)
            #remove the power line footprint
            fps.pop(idx)
        #plot the other footprints
        of_concern = False
        for fp in fps:
            ax.plot(fp.x, fp.y, 'k',
                    alpha = _subcalc_plots_footprint_alpha,
                    zorder = _subcalc_plots_footprint_zorder)
            #mark maximum field if the footprint is 'of concern'
            if(fp.of_concern):
                x = fp.x
                y = fp.y
                B_interp = mod.interp(x, y)
                idx = np.argmax(B_interp)
                m = str(subcalc_funks._sig_figs(B_interp[idx], 3))
                h = ax.plot(x[idx], y[idx], 'yo',
                        markersize = _subcalc_plots_field_marker_size)[0]
                ax.text(x[idx] + xmarg, y[idx] + ymarg, m)
                if(of_concern is False):
                    handles.append(h)
                    labels.append('Points of Potential\nConcern')
                    of_concern = True

        #deal with labeling
        _label_footprint_group(ax, fps)

    #legend
    cvalues = [str(subcalc_funks._check_intable(i)) + ' mG' for i in CS.cvalues]
    ax.legend(CS.collections + handles, cvalues + labels,
            loc = 'center left', bbox_to_anchor = (1.025, 0.5),
            fontsize = 12, numpoints = 1)
    #text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    return(fig, ax, CS)

def ion():
    plt.ion()

def show():
    plt.show()
