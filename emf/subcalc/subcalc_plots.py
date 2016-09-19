from .. import np, mpl, plt

from ..emf_plots import _save_fig

import subcalc_funks

#rcparams for more static global formatting changes
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.figsize'] = (12, 6)
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.color'] = (.2, .2, .2)
mpl.rcParams['axes.labelcolor'] = (.2, .2, .2)
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['legend.borderaxespad'] = 0 #mpl default is None
mpl.rcParams['xtick.color'] = (.2, .2, .2)
mpl.rcParams['ytick.color'] = (.2, .2, .2)

#other more specific/dynamic global formatting variables
_contour_cmap_name = 'viridis_r' #http://matplotlib.org/examples/color/colormaps_reference.html
_contour_linewidths = 2
_contour_alpha = 0.8
_footprint_alpha = .75
_footprint_zorder = -1
_footprint_label_fontsize = 15
_POPC_label_fontsize = 13 #points of potential concern on footprints
_field_marker_size = 8 #marker size for points of potential concern
_field_marker_edgewidth = 1.5
_ax_frameon = True
_ax_ticks_on = False
_leg_edge_on = False

#deal with colormapping
_contour_cmap = mpl.cm.get_cmap(_contour_cmap_name)
if(not hasattr(_contour_cmap, 'colors')):
    _contour_cmap = mpl.colors.ListedColormap(
            mpl.colors.makeMappingArray(256, _contour_cmap),
            name=_contour_cmap_name,
            N=256)

def ion():
    """Call plt.ion() to toggle interactive plotting on"""
    plt.ion()

def show():
    """Call plt.show() to display open figures"""
    plt.show()

def close(*args):
    """Call plt.close() on any Figure objects or lists of Figure objects
    passed in. If nothing is passed, all Figure objects are closed with
    plt.close('all')"""
    if(args):
        for a in args:
            if(hasattr(a, '__len__')):
                for b in a:
                    plt.close(b.number)
            else:
                plt.close(a.number)
    else:
        plt.close('all')

def _format_ax(ax):
    #apply legend formatting
    leg = ax.get_legend()
    if(leg):
        rec = leg.get_frame()
        if(not _leg_edge_on):
            rec.set_edgecolor('white')
    #apply axis formatting
    ax.set_frame_on(_ax_frameon)
    if(not _ax_ticks_on):
        ax.tick_params(axis = 'both', which = 'both',
            bottom = 'off', top = 'off', left = 'off', right = 'off')

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
        ax.text(x, y, fps[0].group, ha = 'center', va = 'center',
                fontsize = _footprint_label_fontsize + 2)
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
        ax.text(x, y, fp.group, fontsize = _footprint_label_fontsize)

def _make_log_color_indexer(Bmin, Bmax, L_cmap):
    """Decorator function to make a lambda function that will map field
    magnitudes to color indices in a colormap object
    args:
        Bmin - the minimum field value, is mapped to 0
        Bmax - the maximum field value, is mapped to the last color in the cmap
        L_cmap - the length of the cmap, or the number of colors in it
    returns:
        f - a lambda function mapping field to cmap index"""
    m = L_cmap/(np.log10(Bmax) - np.log10(Bmin))
    b = - m*np.log10(Bmin)
    def mlci(x):
        idx = int(round(m*np.log10(x) + b))
        if(idx < 0):
            idx = 0
        elif(idx >= L_cmap):
            idx = L_cmap - 1
        return(idx)
    return(mlci)

def plot_contours(mod, **kwargs):
    """Generate a contour plot from the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kwargs:
        levels - iterable, determines contour level values
                    (automatically determined by matplotlib if not used)
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
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
    fig = plt.figure(figsize = (14 + 4, 14*aspect_ratio))
    ax = fig.add_subplot(1,1,1)
    #adjust axis dimensions by transforming between figure and display coords
    box = ax.get_position()
    T = fig.transFigure.transform
    I = fig.transFigure.inverted().transform
    width = I((T((0., box.height))[1]/aspect_ratio, 0.))[0]
    ax.set_position([box.x0, box.y0, width, box.height])
    #plot contours
    if('levels' in kwargs):
        CS = ax.contour(mod.X, mod.Y, mod.B,
                cmap=_contour_cmap,
                levels=np.array(kwargs['levels'], dtype = float),
                locator=mpl.ticker.LogLocator(),
                zorder=-1, alpha=_contour_alpha)
    else:
        CS = ax.contour(mod.X, mod.Y, mod.B,
                cmap=_contour_cmap,
                locator=mpl.ticker.LogLocator(),
                zorder=-1, alpha=_contour_alpha)
    #set the contour linewidths private variable...
    for c in ax.collections:
        c._linewidths = (_contour_linewidths,)
    #plot location of maximum field
    handles, labels = [], []
    peak_B, yidx, xidx = subcalc_funks._2Dmax(mod.B)
    peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
    handles.append(ax.plot(mod.x[xidx], mod.y[yidx], 'o',
            markersize = _field_marker_size,
            markerfacecolor=_contour_cmap.colors[-1],
            markeredgewidth=_field_marker_edgewidth,
            markeredgecolor='r')[0])
    labels.append('Maximum Modeled\nMagnetic Field\n(%s mG)' % peak_B)
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
                    alpha = _footprint_alpha,
                    zorder = _footprint_zorder)[0])
            labels.append(fp.group)
            #remove the power line footprint
            fps.pop(idx)
        #plot the other footprints
        of_concern = False
        mlci = _make_log_color_indexer(min(CS.cvalues),
                max(CS.cvalues), len(_contour_cmap.colors))
        for fp in fps:
            ax.plot(fp.x, fp.y, 'k',
                    alpha = _footprint_alpha,
                    zorder = _footprint_zorder)
            #mark maximum field if the footprint is 'of concern'
            if(fp.of_concern):
                x = fp.x
                y = fp.y
                #interpolate
                B_interp = mod.interp(x, y)
                #store index of max value on footprint
                idx = np.argmax(B_interp)
                #get associated color for marker
                cidx = mlci(B_interp[idx])
                #stringify the max value
                m = str(subcalc_funks._sig_figs(B_interp[idx], 3))
                #plot
                h = ax.plot(x[idx], y[idx], 'o',
                        markeredgecolor='k',
                        markerfacecolor=_contour_cmap.colors[cidx],
                        markeredgewidth = _field_marker_edgewidth,
                        markersize = _field_marker_size)[0]
                ax.text(x[idx] + xmarg, y[idx] + ymarg, m,
                        fontsize = _POPC_label_fontsize)
                if(of_concern is False):
                    handles.append(h)
                    labels.append('Nearest Points of Neighboring\nLots and Buildings')
                    of_concern = True

        #deal with labeling
        _label_footprint_group(ax, fps)

    #legend
    cvalues = [str(subcalc_funks._check_intable(i)) + ' mG' for i in CS.cvalues]
    ax.legend(CS.collections + handles, cvalues + labels,
            loc = 'center left', bbox_to_anchor = (1.025, 0.5),
            numpoints = 1)

    #put North arrow in top right 10% of the axes, if mod.north_angle is set
    if(mod.north_angle is not None):
        #size of box for North arrow, in axes fraction
        scale = 0.15
        #get figure pixel coordinates of axes
        X, Y = mod.grid_limits['xmax'], mod.grid_limits['ymax']
        #use "smaller" dimension for arrow magnitude
        if(X <= Y):
            mag = X*0.05
        else:
            mag = Y*0.05
        #center of arrow
        x, y = X - 2*mag, Y - 2*mag
        #calculate offsets from angle
        x_offset = mag*np.sin(mod.north_angle*(2*np.pi/360.))
        y_offset = mag*np.cos(mod.north_angle*(2*np.pi/360.))
        #draw arrow
        ax.annotate('', xy = (x + x_offset, y + y_offset),
                xycoords = 'data',
                xytext = (x - x_offset, y - y_offset),
                textcoords = 'data',
                zorder = 2,
                arrowprops = {'arrowstyle': 'simple',
                                'facecolor': 'gray',
                                'edgecolor': 'gray',
                                'shrinkB': 10,
                                'alpha': 0.8})
        #write N for north
        ax.annotate('$N$', xy = (x + x_offset, y + y_offset),
                xycoords = 'data', zorder = 2, color = 'gray',
                ha = 'center', va = 'center', fontsize = 16)

    #text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')
    #final formatting
    _format_ax(ax)

    #saving
    _save_fig('contour_plot', fig, **kwargs)

    return(fig, ax, CS)
