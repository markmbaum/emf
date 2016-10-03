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
_default_contour_cmap_name = 'viridis_r' #http://matplotlib.org/examples/color/colormaps_reference.html
_default_pcolormesh_cmap_name = 'magma_r'
_contour_linewidths = 2
_contour_alpha = 0.8
_pcolor_alpha = 0.35
_footprint_alpha = .75
_footprint_zorder = -1
_footprint_label_fontsize = 15
_POPC_label_fontsize = 13 #points of potential concern on footprints
_field_marker_size = 8 #marker size for points of potential concern
_field_marker_edgewidth = 1.5
_ax_frameon = True
_ax_ticks_on = False
_leg_edge_on = False
_be_concerned = True

#get the colormap
def get_cmap(cmap_name):
    """Set the global colormap by passing the name of a matplotlib cmap, see:
    http://matplotlib.org/examples/color/colormaps_reference.html"""

    cmap = mpl.cm.get_cmap(cmap_name)
    if(not hasattr(cmap, 'colors')):
        cmap = mpl.colors.ListedColormap(
                mpl.colors.makeMappingArray(256, cmap),
                name=cmap_name,
                N=256)
    return(cmap)

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
        ax.tick_params(axis='both', which='both',
            bottom='off', top='off', left='off', right='off')

def _set_ax_aspect(fig, ax, aspect):
    """Make the axes comform to a certain aspect ratio (height/width)
    args:
        fig - Figure object containing the Axes
        ax - Axes object to adjust
        aspect - desired aspect ratio (height/width in this case)"""
    #adjust axis dimensions by transforming between figure and display coords
    box = ax.get_position()
    T = fig.transFigure.transform
    I = fig.transFigure.inverted().transform
    width = I((T((0., box.height))[1]/aspect, 0.))[0]
    ax.set_position([box.x0, box.y0, width, box.height])

def _plot_footprints(ax, mod, cmap, mlci):
    """Plot and label footprint outlines
    args:
        ax - Axes object to plot in
        mod - Model object containing the Footprint objects
        cmap - Colormap object for coloring points of concern
        mlci - function for indexing colormap
               (from _make_linear_color_indexer or _make_log_color_indexer)
    returns:
        handles - list of plotted handles for legend
        labels - list of labels for the handles"""

    #get some grid info
    xmax, xmin = np.max(mod.x), np.min(mod.x)
    ymax, ymin = np.max(mod.y), np.min(mod.y)
    xmarg, ymarg = (xmax - xmin)*0.005, (ymax - ymin)*0.005

    handles, labels = [], []
    #plot footprints
    of_concern = False #flag for appending handles to legend list
    power_line_handle_idx = False
    for g in mod.footprint_groups:
        #get footprints in the group and the group name
        fps = [mod.footprints[i] for i in g]
        group_name = fps[0].group
        #plot powerline footprings
        i = 0
        L = len(fps)
        while(i < L):
            fp = fps[i]
            if(fp.power_line):
                h = ax.plot(fp.x, fp.y, 'k--', linewidth=2,
                        alpha=_footprint_alpha,
                        zorder=_footprint_zorder)[0]
                if(power_line_handle_idx is False):
                    handles.append(h)
                    labels.append(fp.group)
                    power_line_handle_idx = len(handles) - 1
                else:
                    labels[power_line_handle_idx] += ';\n' + fp.group
                #remove the power line footprint
                fps.pop(i)
                L -= 1
                i -= 1
            i += 1
        for fp in fps:
            ax.plot(fp.x, fp.y, 'k',
                    alpha=_footprint_alpha,
                    zorder=_footprint_zorder)
            #mark maximum field if the footprint is 'of concern'
            if(fp.of_concern and _be_concerned):
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
                        markerfacecolor=cmap.colors[cidx],
                        markeredgewidth=_field_marker_edgewidth,
                        markersize=_field_marker_size)[0]
                ax.text(x[idx] + xmarg, y[idx] + ymarg, m,
                        fontsize=_POPC_label_fontsize)
                if(of_concern is False):
                    handles.append(h)
                    labels.append('Nearest Points of Neighboring\nLots and Buildings')
                    of_concern = True
        #deal with labeling
        _label_footprint_group(ax, fps)
    #return the legend items
    return(handles, labels)

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
        ax.text(x, y, fps[0].group, ha='center', va='center',
                fontsize=_footprint_label_fontsize + 2)
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
        ax.text(x, y, fp.group, fontsize=_footprint_label_fontsize,
                ha='center', va='center')

def _draw_north_arrow(ax, angle):
    """Draw a north arrow in the upper right corner of the axes
    args:
        ax - Axes object to draw in
        angle - number, direction of north arrow
                (0 is along +y axis and clockwise is increasing)"""
    #size of box for North arrow, in axes fraction
    scale = 0.15
    #get figure pixel coordinates of axes
    X, Y = max(ax.get_xlim()), max(ax.get_ylim())
    #use "smaller" dimension for arrow magnitude
    if(X <= Y):
        mag = X*0.05
    else:
        mag = Y*0.05
    #center of arrow
    x, y = X - 2*mag, Y - 2*mag
    #calculate offsets from angle
    x_offset = mag*np.sin(angle*(2*np.pi/360.))
    y_offset = mag*np.cos(angle*(2*np.pi/360.))
    #draw arrow
    ax.annotate('', xy=(x + x_offset, y + y_offset),
            xycoords='data',
            xytext=(x - x_offset, y - y_offset),
            textcoords='data',
            zorder=2,
            arrowprops={'arrowstyle': 'simple',
                            'facecolor': 'none',
                            'edgecolor': 'black',
                            'shrinkB': 10,
                            'alpha': 0.8})
    #write N for north
    ax.annotate('$N$', xy=(x + x_offset, y + y_offset),
            xycoords='data', zorder=2, color='black',
            ha='center', va='center', fontsize=16)

def _write_Bkey(ax, mod):
    """Place a small text object just outside the upper right corner of the
    axes indicating the component of the magnetic field represented by
    the results
    args:
        ax - Axes to plot in
        mod - Model object"""
    xl, yl = ax.get_xlim(), ax.get_ylim()
    #xmarg, ymarg = 0.005*(xl[1] - xl[0]), 0.005*(yl[1] - yl[0])
    ax.text(xl[1], yl[1],
            'Field Component: %s' % mod.Bkey,
            fontsize=9, ha='right', va='bottom')

def _make_color_indexer(Bmin, Bmax, L_cmap, scale):
    """Decorator function to make a lambda function that will map field
    magnitudes to color indices in a colormap object, on a linear scale
    args:
        Bmin - the minimum field value, is mapped to 0
        Bmax - the maximum field value, is mapped to the last color in the cmap
        L_cmap - the length of the cmap, or the number of colors in it
        scale - string, 'log' or 'linear'
    returns:
        f - a function mapping field values to cmap index"""
    if(scale == 'linear'):
        m = float(L_cmap)/(Bmax - Bmin)
        b = - m*Bmin

        def mlci(x):
            idx = int(round(m*x + b))
            if(idx < 0):
                idx = 0
            elif(idx >= L_cmap):
                idx = L_cmap - 1
            return(idx)

    elif(scale == 'log'):
        m = L_cmap/(np.log10(Bmax) - np.log10(Bmin))
        b = - m*np.log10(Bmin)

        def mlci(x):
            idx = int(round(m*np.log10(x) + b))
            if(idx < 0):
                idx = 0
            elif(idx >= L_cmap):
                idx = L_cmap - 1
            return(idx)

    else:
        raise(EMFError("""
        Illegal scale string "%s" passed to _make_color_indexer.
        Must be 'log' or 'linear'."""))

    return(mlci)

def plot_contour(mod, **kwargs):
    """Generate a contour plot from the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kwargs:
        scale - str, can be 'log' or 'linear' (default is 'linear')
        levels - iterable, determines contour level values
                    (automatically determined by matplotlib if not used)
        north_angle - number, direction of north arrow
                        (0 is along +y axis and clockwise is increasing)
        cmap - str, name of matplotlib colormap, see:
                http://matplotlib.org/examples/color/colormaps_reference.html
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        CS - matplotlib QuadContourSet object"""

    #check kwargs
    if('scale' in kwargs):
        scale = kwargs['scale']
    else:
        scale = 'linear'
    if('cmap' in kwargs):
        _cmap = get_cmap(kwargs['cmap'])
    else:
        _cmap = get_cmap(_default_contour_cmap_name)

    #assess dimensions of grid
    xmax, xmin = np.max(mod.x), np.min(mod.x)
    ymax, ymin = np.max(mod.y), np.min(mod.y)
    aspect_ratio = (ymax - ymin)/(xmax - xmin)

    #generate the figure
    fig = plt.figure(figsize=(14 + 4, 14*aspect_ratio))
    ax = fig.add_subplot(1,1,1)

    #adjust axis dimensions by transforming between figure and display coords
    _set_ax_aspect(fig, ax, aspect_ratio)

    #get contour plotting keyword arguments
    contour_kwargs = {'cmap': _cmap, 'zorder': -1, 'alpha': _contour_alpha}
    #get levels if passed in
    if('levels' in kwargs):
        contour_kwargs['levels'] = np.array(kwargs['levels'], dtype=float)
    #determine contour scaling kwarg
    if(scale == 'log'):
        contour_kwargs['locator'] = mpl.ticker.LogLocator()
    elif(scale == 'linear'):
        contour_kwargs['locator'] = mpl.ticker.MaxNLocator()
    else:
        raise(EMFError("""
        Invalid value passed to keyword argument 'scale.'
        Must be 'log' or 'linear.'"""))

    #plot
    CS = ax.contour(mod.X, mod.Y, mod.B, **contour_kwargs)

    #get color indexer from decorator along the way
    mlci = _make_color_indexer(np.min(CS.cvalues), np.max(CS.cvalues),
                            len(_cmap.colors), scale)

    #set the contour linewidths private variable...
    for c in ax.collections:
        c._linewidths = (_contour_linewidths,)

    #plot location of maximum field
    handles, labels = [], []
    peak_B, yidx, xidx = subcalc_funks._2Dmax(mod.B)
    peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
    handles.append(ax.plot(mod.x[xidx], mod.y[yidx], 'o',
            markersize=_field_marker_size,
            markerfacecolor=_cmap.colors[-1],
            markeredgewidth=_field_marker_edgewidth,
            markeredgecolor='r')[0])
    labels.append('Maximum Modeled\nMagnetic Field\n(%s mG)' % peak_B)

    #plot footprints
    H, L = _plot_footprints(ax, mod, _cmap, mlci)
    handles += H
    labels += L

    #legend
    cvalues = [str(subcalc_funks._check_intable(i)) + ' mG' for i in CS.cvalues]
    ax.legend(CS.collections + handles, cvalues + labels,
            loc='center left', bbox_to_anchor=(1.025, 0.5),
            numpoints=1)

    #north arrow
    if(mod.north_angle is not None):
        _draw_north_arrow(ax, mod.north_angle)
    elif('north_angle' in kwargs):
        _draw_north_arrow(ax, kwargs['north_angle'])

    #write Bkey note
    _write_Bkey(ax, mod)

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #saving
    _save_fig('contour_plot', fig, **kwargs)

    return(fig, ax, CS)

def plot_pcolormesh(mod, **kwargs):
    """Generate a color mesh plot of the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kwargs:
        north_angle - number, direction of north arrow
                        (0 is along +y axis and clockwise is increasing)
        cmap - str, name of matplotlib colormap, see:
                http://matplotlib.org/examples/color/colormaps_reference.html
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        QM - matplotlib QuadMesh object
        cbar - matplotlib Colorbar object (to which legend is attached)"""

    #check kwargs
    if('cmap' in kwargs):
        _cmap = get_cmap(kwargs['cmap'])
    else:
        _cmap = get_cmap(_default_pcolormesh_cmap_name)

    #assess dimensions of grid
    xmax, xmin = np.max(mod.x), np.min(mod.x)
    ymax, ymin = np.max(mod.y), np.min(mod.y)
    aspect_ratio = (ymax - ymin)/(xmax - xmin)

    #generate the figure
    fig = plt.figure(figsize=(14 + 6, 14*aspect_ratio))
    ax = fig.add_subplot(1,1,1)

    #adjust axis dimensions by transforming between figure and display coords
    _set_ax_aspect(fig, ax, aspect_ratio)
    ax.set_xlim(mod.grid_limits['xmin'], mod.grid_limits['xmax'])
    ax.set_ylim(mod.grid_limits['ymin'], mod.grid_limits['ymax'])

    #get plotting keyword arguments
    pcolor_kwargs = {'cmap': _cmap, 'alpha': _pcolor_alpha, 'zorder': -1,
                    'shading': 'gouraud'}

    #plot
    QM = ax.pcolormesh(mod.X, mod.Y, mod.B, **pcolor_kwargs)

    #get color indexer from decorator along the way
    mlci = _make_color_indexer(np.min(mod.B), np.max(mod.B),
                            len(_cmap.colors), 'linear')

    #plot location of maximum field
    handles, labels = [], []
    peak_B, yidx, xidx = subcalc_funks._2Dmax(mod.B)
    peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
    handles.append(ax.plot(mod.x[xidx], mod.y[yidx], 'o',
            markersize=_field_marker_size,
            markerfacecolor=_cmap.colors[-1],
            markeredgewidth=_field_marker_edgewidth,
            markeredgecolor='r')[0])
    labels.append('Maximum Modeled\nMagnetic Field\n(%s mG)' % peak_B)

    #plot footprints
    H, L = _plot_footprints(ax, mod, _cmap, mlci)
    handles += H
    labels += L

    #colorbar and legend
    cbar = fig.colorbar(QM)
    cbar.ax.legend(handles, labels,
            loc='center left', bbox_to_anchor=(3, 0.5),
            numpoints=1)
    _format_ax(cbar.ax)

    #north arrow
    if(mod.north_angle is not None):
        _draw_north_arrow(ax, mod.north_angle)
    elif('north_angle' in kwargs):
        _draw_north_arrow(ax, kwargs['north_angle'])

    #write Bkey note
    _write_Bkey(ax, mod)

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #saving
    _save_fig('contour_plot', fig, **kwargs)

    return(fig, ax, QM, cbar)
