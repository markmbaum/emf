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
_max_fig_width = 18
_max_fig_height = 10
_pcolormesh_legend_padding = 6
_contour_legend_padding = 4
_default_pcolormesh_cmap_name = 'magma_r' #http://matplotlib.org/examples/color/colormaps_reference.html
_default_contour_cmap_name = 'viridis_r'
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

def _equal_ax_objs(x_range, y_range, max_width, max_height, legend_pad):
    """Get a matplotlib Figure with an Axes object having the same
    aspect ratio as the data to be plotted in it.
    args:
        x_range - float, the range of values that will be plotted on the
                  horizontal axis
        y_range - float, the range of values that will be plotted on the
                  vertical axis
        max_width - float, the maximum allowable width of the Figure
                    in inches
        max_height - float, the maximum allowable height of the Figure
                     in inches
        legend_pad - float, inches, padding to the right of the Axes
                     for a legend
    returns:
        fig - Figure object
        ax - Axes object"""
    #calculate aspect ratio
    aspect = abs(float(y_range)/x_range)
    #find which dimension of the figure should be set to max_inches
    if((max_height/aspect + legend_pad) > max_width):
        fig_x = max_width
        fig_y = (max_width - legend_pad)*aspect
    else:
        fig_x = max_height/aspect + legend_pad
        fig_y = max_height
    #generate plot objects
    fig, ax = plt.subplots(1, 1, figsize=(fig_x,fig_y))
    #adjust the aspect ratio of the axes
    _set_ax_aspect(fig, ax, aspect)
    #return
    return(fig, ax)

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

    #get margins for labeling
    xmarg, ymarg = (mod.xmax - mod.xmin)*0.005, (mod.ymax - mod.ymin)*0.005

    handles, labels = [], []
    #plot footprints
    of_concern = False #flag for appending handles to legend list
    power_line_handle_idx = False
    for fps in mod.footprint_groups:

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
                    labels.append('Nearest Points of\nNeighboring Structures')
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

def plot_contour(mod, **kw):
    """Generate a contour plot from the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kw:
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
        max_fig_width - float/int, inches, maximum width of figure
        max_fig_height - float/int, inches, maximum height of figure
        legend_padding - float/int, inches, width left for legend area
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        CS - matplotlib QuadContourSet object"""

    #check kw
    if('scale' in kw):
        scale = kw['scale']
    else:
        scale = 'linear'
    if('cmap' in kw):
        _cmap = get_cmap(kw['cmap'])
    else:
        _cmap = get_cmap(_default_contour_cmap_name)
    if('max_fig_width' in kw):
        max_fig_width = kw['max_fig_width']
    else:
        max_fig_width = _max_fig_width
    if('max_fig_height' in kw):
        max_fig_height = kw['max_fig_height']
    else:
        max_fig_height = _max_fig_height
    if('legend_padding' in kw):
        legend_padding = kw['legend_padding']
    else:
        legend_padding = _contour_legend_padding

    #generate the plot objects
    fig, ax = _equal_ax_objs(mod.xmax - mod.xmin, mod.ymax - mod.ymin,
            max_fig_width, max_fig_height, legend_padding)
    ax.set_xlim(mod.xmin, mod.xmax)
    ax.set_ylim(mod.ymin, mod.ymax)

    #get contour plotting keyword arguments
    contour_kw = {'cmap': _cmap, 'zorder': -1, 'alpha': _contour_alpha}
    #get levels if passed in
    if('levels' in kw):
        contour_kw['levels'] = np.array(kw['levels'], dtype=float)
    #determine contour scaling kwarg
    if(scale == 'log'):
        contour_kw['locator'] = mpl.ticker.LogLocator()
    elif(scale == 'linear'):
        contour_kw['locator'] = mpl.ticker.MaxNLocator()
    else:
        raise(EMFError("""
        Invalid value passed to keyword argument 'scale.'
        Must be 'log' or 'linear.'"""))

    #plot
    CS = ax.contour(mod.X, mod.Y, mod.B, **contour_kw)

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
    elif('north_angle' in kw):
        _draw_north_arrow(ax, kw['north_angle'])

    #write Bkey note
    _write_Bkey(ax, mod)

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #saving
    _save_fig('contour_plot', fig, **kw)

    return(fig, ax, CS)

def plot_pcolormesh(mod, **kw):
    """Generate a color mesh plot of the magnetic field results and
    Footprint objects stored in a Model object
    args:
        mod - Model object
    kw:
        north_angle - number, direction of north arrow
                        (0 is along +y axis and clockwise is increasing)
        cmap - str, name of matplotlib colormap, see:
                http://matplotlib.org/examples/color/colormaps_reference.html
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
        max_fig_width - float/int, inches, maximum width of figure
        max_fig_height - float/int, inches, maximum height of figure
        legend_padding - float/int, inches, width left for legend area
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        QM - matplotlib QuadMesh object
        cbar - matplotlib Axes that the colorbar is drawn in (to which
                legend is attached)"""

    #check kw
    if('cmap' in kw):
        _cmap = get_cmap(kw['cmap'])
    else:
        _cmap = get_cmap(_default_pcolormesh_cmap_name)
    if('max_fig_width' in kw):
        max_fig_width = kw['max_fig_width']
    else:
        max_fig_width = _max_fig_width
    if('max_fig_height' in kw):
        max_fig_height = kw['max_fig_height']
    else:
        max_fig_height = _max_fig_height
    if('legend_padding' in kw):
        legend_padding = kw['legend_padding']
    else:
        legend_padding = _pcolormesh_legend_padding

    #generate the plot objects
    fig, ax = _equal_ax_objs(mod.xmax - mod.xmin, mod.ymax - mod.ymin,
            max_fig_width, max_fig_height, legend_padding)
    ax.set_xlim(mod.xmin, mod.xmax)
    ax.set_ylim(mod.ymin, mod.ymax)

    #get plotting keyword arguments
    pcolor_kw = {'cmap': _cmap, 'alpha': _pcolor_alpha, 'zorder': -1,
                    'shading': 'gouraud'}

    #plot
    QM = ax.pcolormesh(mod.X, mod.Y, mod.B, **pcolor_kw)

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
    cbar = fig.add_subplot(1,2,2)
    box = ax.get_position()
    fig_size = fig.get_size_inches()
    w = 0.35/fig_size[0]
    cbar.set_position([box.x0 + box.width + w, box.y0, w, box.height])
    cbar = fig.colorbar(QM, cax=cbar, drawedges=False)
    cbar.ax.legend(handles, labels,
            loc='center left', bbox_to_anchor=(2.5, 0.5),
            numpoints=1)
    _format_ax(cbar.ax)

    #north arrow
    if(mod.north_angle is not None):
        _draw_north_arrow(ax, mod.north_angle)
    elif('north_angle' in kw):
        _draw_north_arrow(ax, kw['north_angle'])

    #write Bkey note
    _write_Bkey(ax, mod)

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #saving
    _save_fig('contour_plot', fig, **kw)

    return(fig, ax, QM, cbar)
