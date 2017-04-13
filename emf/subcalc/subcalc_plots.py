from .. import os, np, mpl, plt, textwrap, _Rectangle, _Axes3D

from ..emf_plots import _save_fig, _prepare_fig

import subcalc_class
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
_contour_linewidths = 2
_contour_alpha = 0.8
_pcolor_alpha = 0.35
_footprint_alpha = .75
_footprint_zorder = -1
_footprint_label_fontsize = 15
_POPC_label_fontsize = 13 #points of potential concern on footprints
_field_scatter_size = 25 #marker size for points of potential concern
_field_marker_size = 8 #marker size for max field point
_field_marker_edgewidth = 1.5
_ax_frameon = True
_ax_ticks_on = False
_leg_edge_on = False
be_concerned = True

#get the colormap
def _get_cmap(cmap_name):
    """Set the global colormap by passing the name of a matplotlib cmap, see:
    http://matplotlib.org/examples/color/colormaps_reference.html"""

    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError as e:
        print(str(e) + '\n\nFalling back to "YlOrRd"')
        cmap = plt.get_cmap('YlOrRd')
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
    """Call plt.close() on any Figure objects or lists of Figure objects passed in. If nothing is passed, all Figure objects are closed."""
    if(args):
        for a in args:
            if(hasattr(a, '__len__')):
                for b in a:
                    plt.close(b.number)
            else:
                plt.close(a.number)
    else:
        plt.close('all')

def _get_text_alignment(point1, point2):
    """Get the horizontal and vertical text alignment keywords for text placed at the end of a line segment from point1 to point2
    args:
        point1 - x,y pair
        point2 - x,y pair
    returns:
        ha - horizontal alignment string
        va - vertical alignment string"""
    x1, x2, y1, y2 = point1[0], point2[0], point1[1], point2[1]
    if(x1 < x2):
        ha = 'left'
    else:
        ha = 'right'
    if(y1 < y2):
        va = 'bottom'
    else:
        va = 'top'
    return(ha, va)

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
    #set axis equal
    ax.axis('equal')
    #return
    return(fig, ax)

def _plot_footprints(ax, res, cmap, ci, fp_max, fp_text):
    """Plot and label footprint outlines
    args:
        ax - Axes object to plot in
        res - Results object containing the Footprint objects
        cmap - Colormap object for coloring points of concern
        ci - function for indexing colormap (from _make_color_indexer)
        fp_max - bool, toggle whether the maximum fields along footprints
                 "of concern" are labeled
        fp_text - bool, toggle footprint group labeling with text
    returns:
        handles - list of plotted handles for legend
        labels - list of labels for the handles"""

    #get margins for labeling
    xmarg, ymarg = (res.xmax - res.xmin)*0.005, (res.ymax - res.ymin)*0.005

    handles, labels = [], []
    #plot footprints
    of_concern = False #flag for appending handles to legend list
    power_line_handle_idx = False
    for fps in res.footprint_groups:

        group_name = fps[0].group
        #plot powerline footprings
        i = 0
        L = len(fps)
        #count the number of power line footprints
        N_power_lines = sum(1 for i in res.footprints if (i.power_line))
        #only label each power line if there are fewer than 8
        if(N_power_lines < 8):
            label_all_power_lines = True
        else:
            label_all_power_lines = False
        count = 0
        while(i < L):
            fp = fps[i]
            if(fp.power_line):
                count += 1
                h = ax.plot(fp.x, fp.y, 'k--', linewidth=2,
                        alpha=_footprint_alpha,
                        zorder=_footprint_zorder)[0]
                if(power_line_handle_idx is False):
                    handles.append(h)
                    if((not label_all_power_lines) and (count == 1)):
                        labels.append('Power Line Paths')
                    else:
                        labels.append(textwrap.fill(fp.group, width=25))
                    power_line_handle_idx = len(handles) - 1
                elif(label_all_power_lines):
                    labels[power_line_handle_idx] += ';\n' + fp.group
                #remove the power line footprint
                fps.pop(i)
                L -= 1
                i -= 1
            i += 1
        #plot footprint outlines
        for fp in fps:
            ax.plot(fp.x, fp.y, 'k',
                    alpha=_footprint_alpha,
                    zorder=_footprint_zorder)

        #label the footprints themselves
        if(fp_text):
            _label_footprint_group(ax, fps)

    #mark maximum fields on footprints 'of concern'
    if(be_concerned and fp_max):
        #get fields
        x, y, B = res.concern_points()
        #plot points and label
        for i in range(len(x)):
            h, = ax.plot(x[i], y[i], 'o',
                    markersize=_field_marker_size,
                    markerfacecolor=cmap.colors[ci(B[i])],
                    markeredgewidth=_field_marker_edgewidth,
                    markeredgecolor='k')
            ax.text(x[i], y[i], str(subcalc_funks._sig_figs(B[i], 2)),
                    fontsize=_POPC_label_fontsize)
            if(i == 0):
                handles.append(h)
                labels.append('Maximum Fields Along\nStructures of Concern')

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

def _write_Bkey(*args):
    """Place a small text object just outside the upper right corner of the
    axes indicating the component of the magnetic field represented by
    the results
    args:
        fig - figure containing Axes to annotate
        ax - Axes to plot in
        Results object or a string for the key"""
    fig, ax = args[0], args[1]
    xl, yl = ax.get_xlim(), ax.get_ylim()
    #xmarg, ymarg = 0.005*(xl[1] - xl[0]), 0.005*(yl[1] - yl[0])
    if(type(args[2]) is subcalc_class.Results):
        s = args[2].Bkey
    else:
        s = str(args[2])
    box = ax.get_position()
    fig.text(box.x0 + box.width, box.y0 + box.height, 'Component: %s' % s,
            fontsize=8, ha='right', va='bottom')

def _make_color_indexer(Bmin, Bmax, L_cmap, scale):
    """Decorator function to make a lambda function that will map field
    magnitudes to color indices in a colormap object, on a linear scale
    args:
        Bmin - the minimum field value, is mapped to 0
        Bmax - the maximum field value, is mapped to the last color in the cmap
        L_cmap - the length of the cmap, or the number of colors in it
        scale - string, 'log' or 'lin'
    returns:
        f - a function mapping field values to cmap index"""
    if(scale == 'lin'):
        m = float(L_cmap)/(Bmax - Bmin)
        b = - m*Bmin

        def ci(x):
            idx = int(round(m*x + b))
            if(idx < 0):
                idx = 0
            elif(idx >= L_cmap):
                idx = L_cmap - 1
            return(idx)

    elif(scale == 'log'):
        m = L_cmap/(np.log10(Bmax) - np.log10(Bmin))
        b = - m*np.log10(Bmin)

        def ci(x):
            idx = int(round(m*np.log10(x) + b))
            if(idx < 0):
                idx = 0
            elif(idx >= L_cmap):
                idx = L_cmap - 1
            return(idx)

    else:
        raise(subcalc_class.EMFError("""
        Illegal scale string "%s" passed to _make_color_indexer.
        Must be 'log' or 'lin'."""))

    return(ci)

def _make_cbar(ax, min_val, max_val, cmap, ci):
    "fill an axis in with colored rectangles to make a colorbar"
    #number of colors to draw
    L = len(cmap.colors)
    #get the height of the divisions to be filled by rectangles
    h = (max_val - min_val)/L
    for i in range(L):
        color = cmap.colors[ci(min_val + i*h)]
        rect = _Rectangle((0, i*h + min_val), 1, h*1.01,
                facecolor=color, edgecolor=color, zorder=(i - L))
        ax.add_patch(rect)
    #limits
    ax.set_xlim(0, 1)
    ax.set_ylim(min_val, max_val)
    #formatting
    ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.tick_params(axis='y', direction='out')
    ax.set_yticklabels([('%g mG' % i) for i in ax.get_yticks()])

def plot_contour(res, scale='lin', show_max=True, fp_max=True, fp_text=True,
        cmap='viridis_r', max_fig_width=12, max_fig_height=8, legend_padding=3,
        **kw):
    """Generate a contour plot from the magnetic field results and
    Footprint objects stored in a Results object
    args:
        res - Results object
    optional args:
        scale - str, can be 'log' or 'lin' (default is 'lin')
        show_max - bool, toggle labeling of the maximum field location,
                    default is True
        fp_max - bool, toggle whether the maximum fields along footprints
                 "of concern" are labeled
        fp_text - bool, toggle footprint group labeling with text
        cmap - str, name of matplotlib colormap, see:
               http://matplotlib.org/examples/color/colormaps_reference.html
        max_fig_width - float/int, inches, maximum width of figure
        max_fig_height - float/int, inches, maximum height of figure
        legend_padding - float/int, inches, width left for legend area
    kw:
        levels - iterable, determines contour level values
                 (automatically determined by matplotlib if unused)
        north_angle - number, direction of north arrow
                      (0 is along +y axis and clockwise is increasing)

        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        CS - matplotlib QuadContourSet object"""

    #get a colormap object
    cmap = _get_cmap(cmap)

    #generate the plot objects
    fig, ax = _equal_ax_objs(res.xmax - res.xmin, res.ymax - res.ymin,
            max_fig_width, max_fig_height, legend_padding)
    ax.set_xlim(res.xmin, res.xmax)
    ax.set_ylim(res.ymin, res.ymax)

    #get contour plotting keyword arguments
    contour_kw = {'cmap': cmap, 'zorder': -1, 'alpha': _contour_alpha}
    #get levels if passed in
    if('levels' in kw):
        contour_kw['levels'] = np.array(kw['levels'], dtype=float)
    #determine contour scaling kwarg
    if(scale == 'log'):
        contour_kw['locator'] = mpl.ticker.LogLocator()
    elif(scale == 'lin'):
        contour_kw['locator'] = mpl.ticker.MaxNLocator()
    else:
        raise(subcalc_class.EMFError("""
        Invalid value passed to keyword argument 'scale.'
        Must be 'log' or 'lin.'"""))

    #plot
    CS = ax.contour(res.X, res.Y, res.B, **contour_kw)

    #get color indexer from decorator along the way
    ci = _make_color_indexer(np.min(CS.cvalues), np.max(CS.cvalues),
                            len(cmap.colors), scale)

    #set the contour linewidths private variable...
    for c in ax.collections:
        c._linewidths = (_contour_linewidths,)

    #keep list of handles and labels for the legend
    handles, labels = [], []

    #plot location of maximum field
    if(show_max):
        peak_B, yidx, xidx = subcalc_funks._2Dmax(res.B)
        peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
        handles.append(ax.plot(res.x[xidx], res.y[yidx], 'o',
                markersize=_field_marker_size,
                markerfacecolor=cmap.colors[-1],
                markeredgewidth=_field_marker_edgewidth,
                markeredgecolor='r')[0])
        labels.append('Maximum Modeled\nMagnetic Field\n(%s mG)' % peak_B)

    #plot footprints
    H, L = _plot_footprints(ax, res, cmap, ci, fp_max, fp_text)
    handles += H
    labels += L

    #legend
    cvalues = [str(subcalc_funks._check_intable(i)) + ' mG' for i in CS.cvalues]
    ax.legend(CS.collections + handles, cvalues + labels,
            loc='center left', bbox_to_anchor=(1.025, 0.5),
            numpoints=1)

    #north arrow
    if(res.north_angle is not None):
        _draw_north_arrow(ax, res.north_angle)
    elif('north_angle' in kw):
        _draw_north_arrow(ax, kw['north_angle'])

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #write Bkey note
    _write_Bkey(fig, ax, res)

    #saving
    _save_fig('contour-plot', fig, **kw)

    return(fig, ax, CS)

def plot_pcolormesh(res, show_max=True, fp_max=True, fp_text=True,
        cmap='magma_r', max_fig_width=12, max_fig_height=8, legend_padding=5,
        **kw):
    """Generate a color mesh plot of the magnetic field results and
    Footprint objects stored in a Results object
    args:
        res - Results object
    optional args:
        show_max - bool, toggle labeling of the maximum field location,
                    default is True
        fp_max - bool, toggle whether the maximum fields along footprints
                 "of concern" are labeled
        fp_text - bool, toggle footprint group labeling with text
        cmap - str, name of matplotlib colormap, see:
               http://matplotlib.org/examples/color/colormaps_reference.html
        max_fig_width - float/int, inches, maximum width of figure
        max_fig_height - float/int, inches, maximum height of figure
        legend_padding - float/int, inches, width left for legend area
    kw:
        north_angle - number, direction of north arrow
                        (0 is along +y axis and clockwise is increasing)
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                Fthe 'save' keyword
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - matplotlib figure object
        ax - matplotlib axes object
        QM - matplotlib QuadMesh object
        ax_cbar - matplotlib Axes that the colorbar is drawn in (to which
                legend is attached)"""

    #get a colormap object
    cmap = _get_cmap(cmap)

    #generate the plot objects
    fig, ax = _equal_ax_objs(res.xmax - res.xmin, res.ymax - res.ymin,
            max_fig_width, max_fig_height, legend_padding)
    ax.set_xlim(res.xmin, res.xmax)
    ax.set_ylim(res.ymin, res.ymax)

    #get plotting keyword arguments
    pcolor_kw = {'cmap': cmap, 'alpha': _pcolor_alpha, 'zorder': -1,
                    'shading': 'gouraud'}

    #plot
    QM = ax.pcolormesh(res.X, res.Y, res.B, **pcolor_kw)

    #get color indexer from decorator along the way
    ci = _make_color_indexer(np.min(res.B), np.max(res.B),
                            len(cmap.colors), 'lin')

    #keep list of handles and labels for the legend
    handles, labels = [], []

    #plot location of maximum field
    if(show_max):
        peak_B, yidx, xidx = subcalc_funks._2Dmax(res.B)
        peak_B = str(subcalc_funks._sig_figs(peak_B, 3))
        handles.append(ax.plot(res.x[xidx], res.y[yidx], 'o',
                markersize=_field_marker_size,
                markerfacecolor=cmap.colors[-1],
                markeredgewidth=_field_marker_edgewidth,
                markeredgecolor='r')[0])
        labels.append('Maximum Modeled\nMagnetic Field\n(%s mG)' % peak_B)

    #plot footprints
    H, L = _plot_footprints(ax, res, cmap, ci, fp_max, fp_text)
    handles += H
    labels += L

    #colorbar and legend
    ax_cbar = fig.add_subplot(1,2,2)
    box = ax.get_position()
    figsize = fig.get_size_inches()
    w = 0.2/figsize[0]
    ax_cbar.set_position([box.x1 + w, box.y0, w, box.height])
    _make_cbar(ax_cbar, res.Bmin, res.Bmax, cmap, ci)
    ax_cbar.legend(handles, labels, loc='center left',
            bbox_to_anchor=(6, 0.5), numpoints=1)
    _format_ax(ax_cbar)

    #north arrow
    if(res.north_angle is not None):
        _draw_north_arrow(ax, res.north_angle)
    elif('north_angle' in kw):
        _draw_north_arrow(ax, kw['north_angle'])

    #axes text
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')

    #final formatting
    _format_ax(ax)

    #write Bkey note
    _write_Bkey(fig, ax, res)

    #saving
    _save_fig('pcolormesh-plot', fig, **kw)

    return(fig, ax, QM, ax_cbar)

def plot_path(obj, points, n=101, x_labeling='distance', scale='lin',
    cmap='viridis_r', edgecolors='gray', Bkey='Bmax', **kw):
    """Plot a Results object's fields along a line segment in the results domian,
    essentially a cross section of the fields
    args:
        obj - Results object or Model object
        points - an iterable of x,y pairs representing a path through the results
                 domain to plot, for example: [(1,2), (1,3), (2,4)]
    optional args:
        n - integer, number of points sampled (default 101)
        x_labeling - 'distance' or 'location', for x axis ticks labeled
                      according to the distance along the segment or with
                      the (x,y) coordinates of sample points,
                      default is 'distance'
        scale - str, color scaling, can be 'log' or 'lin' (default is 'lin')
        cmap - str, name of matplotlib colormap, see:
               http://matplotlib.org/examples/color/colormaps_reference.html
        edgecolors - color or sequence of colors, to set facecolors to
                     edgecolors, use 'face'
        Bkey - str, the Bkey (component of fields) to plot
    kw:
        ax - target Axes
        fig - Figure object, target figure for plotting, overridden by 'ax'
        save - bool, toggle whether the figure is saved
        path - string, path to directory for saved figure. If set, overrides
                the 'save' keyword
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - Figure object
        ax - Axes object"""

    #get a colormap object
    cmap = _get_cmap(cmap)

    #compute the path
    if(type(obj) is subcalc_class.Results):
        #store the original Bkey of the Results, then switch it to the
        #selected key for plotting, to be switched back after interpolation
        init_Bkey = obj.Bkey
        obj.Bkey = Bkey
        #interpolate
        x, y, B_interp = obj.path(points, n)
        #switch Bkey back to original
        obj.Bkey = init_Bkey
        #get the maxima and minima for coloring points
        bmin, bmax = obj.Bmin, obj.Bmax
        #compute distance along path
        dist = subcalc_funks.cumulative_distance(x, y)
    elif(type(obj) is subcalc_class.Model):
        #unpack points
        x, y = zip(*points)
        #sample
        df = obj.sample(x, y, n=n)
        B_interp, dist = df[Bkey].values, df['dist'].values
        #store minima and maxima for coloring
        bmin, bmax = min(B_interp), max(B_interp)

    #get plotting objects
    fig, ax = _prepare_fig(**kw)

    #plot
    ci = _make_color_indexer(bmin, bmax, len(cmap.colors), scale)
    colors = [cmap.colors[ci(i)] for i in B_interp]
    ax.scatter(dist, B_interp, s=15, c=colors, edgecolors=edgecolors)

    #label
    ax.set_ylabel('Magnetic Field ($mG$)')
    sf = subcalc_funks._sig_figs
    dist = subcalc_funks.cumulative_distance(points)
    x, y = zip(*points)
    point_strings = [str((sf(x[i],3), sf(y[i],3))) for i in range(len(points))]
    if(x_labeling == 'location'):
        ax.set_xticks(dist)
        ax.set_xticklabels(point_strings, rotation=45, ha='right')
        ax.set_xlabel('X, Y Location ($ft$)')
    else:
        ax.set_xlabel('Distance Along Path ($ft$)')
    #ax.set_title('Magnetic Field along %s ft Path from %s' %
    #        (str(sf(dist[-1], 3)), ' to '.join(point_strings)))

    #format
    ax.margins(0.025)
    ax.grid(b=True)
    _format_ax(ax)
    #tight layout
    fig.tight_layout(pad=3)
    #write Bkey
    _write_Bkey(fig, ax, Bkey)
    #save or not
    _save_fig('segment-plot', fig, **kw)

    return(fig, ax)

def plot_cross_sections(obj, paths, xs_label_size=12, xs_color='black',
        map_style='contour', scale='lin', n=101, x_labeling='distance',
        show_max=True, fp_max=True, fp_text=True, cmap='viridis_r',
        max_fig_width=12, max_fig_height=8, legend_padding=6, **kw):
    """Generate a map style plot (either contour or pcolormesh) with
    cross sections labeled on it and generate plots of the fields corresponding
    to the cross sections
    args:
        obj - Results object or Model object, because the grid of fields
              in a Results object must be interpolated to compute fields
              along each cross section, passing a Model object instead will
              yield smoother profiles of the fields along each cross section.
              The Model object will, however, be used to compute a grid of
              fields for the map style plot (contour or pcolormesh), so
              passing a Model object could be slower.
        paths - An iterable of iterables of x,y pairs representing paths through
                the results domain to plot as cross sections. For example,
                    ([(1,2), (3,5)], [(2,5), (9,3), (4,7)], [(5,3), (9,2)])
    optional args:
        xs_label_size - int, fontsize of text labels on the map style figure
        xs_color - any matplotlib compatible color definition
        map_style - str, 'contour' or 'pcolormesh', determines which map style
                    plot is generated with the cross sections labeled on it,
                    default is 'contour'
        scale - str, can be 'log' or 'lin', only applies if map_style is
                'contour' (default is 'lin')
        n - integer, number of points sampled along the sections (default 101)
        x_labeling - 'distance' or 'location', for x axis ticks on the cross
                      section plots labeled according to the distance along
                      the segment or with the (x,y) coordinates of sample
                      points, default is 'distance'
        show_max - bool, toggle labeling of the maximum field location,
                    default is True
        fp_max - bool, toggle whether the maximum fields along footprints
                 "of concern" are labeled
        fp_text - bool, toggle footprint group labeling with text
        cmap - str, name of matplotlib colormap, see:
               http://matplotlib.org/examples/color/colormaps_reference.html
        max_fig_width - float/int, inches, maximum width of figure
        max_fig_height - float/int, inches, maximum height of figure
        legend_padding - float/int, inches, width left for legend area
    kw:
        prefix - string prepended to the file names of saved plots
        suffix - string appended to the file names of saved plots

            and

        any keyword arguments that can be passed to plot_contour(),
        plot_pcolormesh(), or plot_segment()

        note: Only a directory name can be passed to the 'path' keyword.
              File names aren't accepted, which prevents saved plots from
              overwriting each other. File names are created automatically.
    returns:
        A tuple of tuples of plotting objects. The first tuple contains the
        return arguments of the map plot (contour or pcolormesh) and the
        other tuples contain the return arguments of plot_path, for however
        many cross sections are created."""

    #deal with Model vs Results input
    if(type(obj) is subcalc_class.Results):
        res = obj
    elif(type(obj) is subcalc_class.Model):
        res = obj.calculate()

    #separate saving kw from others
    save_kw = {}
    for k in ['save', 'path', 'format', 'prefix', 'suffix']:
        if(k in kw):
            save_kw[k] = kw[k]
            kw.pop(k)

    #deal with the saving kw
    if('prefix' in save_kw):
        save_kw['save'] = True
        fn_prefix = save_kw['prefix']
        if(fn_prefix[-1] != '-'):
            fn_prefix = fn_prefix + '-'
    else:
        fn_prefix = ''

    if('suffix' in save_kw):
        save_kw['save'] = True
        fn_suffix = save_kw['suffix']
        if(fn_suffix[0] != '-'):
            fn_suffix = '-' + fn_suffix
    else:
        fn_suffix = ''

    if('save' in save_kw):
        save = save_kw['save']
    elif('path' in save_kw):
        if(not os.path.isdir(save_kw['path'])):
            raise(subcalc_class.EMFError('The path keyword argument to plot_cross_sections must be a directory path. Plot names are created automatically, with some control available through the prefix and suffix keyword arguments.'))
        save_kw['save'] = True
        save = True
    else:
        save = False

    #check inputs
    if(len(paths) > 26):
        raise(subcalc_class.EMFError('There cannot be more than 26 cross sections on a single figure. Make sure that your input for the "points" argument has the correct number of levels (sublists).'))

    #list of return arguments
    R = []

    #plot the map style figure
    if(map_style == 'contour'):
        r = plot_contour(res, scale, show_max, fp_max, fp_text, cmap,
                max_fig_width, max_fig_height, legend_padding, **kw)
        R.append(r)
        fig, ax = r[0], r[1]
        fn = fn_prefix + 'contour-with-cross-sections' + fn_suffix
    else:
        r = plot_pcolormesh(res, show_max, fp_max, fp_text, cmap, max_fig_width,
                max_fig_height, legend_padding, **kw)
        R.append(r)
        fig, ax = r[0], r[1]
        fn = fn_prefix + 'pcolormesh-with-cross-sections' + fn_suffix
    #draw cross section traces on the figure
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for i, path in enumerate(paths):
        #get x,y
        x, y = zip(*path)
        xb, xe, yb, ye = x[0], x[-1], y[0], y[-1]
        #plot the trace
        ax.plot(x, y, color=xs_color)
        #label the trace
        hab, vab = _get_text_alignment(path[1], path[0])
        hae, vae = _get_text_alignment(path[-2], path[-1])
        ax.text(xb, yb, alphabet[i], ha=hab, va=vab,
                color=xs_color, fontsize=xs_label_size)
        ax.text(xe, ye, alphabet[i] + "'", ha=hae, va=vae,
                color=xs_color, fontsize=xs_label_size)
    #save or don't
    if(save):
        _save_fig(fn, fig, **save_kw)

    #plot the cross sections
    for i, path in enumerate(paths):
        r = plot_path(obj, path, n, x_labeling, scale, cmap, **kw)
        R.append(r)
        fig, ax = r
        c = alphabet[i]
        ax.set_title("Cross Section %s-%s'" % (c, c))
        if(save):
            fn = '%scross-section%s' % (fn_prefix, c + fn_suffix)
            _save_fig(fn, fig, **save_kw)

    return(tuple(R))

def plot_wires_3D(mod, include_fields=False, cmap='viridis_r', Bkey='Bmax'):
    """Create a 3D plot showing wire segments in space.
    args:
        mod - Model object
    optional args:
        include_fields - bool, if True, a contour slice is drawn at the
                         modeled grid height
        cmap - string, colormap used to draw contours
               http://matplotlib.org/examples/color/colormaps_reference.html
        Bkey - if include_fields is True, Bkey can be used to specify which
               component of the fields is used to draw contours
    returns:
        fig - Figure object
        ax - matplotlib Axes3D object"""

    #make the plot
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = [], [], []
    for seg in mod.segments:
        x.append(seg[0][0])
        x.append(seg[1][0])
        y.append(seg[0][1])
        y.append(seg[1][1])
        z.append(seg[0][2])
        z.append(seg[1][2])
        ax.plot(x[-2:], y[-2:], z[-2:], 'k')
    #set limits to span equal ranges
    xlim = min(x), max(x)
    ylim = min(y), max(y)
    zlim = min(z), max(z)
    xr = xlim[1] - xlim[0]
    yr = ylim[1] - ylim[0]
    zr = zlim[1] - zlim[0]
    xm = sum(xlim)/2.0
    ym = sum(ylim)/2.0
    zm = sum(zlim)/2.0
    hmr = max([xr, yr, zr])/2.0
    ax.set_xlim(xm - hmr, xm + hmr)
    ax.set_ylim(ym - hmr, ym + hmr)
    ax.set_zlim(zm - hmr, zm + hmr)
    #add text
    ax.set_xlabel('x (ft)')
    ax.set_ylabel('y (ft)')
    ax.set_zlabel('z (ft)')
    ax.axis('equal')
    ax.set_title(mod.name)

    if(include_fields):
        res = mod.calculate()
        res.Bkey = Bkey
        ax.contour(res.X, res.Y, res.B, offset=mod.z, cmap=cmap)
        ax.set_title(mod.name + '-' + Bkey)

    return(fig, ax)
