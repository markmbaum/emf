from .. import np, mpl, plt, textwrap

from ..emf_plots import _save_fig, _prepare_fig

import fields_funks

#rcparams for more static global formatting changes
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.figsize'] = (14, 6)
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.color'] = (.2, .2, .2)
mpl.rcParams['axes.labelcolor'] = (.2, .2, .2)
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['legend.borderaxespad'] = 0 #mpl default is None
mpl.rcParams['xtick.color'] = (.2, .2, .2)
mpl.rcParams['ytick.color'] = (.2, .2, .2)

#other more specific/dynamic global formatting variables
_B_color = 'darkgreen'
_E_color = 'midnightblue'
_fields_linewidth = 1.75
_ROW_linewidth = 0.75
_ROW_color = 'gray'
_ground_surface_linewidth = 1
_ground_surface_color = 'gray'
_ax_frameon = False
_ax_ticks_on = False
_leg_edge_on = False
_colormap = [(0, 0.4470, 0.7410),(0.8500, 0.3250,0.0980),
            (0.9290, 0.6940, 0.1250),(0.4940, 0.1840, 0.5560),
            (0.4660, 0.6740, 0.1880),(0.3010, 0.7450, 0.9330),
            (0.6350, 0.0780, 0.1840)]
#useful globals for the CrossSection plotting routines, unlikely to collide
#with other variables of the same name
_fields_plots_xs_headspace = 0.45 #space at the top of plots for legend
_include_headspace = True #toggle scaling to make legend overlap unlikely
_fields_plots_xs_wireperc = 0.3 #percent of max field value to scale wire height

#-------------------------------------------------------------------------------
#general plotting support functions

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

def _format_bar_axes_legends(*args):
    """Apply axis formatting commands to axes objects
    args:
        some number of axes objects"""
    for i in range(len(args)):
        ax = args[i]
        #apply legend formatting
        leg = ax.get_legend()
        if(leg):
            rec = leg.get_frame()
            if(not _leg_edge_on):
                rec.set_edgecolor('white')
        #apply axis formatting
        if(i == 0):
            _spines_off(ax, 'top', 'right')
        else:
            _spines_off(ax, 'top', 'right')
        if(not _ax_ticks_on):
            ax.tick_params(axis = 'both', which = 'both',
                bottom = 'off', top = 'off', left = 'off', right = 'off')
    #make limits equal
    minylim, maxylim = args[0].get_ylim()
    for ax in args:
        yl = ax.get_ylim()
        if(yl[0] < minylim):
            minylim = yl[0]
        if(yl[1] > maxylim):
            maxylim = yl[1]
    for ax in args:
        ax.set_ylim(minylim, maxylim)

def _format_line_axes_legends(*args):
    """Apply axis formatting commands to axes objects
    args:
        some number of axes objects with twin x axes"""
    for ax in args:
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

def _format_twin_axes(*args):
    """take care of scaling problems"""
    if(len(args) > 1):
        #get minimum y limit
        ylow, yhigh = 0., 0.
        for ax in args:
            yl = ax.get_ylim()
            if(yl[0] < ylow):
                ylow = yl[0]
                yhigh = yl[1]
        #scale all axes identically so that they overlap at y = 0
        if(yhigh != 0):
            frac = ylow/yhigh
            for ax in args:
                yl = ax.get_ylim()
                ax.set_ylim(frac*yl[1], yl[1])


def _check_und_conds(xss, axs, **kw):
    """Check for underground conductors, adding some padding at the bottom
    of the axes and drawing a ground surface line if there are underground
    lines, or setting the lower y axis limit to zero if not
    args:
        xss - list of CrossSection objects containing Conductors to check
        axs - list of Axes object to potentially apply padding to"""
    und = []
    for xs in xss:
        und += [(i.y <= 0) for i in xs.conds]
    if(any(und)):
        #draw ground surface line
        xl = axs[0].get_xlim()
        kw['H'].append(axs[0].plot(xl, [0]*2, ':', color=_ground_surface_color,
                linewidth=_ground_surface_linewidth)[0])
        kw['L'].append('Ground Surface')
        for ax in axs:
            #apply padding
            yl = ax.get_ylim()
            r = yl[1] - yl[0]
            ax.set_ylim(yl[0] - r*0.05, yl[1])
    else:
        for ax in axs:
            yl = ax.get_ylim()
            ax.set_ylim(0, yl[1])

def _spines_on(ax, *args):
    """Turn axes spines on by passing strings indicating which spines to
    make visible, like 'left', 'bottom', 'top', 'right'
    args:
        ax - target Axes object
        strings indicating spines to turn on"""
    for s in args:
        ax.spines[s].set_visible(True)

def _spines_off(ax, *args):
    """Turn axes spines on by passing strings indicating which spines to
    make invisible, like 'left', 'bottom', 'top', 'right'
    args:
        ax - target Axes object
        strings indicating spines to turn off"""
    for s in args:
        ax.spines[s].set_visible(False)

def _color_twin_axes(ax1, color1, ax2, color2):
    """Assign colors to split y axes"""
    #spines
    ax1.spines['left'].set_color(color1)
    ax1.spines['right'].set_color(color2)
    ax2.spines['left'].set_color(color1)
    ax2.spines['right'].set_color(color2)
    #text
    ax1.yaxis.label.set_color(color1)
    ax2.yaxis.label.set_color(color2)
    #ticks
    ax1.tick_params(axis = 'y', colors = color1)
    ax2.tick_params(axis = 'y', colors = color2)

#-------------------------------------------------------------------------------
#plotting routines working with a CrossSection object

def _find_xmax(xs, **kw):
    #get x cutoff, if any
    if('xmax' in kw):
        xmax = abs(kw['xmax'])
    else:
        xmax = max(abs(xs.fields.index))
    return(xmax)

def _plot_wires(ax, hot, gnd, v, **kw):
    """Plot conductor symbols in ax. Returns handles for the hot and gnd
    Conductors.
    args:
        ax - target axis
        hot - list of non-grounded conductors
        gnd - list of grounded conductors
        v - iterable of calculated field results used to scale conductors
    kw:
        scale - if False, will not scale conductor heights to false units"""
    #get x and y coordinates
    L = len(hot)
    x = np.array([c.x for c in hot + gnd])
    y = np.array([c.y for c in hot + gnd])
    #calculate the scaling factor
    scale = _fields_plots_xs_wireperc*max(np.absolute(v))/max(np.absolute(y))
    if('scale' in kw):
        if(kw['scale'] is False):
            scale = 1.0
    #plot
    if(hot):
        kw['H'].append(ax.plot(x[:L], scale*y[:L], 'kd')[0])
        kw['L'].append('Conductors')
    if(gnd):
        kw['H'].append(ax.plot(x[L:], scale*y[L:], 'd', color = 'gray')[0])
        kw['L'].append('Grounded Conductors')

def _plot_ROW_edges(ax, lROW, rROW, **kw):
    """Plot dashed lines marking the left and right edges of the
    Right-of-Way. Axis limits are also adjusted to allow extra space on the
    sides to make the ROW edge lines visible if needed. The ROW edge line
    handles are returned in a list.
    args:
        ax - target axis
        lROW - x value of the left edge of the ROW
        rROW - x value of the right edge of the ROW"""
    yl = list(ax.get_ylim())
    if(yl[0] > 0):
        yl[0] = 0
    kw['H'].append(ax.plot([lROW]*2, yl, '--', color = _ROW_color,
                        linewidth = _ROW_linewidth, zorder = -1)[0])
    ax.plot([rROW]*2, yl, '--', color = _ROW_color,
                linewidth = _ROW_linewidth, zorder = -1)
    kw['L'].append('ROW Edges')
    xl = ax.get_xlim()
    if((xl[0] == lROW) or (xl[1] == rROW)):
        ax.set_xlim((xl[0]*1.15, xl[1]*1.15))

def plot_Bmax(xs, **kw):
    """Plot the maximum magnetic field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xs - a CrossSection object
    kw:
        ax - target Axes
        fig - Figure object, target figure for plotting, overridden by 'ax'
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - Figure object
        ax - Axes object"""

    #get plotting objects
    fig, ax = _prepare_fig(**kw)
    #get plotting specs
    xmax = _find_xmax(xs)
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    #plot the field curve
    kw['H'].append(ax.plot(xs.fields['Bmax'][-xmax:xmax],
                color=_B_color, linewidth=_fields_linewidth)[0])
    kw['L'].append(r'Magnetic Field $(mG)$')
    #plot wires
    _plot_wires(ax, xs.hot, xs.gnd, xs.fields['Bmax'], **kw)
    _check_und_conds([xs], [ax], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xs.lROW, xs.rROW, **kw)
    #set axis text and legend
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Maximum Magnetic Field $(mG)$')
    ax.set_title('Maximum Magnetic Field - %s' % xs.title)
    ax.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax)
    #save the fig or don't, depending on keywords
    _save_fig(xs.sheet, fig, **kw)
    #return
    return(fig, ax)

def plot_Emax(xs, **kw):
    """Plot the maximum electric field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xs - a CrossSection object
    kw:
        ax - target Axes
        fig - Figure object, target figure for plotting, overridden by 'ax'
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - Figure object"""

    #get plotting objects
    fig, ax = _prepare_fig(**kw)
    #get plotting specs
    xmax = _find_xmax(xs)
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    #plot the field curve
    kw['H'].append(ax.plot(xs.fields['Emax'][-xmax:xmax], color=_E_color,
                linewidth=_fields_linewidth)[0])
    kw['L'].append(r'Electric Field $(kV/m)$')
    #plot wires
    _plot_wires(ax, xs.hot, xs.gnd, xs.fields['Emax'], **kw)
    _check_und_conds([xs], [ax], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xs.lROW, xs.rROW, **kw)
    #set axis text and legend
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Maximum Electric Field $(kV/m)$')
    ax.set_title('Maximum Electric Field - %s' % xs.title)
    ax.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax)
    #save the fig or don't, depending on keywords
    _save_fig(xs.sheet, fig, **kw)
    #return
    return(fig, ax)

def plot_max_fields(xs, **kw):
    """Plot the maximum magnetic and electric field on split vertical axes,
    along the ROW with conductor locations shown in artificial coordinates.
    args:
        xs - a CrossSection object
    kw:
        ax - target Axes
        fig - Figure object, target figure for plotting, overridden by 'ax'
        title - string, exact plot title
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
    returns:
        fig - Figure object
        ax_B - B field Axes
        ax_E - E field Axes"""

    #get plotting objects
    fig, ax_B = _prepare_fig(**kw)
    ax_E = ax_B.twinx()
    #get plotting specs
    xmax = _find_xmax(xs)
    #plot the field curves
    kw['H'] = [ax_B.plot(xs.fields['Bmax'][-xmax:xmax], color=_B_color,
                        linewidth=_fields_linewidth)[0],
                ax_E.plot(xs.fields['Emax'][-xmax:xmax], color=_E_color,
                        linewidth=_fields_linewidth)[0]]
    kw['L'] = [r'Magnetic Field $(mG)$',
            r'Electric Field $(kV/m)$']
    #plot wires
    _plot_wires(ax_B, xs.hot, xs.gnd, xs.fields['Bmax'], **kw)
    _check_und_conds([xs], [ax_B, ax_E], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax_B.get_ylim()
        ax_B.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
        yl = ax_E.get_ylim()
        ax_E.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax_B, xs.lROW, xs.rROW, **kw)
    #set axis text
    ax_B.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax_B.set_ylabel(r'Maximum Magnetic Field $(mG)$', color=_B_color)
    ax_E.set_ylabel(r'Maximum Electric Field $(kV/m)$', color=_E_color)
    ax_B.set_title('Maximum Magnetic and Electric Fields - %s' % xs.title)
    #set color of axis spines and ticklabels
    _color_twin_axes(ax_B, _B_color, ax_E, _E_color)
    #legend
    ax_B.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax_B, ax_E)
    _format_twin_axes(ax_B, ax_E)
    #save the fig or don't, depending on keywords
    _save_fig(xs.sheet, fig, **kw)
    #return
    return(fig, ax_B, ax_E)

def plot_xs(xs, **kw):
    """Plot CrossSection conductors and ROW edges without any fields
    args:
        xs - CrossSection objects
    kw:
        ax - target Axes
        fig - Figure object, target figure for plotting, overridden by 'ax'
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""
    #get plotting objects
    fig, ax = _prepare_fig(**kw)
    size = fig.get_size_inches()
    fig.set_size_inches(size[0] + 3, size[1])
    #init handles and labels lists for legend
    kw['H'], kw['L'], kw['scale'] = [], [], False
    #plot wires
    _plot_wires(ax, xs.hot, xs.gnd, [c.y for c in xs.conds], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xs.lROW, xs.rROW, **kw)
    #check underground conductors
    _check_und_conds([xs], [ax], **kw)
    #set axis text and legend
    ax.set_title('Cross Section Configuration - %s' % xs.title)
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Height Above Ground $(ft)$')
    ax.legend(kw['H'], kw['L'], numpoints=1)
    #format
    _format_line_axes_legends(ax)
    #save or don't
    _save_fig(xs.sheet, fig, **kw)
    #return
    return(fig, ax)

def _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, field, unit, **kw):
    """Handle plotting of DAT comparison features that don't require unique
    strings
    args:
        ax_abs - axis of absolute error plot
        ax_per - axis of percentage error plot
        ax_mag - axis of field magnitude plot
        pan - pandas.Panel object containing results and errors
        field - string, column label of field to be plotted (Bmax/Emax)"""

    #plot absolute error
    h_abs = ax_abs.plot(pan['Absolute Difference'][field],
            color=mpl.rcParams['axes.labelcolor'], zorder=-2)
    ax_abs.set_ylabel('Absolute Difference' + unit)
    #plot percentage error
    h_per = ax_per.plot(pan['Percent Difference'][field],
            color='firebrick', zorder=-1)
    ax_per.set_ylabel('Percent Difference', color='firebrick')
    #set error axes legend
    ax_per.legend(h_abs + h_per, ['Absolute Difference','Percent Difference'])
    ax_per.get_legend().set_zorder(1)
    #plot full results profiles
    kw['H'] += [ax_mag.plot(pan['FIELDS_DAT_results'][field],
                    color=_colormap[1])[0],
                ax_mag.plot(pan['python_results'][field],
                    color=_colormap[0])[0]]
    kw['L'] += ['FIELDS', 'New Code']
    ax_mag.set_xlabel('Distance from ROW Center $(ft)$')

def _plot_DAT_comparison(xs, pan, **kw):
    """Generate 2 subplots showing the FIELDS results (from a .DAT file)
    compared to the results of this code and the error.
    args:
        xs - CrossSection object
        pan - pandas.Panel object containing results and errors
    kw:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Bmax
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Bmax', '$(mG)$', **kw)
    _plot_wires(ax_mag, xs.hot, xs.gnd, pan['python_results']['Bmax'], **kw)
    _check_und_conds([xs], [ax_mag], **kw)
    ax_abs.set_title('Absolute and Percent Difference, Max Magnetic Field')
    ax_mag.set_ylabel('Bmax $(mG)$')
    ax_mag.set_title('Model Results, Magnetic Field')
    ax_mag.legend(kw['H'], kw['L'], numpoints=1)
    _color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'firebrick')
    _format_line_axes_legends(ax_abs, ax_per, ax_mag)
    _format_twin_axes(ax_abs, ax_per)
    _save_fig(xs.sheet + '-DAT_comparison_Bmax', fig, **kw)
    plt.close(fig)

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Emax
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Emax', '$(kV/m)$', **kw)
    _plot_wires(ax_mag, xs.hot, xs.gnd, pan['python_results']['Emax'], **kw)
    _check_und_conds([xs], [ax_mag], **kw)
    ax_abs.set_title('Absolute and Percent Difference, Max Electric Field')
    ax_mag.set_ylabel('Emax $(kV/m)$')
    ax_mag.set_title('Model Results, Electric Field')
    ax_mag.legend(kw['H'], kw['L'], numpoints=1)
    _color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'firebrick')
    _format_line_axes_legends(ax_abs, ax_per, ax_mag)
    _format_twin_axes(ax_abs, ax_per)
    _save_fig(xs.sheet + '-DAT_comparison_Emax', fig, **kw)
    plt.close(fig)

#-------------------------------------------------------------------------------
#plotting routines working primarily with a SectionBook object

#useful globals for the section book plotting routine(s), unlikely to collide
#with other variables of the same name
_fields_plots_sb_headspace = 0.4 #space at the top of plots for legend
_fields_plots_sb_wireperc = 0.3 #percent of max field value to scale wire heights

def _plot_group_fields(ax, xss, field, **kw):
    """Plot the results of fields calculations
    args:
        ax - target axis
        xss - list of CrossSection objects to plot results from
        field - column label of the results to plot
    kw:
        xmax - cutoff distance from ROW center"""
    #check for an xmax keyword
    if('xmax' in kw):
        xmax = kw['xmax']
    else:
        xmax = False
    #plot the fields, keeping handles and finding max of all
    fields_list = [xs.fields[field] for xs in xss]
    for i in range(len(fields_list)):
        #plot
        if(xmax):
            kw['H'].append(ax.plot(fields_list[i][-xmax:xmax],
                    color = _colormap[i%7],
                    linewidth = _fields_linewidth)[0])
        else:
            kw['H'].append(ax.plot(fields_list[i], color = _colormap[i%7],
                    linewidth = _fields_linewidth)[0])
        kw['L'].append(xss[i].sheet)

def _plot_group_wires(ax, xss, max_field, **kw):
    """Plot the conductors of 1 or 2 CrossSections, using split color
    for conductor locations shared by 2 CrossSections
    args:
        ax - target axis
        xss - list of CrossSection objects to plot results from
        max_field - maximum value of all fields plotted in ax"""

    #---only hot conductors are plotted
    #---use sets to group shared and unshared conductor locations
    #check if all xss have the same conductor locations
    all_same = True
    xy_0 = set([(c.x,c.y) for c in xss[0].hot])
    for i in range(1, len(xss)):
        xy_1 = set([(c.x,c.y) for c in xss[i].hot])
        if(xy_0 != xy_1):
            all_same = False
            break
        else:
            xy_0 = xy_1
    if(all_same):
        x, y = [], []
        for i in xy_0:
            x.append(i[0])
            y.append(i[1])
        y = np.array(y)
        #plot conductors
        scale = _fields_plots_xs_wireperc*max_field/np.max(np.absolute(y))
        kw['H'].append(ax.plot(x, scale*y, 'd', color = 'k')[0])
        kw['L'].append('Conductors')
    elif(len(xss) == 2):
        #get x and y pairs of Conductors in each group
        xy_0 = [(c.x,c.y) for c in xss[0].hot]
        xy_1 = [(c.x,c.y) for c in xss[1].hot]
        #grab all x and y coordinates while they're available
        all_x, all_y = zip(*(xy_0 + xy_1))
        all_y = np.array(all_y)
        all_x = np.array(all_x)
        #use sets to form the groups of Conductors
        xy_0 = set(xy_0)
        xy_1 = set(xy_1)
        #assemble shared x,y pairs
        shared = xy_0 & xy_1
        #remove shared pairs from xy_0 and xy_1
        xy_0 -= shared
        xy_1 -= shared
        scale = _fields_plots_xs_wireperc*max_field/np.max(np.absolute(all_y))
        #cross section 0 conductors only
        if(len(xy_0) > 0):
            x,y = zip(*xy_0)
            kw['H'].append(ax.plot(x, scale*np.array(y), 'd',
                            color = _colormap[0])[0])
        else:
            #still need a handle for the legend
            kw['H'].append(mpl.lines.Line2D([], [], marker = 'd',
                    linestyle = '', color = _colormap[0]))
        kw['L'].append('Conductors - ' + xss[0].sheet)
        #cross section 1 conductors only
        if(len(xy_1) > 0):
            x,y = zip(*xy_1)
            kw['H'].append(ax.plot(x, scale*np.array(y), 'd',
                            color = _colormap[1])[0])
        else:
            #still need a handle for the legend
            kw['H'].append(mpl.lines.Line2D([], [], marker = 'd',
                    linestyle = '', color = _colormap[1]))
        kw['L'].append('Conductors - ' + xss[1].sheet)
        #shared conductors
        if(len(shared) > 0):
            x,y = zip(*shared)
            ax.plot(x, scale*np.array(y), 'd', color = _colormap[0],
                        fillstyle = 'left')
            ax.plot(x, scale*np.array(y), 'd', color = _colormap[1],
                        fillstyle = 'right')
        #underground line adjustments
        _check_und_conds(xss, [ax], **kw)

def _plot_group_ROW_edges(ax, xss, **kw):
    """Plot lines delineating the Right-of-Way boundaries
    args:
        ax - target axis
        xss - list of CrossSection objects to plot results from"""

    #if the CrossSections all have the same ROW edges, plot a single pair
    l, r = [xs.lROW for xs in xss], [xs.rROW for xs in xss]
    if(all([i == l[0] for i in l[1:]]) and all([i == r[0] for i in r[1:]])):
        _plot_ROW_edges(ax, l[0], r[0], **kw)

    #if there are only two CrossSections and they have different ROW edges,
    #handle ROW edges independently
    elif(len(xss) == 2):
        yl = ax.get_ylim()
        ax.plot([xss[0].lROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1)
        ax.plot([xss[1].lROW]*2, yl, '--', color = _colormap[1],
                    linewidth = _ROW_linewidth, zorder = -1)
        #check if the left ROW edges are in the same place for both sections
        ax.plot([xss[0].rROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1)
        ax.plot([xss[1].rROW]*2, yl, '--', color = _colormap[1],
                    linewidth = _ROW_linewidth, zorder = -1)
        #create a line handle and label for each color of the dashed line
        kw['H'] += [
                mpl.lines.Line2D([], [], linestyle = '--',color = _colormap[0]),
                mpl.lines.Line2D([], [], linestyle = '--', color = _colormap[1])
                ]
        kw['L'] += ['ROW Edges - ' + xss[0].sheet, 'ROW Edges - ' + xss[1].sheet]

def plot_groups(sb, **kw):
    """Plot the fields of grouped CrossSections in the same axis, a plot for
    both fields. Only plots groups of more than one CrossSection.
    args:
        sb - SectionBook object to pull plotting groups from
    kw:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
        xmax - cutoff distance from ROW center
        B - bool, toggle magnetic field plots, default is True
        E - bool, toggle electric field plots, default is True
        groups - a list of group names to plot, default is all groups
        return_figs - toggle whether a list of figure objects is returned
                      instead of closing the figures to clear memory,
                      default is False
    returns:
        figs - dict of dicts, keys are 'E' and 'B', which are each keyed
                by group tags, leading to Figure objects"""

    #check kws
    return_figs = False
    if('return_figs' in kw):
        if(kw['return_figs']):
            return_figs = True
            figs = {'E': {}, 'B': {}}
    B_flag = True
    if('B' in kw):
        B_flag = kw['B']
    E_flag = True
    if('E' in kw):
        E_flag = kw['E']
    groups = sb.tags
    if('groups' in kw):
        groups = set(kw['groups'])

    #iterate over groups with more than 1 CrossSection
    for xss in sb.tag_groups:

        if(xss[0].tag in groups):

            #BMAX
            if(B_flag):
                #get plotting objects
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                #init handles and labels lists for legend
                kw['H'], kw['L'] = [], []
                #plot the Bmax results for each xs in the group
                _plot_group_fields(ax, xss, 'Bmax', **kw)
                #plot wires
                max_field = max([xs.fields['Bmax'].max() for xs in xss])
                _plot_group_wires(ax, xss, max_field, **kw)
                #adjust axis limits if called for
                if(_include_headspace):
                    yl = ax.get_ylim()
                    ax.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
                #plot ROW lines
                _plot_group_ROW_edges(ax, xss, **kw)
                #set axis text and legend
                ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
                ax.set_ylabel(r'Maximum Magnetic Field $(mG)$')
                t = 'Maximum Magnetic Field - %s' % str(xss[0].tag)
                ax.set_title(t)
                ax.legend(kw['H'], kw['L'], numpoints=1)
                _format_line_axes_legends(ax)
                #save the figure if keyword 'save' == True, and append fig
                _save_fig('group_%s-Bmax' % str(xss[0].tag), fig, **kw)
                #store the fig or close it
                if(return_figs):
                    figs['B'][xss[0].tag] = fig
                else:
                    plt.close(fig)

            #EMAX
            if(E_flag):
                #get plotting objects
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                #init handles and labels lists for legend
                kw['H'], kw['L'] = [], []
                #plot the Bmax results for each xs in the group
                _plot_group_fields(ax, xss, 'Emax', **kw)
                #plot wires
                max_field = max([xs.fields['Emax'].max() for xs in xss])
                _plot_group_wires(ax, xss, max_field, **kw)
                #adjust axis limits if called for
                if(_include_headspace):
                    yl = ax.get_ylim()
                    ax.set_ylim(yl[0], (1 + _fields_plots_xs_headspace)*yl[1])
                #plot ROW lines
                _plot_group_ROW_edges(ax, xss, **kw)
                #set axis text and legendf
                ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
                ax.set_ylabel(r'Maximum Electric Field $(kV/m)$')
                t = 'Maximum Electric Field - %s' % str(xss[0].tag)
                ax.set_title(t)
                ax.legend(kw['H'], kw['L'], numpoints=1)
                _format_line_axes_legends(ax)

                #save the figure if keyword 'save' == True or a path string is passed
                _save_fig('group_%s-Emax' % str(xss[0].tag), fig, **kw)
                #store the fig or close it
                if(return_figs):
                    figs['E'][xss[0].tag] = fig
                else:
                    plt.close(fig)

    if(return_figs):
        return(figs)

def _reorder_xss(xss, **kw):
    """Use the xs_order kw, if it exists, to create a list of
    CrossSection sheet strings that determines the order bars are plotted
    in
    args:
        xss - list of CrossSections corresponding to a group
    kw:
        xs_order - dict, if any keys are the same as the tag of the
                   CrossSection group represented by xss, they should
                   map to a list of strings with CrossSection sheet names
                   specifiying an order to plot in.
    returns:
        reorder - a reordered list of CrossSections"""

    #see if there is an xs_order kwarg
    if('xs_order' in kw):
        xs_order = kw['xs_order']
        tag = xss[0].tag
        #see if it corresponds to the group represented by xss
        if(tag in xs_order):
            #get the order
            order = xs_order[tag]
            #compile xs sheet names
            sheets = [xs.sheet for xs in xss]
            #make a dict for later
            d = dict(zip(sheets, xss))
            #reorder the sheets
            reorder = []
            for sh in order:
                if(sh in sheets):
                    reorder.append(sheets.pop(sheets.index(sh)))
            #put leftovers at the end of the list
            reorder = reorder + sheets
            #generate reodered list and return
            return([d[k] for k in reorder])

    #return the order of xss
    return(xss)

def _plot_group_bars(ax, xss, field, side):
    """Plot bars in ax for each xs in xss
    args:
        ax - target Axes object
        xss - list of CrossSection objects corresponding to a group
        field - string, field component to plot (Ex, Bx, Bprod, etc.)
        side - string, 'left' or 'right', selects which ROW edge side"""

    #translate side to index
    if(side == 'left'):
        side = 0
    else:
        side = 1
    #plot the bars
    x = range(len(xss))
    values = [xs.ROW_edge_fields[field].values[side] for xs in xss]
    ax.bar(x, values, color=[_colormap[i%7] for i in range(len(xss))],
            bottom=0.0, align='center', alpha=0.8, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([textwrap.fill(xs.sheet, 15) for xs in xss],
            rotation='vertical', fontsize=10)

def _generate_ROW_value_plot_objects(xss):
    """generate figure and axes for plot_groups_at_ROW
    args:
        xss - list of CrossSections in group
    returns:
        fig - figure object
        axl - left axes
        axr - right axes"""

    #see if ROW edge distances are unique
    lROW, rROW = set([xs.lROW for xs in xss]), set([xs.rROW for xs in xss])
    if(len(lROW) == 1):
        lROW_title_add = (r' (%s $ft$)'
                    % str(fields_funks._sig_figs(lROW.pop(), 3)))
    else:
        lROW_title_add = ''
    if(len(rROW) == 1):
        rROW_title_add = (r' (%s $ft$)'
                    % str(fields_funks._sig_figs(rROW.pop(), 3)))
    else:
        rROW_title_add = ''

    #scaling parameters
    ax_y_up_frac = 0.9
    ax_y_shrink_frac = 0.85

    #get plotting objects
    fig = plt.figure()
    #left ROW axis
    axl = fig.add_subplot(1,2,1)
    axl.set_title('Left ROW Edge' + lROW_title_add)
    pos = axl.get_position()
    axl.set_position([pos.x0,
                    pos.y0 + pos.height*(1 - ax_y_up_frac),
                    pos.width,
                    pos.height*ax_y_shrink_frac])
    #right ROW axis
    axr = fig.add_subplot(1,2,2)
    axr.set_title('Right ROW Edge' + rROW_title_add)
    pos = axr.get_position()
    axr.set_position([pos.x0,
                    pos.y0 + pos.height*(1 - ax_y_up_frac),
                    pos.width,
                    pos.height*ax_y_shrink_frac])
    return(fig, axl, axr)

def plot_groups_at_ROW(sb, **kw):
    """Create bar charts showing the field values at ROW edges for each
    CrossSection group in a SectionBook
    args:
        sb - SectionBook object to pull plotting groups from
    kw:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
        B - bool, toggle magnetic field plots, default is True
        E - bool, toggle electric field plots, default is True
        groups - a list of group names to plot, default is all groups
        xs_order - dict, keys are CrossSection group tags, which map to
                   lists of CrossSection sheets specifying the order of the
                   plotted CrossSection bars (left to right). Not all
                   CrossSecitons in a group must be listed. Any/all can be
                   left out.
        return_figs - toggle whether a list of figure objects is returned
                      instead of closing the figures to clear memory,
                      default is False
    returns:
        figs - dict of dicts, keys are 'E' and 'B', which are each keyed
                by group tags, leading to Figure objects"""

    #check kws
    return_figs = False
    if('return_figs' in kw):
        if(kw['return_figs']):
            return_figs = True
            figs = {'E': {}, 'B': {}}
    B_flag = True
    if('B' in kw):
        B_flag = kw['B']
    E_flag = True
    if('E' in kw):
        E_flag = kw['E']
    groups = sb.tags
    if('groups' in kw):
        groups = set(kw['groups'])

    #iterate over groups with more than 1 CrossSection
    for xss in sb.tag_groups:

        if(xss[0].tag in groups):

            #get reordered CrossSection list
            xss = _reorder_xss(xss, **kw)

            #plot Bmax
            if(B_flag):
                fig, axl, axr = _generate_ROW_value_plot_objects(xss)
                _plot_group_bars(axl, xss, 'Bmax', 'left')
                _plot_group_bars(axr, xss, 'Bmax', 'right')
                #format
                _format_bar_axes_legends(axl, axr)
                #apply text
                axl.set_ylabel(r'Maximum Magnetic Field $(mG)$')
                fig.suptitle('Maximum Magnetic Fields at ROW Edges - %s' %
                        xss[0].tag, fontsize=18)
                #save?
                _save_fig('group_%s-ROW-Bmax' % str(xss[0].tag), fig, **kw)
                #store the fig or close it
                if(return_figs):
                    figs['B'][xss[0].tag] = fig
                else:
                    plt.close(fig)

            #plot Emax
            if(E_flag):
                fig, axl, axr = _generate_ROW_value_plot_objects(xss)
                _plot_group_bars(axl, xss, 'Emax', 'left')
                _plot_group_bars(axr, xss, 'Emax', 'right')
                #format
                _format_bar_axes_legends(axl, axr)
                #apply text
                axl.set_ylabel(r'Maximum Electric Field $(kV/m)$')
                fig.suptitle('Maximum Electric Fields at ROW Edges - %s' %
                        xss[0].tag, fontsize=18)
                #save?
                _save_fig('group_%s-ROW-Emax' % str(xss[0].tag), fig, **kw)
                #store the fig or close it
                if(return_figs):
                    figs['E'][xss[0].tag] = fig
                else:
                    plt.close(fig)

    if(return_figs):
        return(figs)
