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
_fields_linewidth = 2
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
_fields_plots_xc_headspace = 0.4 #space at the top of plots for legend
_include_headspace = True #toggle scaling to make legend overlap unlikely
_fields_plots_xc_wireperc = 0.3 #percent of max field value to scale wire height

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


def _check_und_conds(xcs, axs, **kw):
    """Check for underground conductors, adding some padding at the bottom
    of the axes and drawing a ground surface line if there are underground
    lines, or setting the lower y axis limit to zero if not
    args:
        xcs - list of CrossSection objects containing Conductors to check
        axs - list of Axes object to potentially apply padding to"""
    und = []
    for xc in xcs:
        und += [(i.y <= 0) for i in xc.conds]
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
            ax.set_ylim(yl[0] - r*0.05)
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

def _find_xmax(xc, **kw):
    #get x cutoff, if any
    if('xmax' in kw):
        xmax = abs(kw['xmax'])
    else:
        xmax = max(abs(xc.fields.index))
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
    scale = _fields_plots_xc_wireperc*max(np.absolute(v))/max(np.absolute(y))
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
    yl = ax.get_ylim()
    kw['H'].append(ax.plot([lROW]*2, yl, '--', color = _ROW_color,
                        linewidth = _ROW_linewidth, zorder = -1)[0])
    ax.plot([rROW]*2, yl, '--', color = _ROW_color,
                linewidth = _ROW_linewidth, zorder = -1)
    kw['L'].append('ROW Edges')
    xl = ax.get_xlim()
    if((xl[0] == lROW) or (xl[1] == rROW)):
        ax.set_xlim((xl[0]*1.15, xl[1]*1.15))

def plot_Bmax(xc, **kw):
    """Plot the maximum magnetic field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
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
    xmax = _find_xmax(xc)
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    #plot the field curve
    kw['H'].append(ax.plot(xc.fields['Bmax'][-xmax:xmax],
                color=_B_color, linewidth=_fields_linewidth)[0])
    kw['L'].append(r'Magnetic Field $(mG)$')
    #plot wires
    _plot_wires(ax, xc.hot, xc.gnd, xc.fields['Bmax'], **kw)
    _check_und_conds([xc], [ax], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xc.lROW, xc.rROW, **kw)
    #set axis text and legend
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Maximum Magnetic Field $(mG)$')
    ax.set_title('Maximum Magnetic Field - %s' % xc.title)
    ax.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax)
    #save the fig or don't, depending on keywords
    _save_fig(xc.sheet, fig, **kw)
    #return
    return(fig, ax)

def plot_Emax(xc, **kw):
    """Plot the maximum electric field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
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
    xmax = _find_xmax(xc)
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    #plot the field curve
    kw['H'].append(ax.plot(xc.fields['Emax'][-xmax:xmax], color=_E_color,
                linewidth=_fields_linewidth)[0])
    kw['L'].append(r'Electric Field $(kV/m)$')
    #plot wires
    _plot_wires(ax, xc.hot, xc.gnd, xc.fields['Emax'], **kw)
    _check_und_conds([xc], [ax], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xc.lROW, xc.rROW, **kw)
    #set axis text and legend
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Maximum Electric Field $(kV/m)$')
    ax.set_title('Maximum Electric Field - %s' % xc.title)
    ax.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax)
    #save the fig or don't, depending on keywords
    _save_fig(xc.sheet, fig, **kw)
    #return
    return(fig, ax)

def plot_max_fields(xc, **kw):
    """Plot the maximum magnetic and electric field on split vertical axes,
    along the ROW with conductor locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
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
    xmax = _find_xmax(xc)
    #plot the field curves
    kw['H'] = [ax_B.plot(xc.fields['Bmax'][-xmax:xmax], color=_B_color,
                        linewidth=_fields_linewidth)[0],
                ax_E.plot(xc.fields['Emax'][-xmax:xmax], color=_E_color,
                        linewidth=_fields_linewidth)[0]]
    kw['L'] = [r'Magnetic Field $(mG)$',
            r'Electric Field $(kV/m)$']
    #plot wires
    _plot_wires(ax_B, xc.hot, xc.gnd, xc.fields['Bmax'], **kw)
    _check_und_conds([xc], [ax_B, ax_E], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax_B.get_ylim()
        ax_B.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
        yl = ax_E.get_ylim()
        ax_E.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax_B, xc.lROW, xc.rROW, **kw)
    #set axis text
    ax_B.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax_B.set_ylabel(r'Maximum Magnetic Field $(mG)$', color=_B_color)
    ax_E.set_ylabel(r'Maximum Electric Field $(kV/m)$', color=_E_color)
    ax_B.set_title('Maximum Magnetic and Electric Fields - %s' % xc.title)
    #set color of axis spines and ticklabels
    _color_twin_axes(ax_B, _B_color, ax_E, _E_color)
    #legend
    ax_B.legend(kw['H'], kw['L'], numpoints=1)
    _format_line_axes_legends(ax_B, ax_E)
    _format_twin_axes(ax_B, ax_E)
    #save the fig or don't, depending on keywords
    _save_fig(xc.sheet, fig, **kw)
    #return
    return(fig, ax_B, ax_E)

def plot_xc(xc, **kw):
    """Plot CrossSection conductors and ROW edges without any fields
    args:
        xc - CrossSection objects
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
    _plot_wires(ax, xc.hot, xc.gnd, [c.y for c in xc.conds], **kw)
    _check_und_conds([xc], [ax], **kw)
    #adjust axis limits if headspace is called for
    if(_include_headspace):
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
    #plot ROW lines
    _plot_ROW_edges(ax, xc.lROW, xc.rROW, **kw)
    #set axis text and legend
    ax.set_title('Cross Section Configuration - %s' % xc.title)
    ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
    ax.set_ylabel(r'Height Above Ground $(ft)$')
    ax.legend(kw['H'], kw['L'], numpoints=1)
    #format
    _format_line_axes_legends(ax)
    #save or don't
    _save_fig(xc.sheet, fig, **kw)
    #return
    return(fig, ax)

def _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, field, **kw):
    """Handle plotting of DAT comparison features that don't require unique
    strings
    args:
        ax_abs - axis of absolute error plot
        ax_per - axis of percentage error plot
        ax_mag - axis of field magnitude plot
        pan - pandas.Panel object containing results and errors
        field - string, column label of field to be plotted (Bmax/Emax)"""

    #plot absolute error
    h_abs = ax_abs.plot(pan['Absolute Difference'][field], 'k', zorder=-1)
    ax_abs.set_ylabel(r'Absolute Difference $(kV/m)$')
    #plot percentage error
    h_per = ax_per.plot(pan['Percent Difference'][field], 'r', zorder=-1)
    ax_per.set_ylabel('Percent Difference', color = 'r')
    #set error axes legend
    ax_abs.legend(h_abs + h_per, ['Absolute Difference','Percent Difference'])
    #plot results
    h_fld = ax_mag.plot(pan['FIELDS_DAT_results'][field], 'k', zorder=-1)
    h_nm = ax_mag.plot(pan['python_results'][field], 'b', zorder=-1)
    ax_mag.set_xlabel(r'Distance from ROW Center $(ft)$')
    #set results legend
    ax_mag.legend(h_fld + h_nm + kw['H'], ['FIELDS','New Code'] + kw['L'],
            numpoints = 1)

def _plot_DAT_comparison(xc, pan, **kw):
    """Generate 2 subplots showing the FIELDS results (from a .DAT file)
    compared to the results of this code and the error.
    args:
        xc - CrossSection object
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
    _plot_wires(ax_mag, xc.hot, xc.gnd, pan['python_results']['Bmax'], **kw)
    _check_und_conds([xc], [ax_mag], **kw)
    _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Bmax', **kw)
    ax_abs.set_title('Absolute and Percent Difference, Max Magnetic Field')
    ax_mag.set_ylabel(r'Bmax $(mG)$')
    ax_mag.set_title('Model Results, Magnetic Field')
    _color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'r')
    _format_line_axes_legends(ax_abs)
    plt.tight_layout()
    _format_line_axes_legends(ax_abs, ax_per, ax_mag)
    _format_twin_axes(ax_abs, ax_per, ax_mag)
    _save_fig(xc.sheet + '-DAT_comparison_Bmax', fig, **kw)
    plt.close(fig)

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Emax
    #init handles and labels lists for legend
    kw['H'], kw['L'] = [], []
    _plot_wires(ax_mag, xc.hot, xc.gnd, pan['python_results']['Emax'], **kw)
    _check_und_conds([xc], [ax_mag], **kw)
    _plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Emax', **kw)
    ax_abs.set_title('Absolute and Percent Difference, Max Electric Field')
    ax_mag.set_ylabel(r'Emax $(kV/m)$')
    ax_mag.set_title('Model Results, Electric Field')
    _color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'r')
    plt.tight_layout()
    _format_line_axes_legends(ax_abs, ax_per, ax_mag)
    _format_twin_axes(ax_abs, ax_per, ax_mag)
    _save_fig(xc.sheet + '-DAT_comparison_Emax', fig, **kw)
    plt.close(fig)

#-------------------------------------------------------------------------------
#plotting routines working primarily with a SectionBook object

#useful globals for the section book plotting routine(s), unlikely to collide
#with other variables of the same name
_fields_plots_sb_headspace = 0.4 #space at the top of plots for legend
_fields_plots_sb_wireperc = 0.3 #percent of max field value to scale wire heights

def _plot_group_fields(ax, xcs, field, **kw):
    """Plot the results of fields calculations
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from
        field - column label of the results to plot
    kw:
        xmax - cutoff distance from ROW center"""
    #check for an xmax keyword
    if('xmax' in kw):
        xmax = kw['xmax']
    else:
        xmax = False
    #plot the fields, keeping handles and finding max of all
    fields_list = [xc.fields[field] for xc in xcs]
    for i in range(len(fields_list)):
        #plot
        if(xmax):
            kw['H'].append(ax.plot(fields_list[i][-xmax:xmax],
                    color = _colormap[i%7],
                    linewidth = _fields_linewidth)[0])
        else:
            kw['H'].append(ax.plot(fields_list[i], color = _colormap[i%7],
                    linewidth = _fields_linewidth)[0])
        kw['L'].append(xcs[i].sheet)

def _plot_group_wires(ax, xcs, max_field, **kw):
    """Plot the conductors of 1 or 2 CrossSections, using split color
    for conductor locations shared by 2 CrossSections
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from
        max_field - maximum value of all fields plotted in ax"""
    if(len(xcs) == 2):
        #---only hot conductors are plotted
        #---use sets to group shared and unshared conductor locations
        #get x and y pairs of Conductors in each group
        xy_0 = [(c.x,c.y) for c in xcs[0].hot]
        xy_1 = [(c.x,c.y) for c in xcs[1].hot]
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
        scale = _fields_plots_xc_wireperc*max_field/np.max(np.absolute(all_y))
        #cross section 0 conductors only
        if(len(xy_0) > 0):
            x,y = zip(*xy_0)
            kw['H'].append(ax.plot(x, scale*np.array(y), 'd',
                            color = _colormap[0])[0])
        else:
            #still need a handle for the legend
            kw['H'].append(mpl.lines.Line2D([], [], marker = 'd',
                    linestyle = '', color = _colormap[0]))
        kw['L'].append('Conductors - ' + xcs[0].sheet)
        #cross section 1 conductors only
        if(len(xy_1) > 0):
            x,y = zip(*xy_1)
            kw['H'].append(ax.plot(x, scale*np.array(y), 'd',
                            color = _colormap[1])[0])
        else:
            #still need a handle for the legend
            kw['H'].append(mpl.lines.Line2D([], [], marker = 'd',
                    linestyle = '', color = _colormap[1]))
        kw['L'].append('Conductors - ' + xcs[1].sheet)
        #shared conductors
        if(len(shared) > 0):
            x,y = zip(*shared)
            ax.plot(x, scale*np.array(y), 'd', color = _colormap[0],
                        fillstyle = 'left')
            ax.plot(x, scale*np.array(y), 'd', color = _colormap[1],
                        fillstyle = 'right')
        #underground line adjustments
        _check_und_conds(xcs, [ax], **kw)

def _plot_group_ROW_edges(ax, xcs, **kw):
    """Plot lines delineating the Right-of-Way boundaries
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from"""

    #if there are only two CrossSections, handle ROW edges independently
    if(len(xcs) == 2):
        yl = ax.get_ylim()
        #check if the left ROW edges are in the same place for both sections
        if(xcs[0].lROW == xcs[1].lROW):
            #plot overlapping dashed lines to create double colored dashes
            ax.plot([xcs[0].lROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1,
                        dashes = [6,6,6,6])
            ax.plot([xcs[1].lROW]*2, yl, '--', color = _colormap[1],
                        linewidth = _ROW_linewidth, zorder = -1,
                        dashes = [6,18,6,18])
        else:
            ax.plot([xcs[0].lROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1)
            ax.plot([xcs[1].lROW]*2, yl, '--', color = _colormap[1],
                        linewidth = _ROW_linewidth, zorder = -1)
        #check if the left ROW edges are in the same place for both sections
        if(xcs[0].rROW == xcs[1].rROW):
            #plot overlapping dashed lines to create double colored dashes
            ax.plot([xcs[0].rROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1,
                        dashes = [6,6,6,6])
            ax.plot([xcs[1].rROW]*2, yl, '--', color = _colormap[1],
                        linewidth = _ROW_linewidth, zorder = -1,
                        dashes = [6,18,6,18])
        else:
            ax.plot([xcs[0].rROW]*2, yl, '--', color = _colormap[0],
                        linewidth = _ROW_linewidth, zorder = -1)
            ax.plot([xcs[1].rROW]*2, yl, '--', color = _colormap[1],
                        linewidth = _ROW_linewidth, zorder = -1)
        #create a line handle and label for each color of the dashed line
        kw['H'] += [
                mpl.lines.Line2D([], [], linestyle = '--',color = _colormap[0]),
                mpl.lines.Line2D([], [], linestyle = '--', color = _colormap[1])
                ]
        kw['L'] += ['ROW Edges - ' + xcs[0].sheet,
                        'ROW Edges - ' + xcs[1].sheet]
    elif(len(xcs) > 2):
        #if there are more than two CrossSections and they all have the same
        #ROW edges, plot the single set of edge lines
        l, r = [xc.lROW for xc in xcs], [xc.rROW for xc in xcs]
        if(all([i == l[0] for i in l[1:]]) and all([i == r[0] for i in r[1:]])):
            _plot_ROW_edges(ax, l[0], r[0], **kw)

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
        return_figs - toggle whether a list of figure objects is returned
                      instead of closing the figures to clear memory,
                      default is False
    returns:
        figs - dict of dicts, keys are 'E' and 'B', which are each keyed
                by group tags, leading to Figure objects"""
    #check return kwarg
    return_figs = False
    if('return_figs' in kw):
        if(kw['return_figs']):
            return_figs = True
            figs = {'E': {}, 'B': {}}

    #iterate over groups with more than 1 CrossSection
    for xcs in sb.tag_groups:

        #BMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #init handles and labels lists for legend
        kw['H'], kw['L'] = [], []
        #plot the Bmax results for each xc in the group
        _plot_group_fields(ax, xcs, 'Bmax', **kw)
        #plot wires
        max_field = max([xc.fields['Bmax'].max() for xc in xcs])
        _plot_group_wires(ax, xcs, max_field, **kw)
        #adjust axis limits if called for
        if(_include_headspace):
            yl = ax.get_ylim()
            ax.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
        #plot ROW lines
        _plot_group_ROW_edges(ax, xcs, **kw)
        #set axis text and legend
        ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
        ax.set_ylabel(r'Maximum Magnetic Field $(mG)$')
        t = 'Maximum Magnetic Field - %s' % str(xcs[0].tag)
        ax.set_title(t)
        ax.legend(kw['H'], kw['L'], numpoints = 1)
        _format_line_axes_legends(ax)
        #save the figure if keyword 'save' == True, and append fig to figs
        _save_fig('group_%s-Bmax' % str(xcs[0].tag), fig, **kw)
        #store the fig or close it
        if(return_figs):
            figs['B'][xcs[0].tag] = fig
        else:
            plt.close(fig)

        #EMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #init handles and labels lists for legend
        kw['H'], kw['L'] = [], []
        #plot the Bmax results for each xc in the group
        _plot_group_fields(ax, xcs, 'Emax', **kw)
        #plot wires
        max_field = max([xc.fields['Emax'].max() for xc in xcs])
        _plot_group_wires(ax, xcs, max_field, **kw)
        #adjust axis limits if called for
        if(_include_headspace):
            yl = ax.get_ylim()
            ax.set_ylim(yl[0], (1 + _fields_plots_xc_headspace)*yl[1])
        #plot ROW lines
        _plot_group_ROW_edges(ax, xcs, **kw)
        #set axis text and legendf
        ax.set_xlabel(r'Distance from Center of ROW $(ft)$')
        ax.set_ylabel(r'Maximum Electric Field $(kV/m)$')
        t = 'Maximum Electric Field - %s' % str(xcs[0].tag)
        ax.set_title(t)
        ax.legend(kw['H'], kw['L'], numpoints = 1)
        _format_line_axes_legends(ax)

        #save the figure if keyword 'save' == True or a path string is passed
        _save_fig('group_%s-Emax' % str(xcs[0].tag), fig, **kw)
        #store the fig or close it
        if(return_figs):
            figs['E'][xcs[0].tag] = fig
        else:
            plt.close(fig)

    if(return_figs):
        return(figs)

def _reorder_xcs(xcs, **kw):
    """Use the xc_order kw, if it exists, to create a list of
    CrossSection sheet strings that determines the order bars are plotted
    in
    args:
        xcs - list of CrossSections corresponding to a group
    kw:
        xc_order - dict, if any keys are the same as the tag of the
                   CrossSection group represented by xcs, they should
                   map to a list of strings with CrossSection sheet names
                   specifiying an order to plot in.
    returns:
        reorder - a reordered list of CrossSections"""

    #see if there is an xc_order kwarg
    if('xc_order' in kw):
        xc_order = kw['xc_order']
        tag = xcs[0].tag
        #see if it corresponds to the group represented by xcs
        if(tag in xc_order):
            #get the order
            order = xc_order[tag]
            #compile xc sheet names
            sheets = [xc.sheet for xc in xcs]
            #make a dict for later
            d = dict(zip(sheets, xcs))
            #reorder the sheets
            reorder = []
            for sh in order:
                if(sh in sheets):
                    reorder.append(sheets.pop(sheets.index(sh)))
            #put leftovers at the end of the list
            reorder = reorder + sheets
            #generate reodered list and return
            return([d[k] for k in reorder])

    #return the order of xcs
    return(xcs)

def _plot_group_bars(ax, xcs, field, side):
    """Plot bars in ax for each xc in xcs
    args:
        ax - target Axes object
        xcs - list of CrossSection objects corresponding to a group
        field - string, field component to plot (Ex, Bx, Bprod, etc.)
        side - string, 'left' or 'right', selects which ROW edge side"""

    #translate side to index
    if(side == 'left'):
        side = 0
    else:
        side = 1
    #plot the bars
    x = range(len(xcs))
    values = [xc.ROW_edge_fields[field].values[side] for xc in xcs]
    ax.bar(x, values, color=[_colormap[i%7] for i in range(len(xcs))],
            bottom=0.0, align='center', alpha=0.8, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([textwrap.fill(xc.sheet, 15) for xc in xcs],
            rotation='vertical', fontsize=10)

def _generate_ROW_value_plot_objects(xcs):
    """generate figure and axes for plot_groups_at_ROW
    args:
        xcs - list of CrossSections in group
    returns:
        fig - figure object
        axl - left axes
        axr - right axes"""

    #see if ROW edge distances are unique
    lROW, rROW = set([xc.lROW for xc in xcs]), set([xc.rROW for xc in xcs])
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
    kw"
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
        xc_order - dict, keys can be CrossSection group tags, which map to
                   lists of strings specifying the order of the plotted
                   CrossSection bars (left to right). Not all CrossSecitons
                   in a group must be listed. Any/all can be left out.
        return_figs - toggle whether a list of figure objects is returned
                      instead of closing the figures to clear memory,
                      default is False
    returns:
        figs - dict of dicts, keys are 'E' and 'B', which are each keyed
                by group tags, leading to Figure objects"""
    #check return kwarg
    return_figs = False
    if('return_figs' in kw):
        if(kw['return_figs']):
            return_figs = True
            figs = {'E': {}, 'B': {}}

    #iterate over groups with more than 1 CrossSection
    for xcs in sb.tag_groups:

        #get reordered CrossSection list
        xcs = _reorder_xcs(xcs, **kw)

        #plot Bmax
        fig, axl, axr = _generate_ROW_value_plot_objects(xcs)
        _plot_group_bars(axl, xcs, 'Bmax', 'left')
        _plot_group_bars(axr, xcs, 'Bmax', 'right')
        #format
        _format_bar_axes_legends(axl, axr)
        #apply text
        axl.set_ylabel(r'Maximum Magnetic Field $(mG)$')
        fig.suptitle('Maximum Magnetic Fields at ROW Edges - %s' % xcs[0].tag,
                fontsize=18)
        #save?
        _save_fig('group_%s-ROW-Bmax' % str(xcs[0].tag), fig, **kw)
        #store the fig or close it
        if(return_figs):
            figs['B'][xcs[0].tag] = fig
        else:
            plt.close(fig)

        #plot Emax
        fig, axl, axr = _generate_ROW_value_plot_objects(xcs)
        _plot_group_bars(axl, xcs, 'Emax', 'left')
        _plot_group_bars(axr, xcs, 'Emax', 'right')
        #format
        _format_bar_axes_legends(axl, axr)
        #apply text
        axl.set_ylabel(r'Maximum Electric Field $(kV/m)$')
        fig.suptitle('Maximum Electric Fields at ROW Edges - %s' % xcs[0].tag,
                fontsize=18)
        #save?
        _save_fig('group_%s-ROW-Emax' % str(xcs[0].tag), fig, **kw)
        #store the fig or close it
        if(return_figs):
            figs['E'][xcs[0].tag] = fig
        else:
            plt.close(fig)

    if(return_figs):
        return(figs)
