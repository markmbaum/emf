import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
lines = mpl.lines

import emf_funks

#rcparams for more static global formatting changes
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.figsize'] = (12, 6)
mpl.rcParams['font.family'] = 'calibri'
mpl.rcParams['text.color'] = (.2, .2, .2)
mpl.rcParams['axes.labelcolor'] = (.2, .2, .2)
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.borderaxespad'] = 0 #default is None
mpl.rcParams['xtick.color'] = (.2, .2, .2)
mpl.rcParams['ytick.color'] = (.2, .2, .2)

#other more specific/dynamic global formatting variables
emf_plots_B_color = 'darkgreen'
emf_plots_E_color = 'midnightblue'
emf_plots_fields_linewidth = 2
emf_plots_ROW_linewidth = 0.75
emf_plots_ROW_color = 'gray'
emf_plots_ax_frameon = False
emf_plots_ax_ticks_on = False
emf_plots_leg_edge_on = False
emf_plots_colormap = [(0, 0.4470, 0.7410),(0.8500, 0.3250,0.0980),
                    (0.9290, 0.6940, 0.1250),(0.4940, 0.1840, 0.5560),
                    (0.4660, 0.6740, 0.1880),(0.3010, 0.7450, 0.9330),
                    (0.6350, 0.0780, 0.1840)]

#-------------------------------------------------------------------------------
#general plotting support functions

def ion():
    plt.ion()

def show():
    plt.show()

def format_axes_legends(*args):
    """Apply axis formatting commands to axes objects
    args: some number of axis objects"""
    for ax in args:
        #apply axis formatting
        ax.set_frame_on(emf_plots_ax_frameon)
        if(not emf_plots_ax_ticks_on):
            ax.tick_params(axis = 'both', which = 'both',
                    bottom = 'off', top = 'off', left = 'off', right = 'off')
        #apply legend formatting
        leg = ax.get_legend()
        if(leg):
            rec = leg.get_frame()
            if(not emf_plots_leg_edge_on):
                rec.set_edgecolor('white')

def save_fig(filename_if_needed, fig, **kwargs):
    """Snippet executed at the end of plotting methods to handle saving
    args:
        filename_if_needed - string used for filename if it's not in 'path'
        fig - figure to save
    kwargs:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""
    #force saving if a path is passed in
    if('path' in kwargs):
        kwargs['save'] = True
    #condition filename and format strings for saving
    if('save' in kwargs):
        if(kwargs['save']):
            #get filename
            fn = emf_funks.path_manage(filename_if_needed, '', **kwargs)
            #get format/extension
            if('format' in kwargs):
                fmt = kwargs['format']
                if('.' in fmt):
                    fmt = fmt[fmt.index('.')+1:]
            else:
                fmt = 'png'
            #save the plot
            fn += '.' + fmt
            plt.savefig(fn, format = fmt)
            print('plot saved to: "%s"' % fn)

def color_twin_axes(ax1, color1, ax2, color2):
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

#useful globals for the CrossSection plotting routines, unlikely to collide
#with other variables of the same name
emf_plots_xc_headspace = 0.4 #space at the top of plots for legend
emf_plots_xc_wireperc = 0.3 #percent of max field value to scale wire heights

def prepare_fig(xc, **kwargs):
    """Snippet executed at the beginning of plotting methods to handle figure
    object generation and some keywords initializing other params.
    args:
        xc - CrossSection objects
    kwargs:
        figure - figure object to recycle, if needed
        xmax - cutoff distance from the ROW center"""
    #prepare figure and axis
    if('figure' in kwargs):
        fig = kwargs['figure']
    else:
        fig = plt.figure()
    ax = plt.gca()
    #get x cutoff, if any
    if('xmax' in kwargs):
        xmax = abs(kwargs['xmax'])
    else:
        xmax = max(abs(xc.fields.index))
    #get appropriate linestyle based on number of sample points
    n = xc.fields[-xmax:xmax].shape[0]
    if(n > 100):
        linesym = '-'
    else:
        linesym = '.-'
    return(fig, ax, xmax, linesym)

def plot_wires(ax, hot, gnd, v):
    """Plot conductor symbols in ax.Returns handles for the hot and gnd
    Conductors.
    args:
        ax - target axis
        hot - list of non-grounded conductors
        gnd - list of grounded conductors
        v - iterable of calculated field results used to scale conductors"""
    #get x and y coordinates
    x = np.array([c.x for c in hot + gnd])
    y = np.array([c.y for c in hot + gnd])
    #calculate the scaling factor
    scale = emf_plots_xc_wireperc*np.max(v)/np.max(np.absolute(y))
    #bring underground lines to zero
    y[y < 0.] = 0.
    #plot
    h, l = [], []
    if(hot):
        h.append(ax.plot(x[:len(hot)], scale*y[:len(hot)], 'kd')[0])
        l.append('Conductors')
    if(gnd):
        h.append(ax.plot(x[len(hot):], scale*y[len(hot):], 'd', color = 'gray')[0])
        l.append('Grounded Conductors')
    return(h, l)

def plot_ROW_edges(ax, lROW, rROW):
    """Plot dashed lines marking the left and right edges of the
    Right-of-Way. Axis limits are also adjusted to allow extra space on the
    sides to make the ROW edge lines visible if needed. The ROW edge line
    handles are returned in a list.
    args:
        ax - target axis
        lROW - x value of the left edge of the ROW
        rROW - x value of the right edge of the ROW"""
    yl = ax.get_ylim()
    hROW = ax.plot([lROW]*2, yl, '--', color = emf_plots_ROW_color,
                    linewidth = emf_plots_ROW_linewidth, zorder = -1)
    ax.plot([rROW]*2, yl, '--', color = emf_plots_ROW_color,
            linewidth = emf_plots_ROW_linewidth, zorder = -1)
    l = ['ROW Edges']
    xl = ax.get_xlim()
    if((xl[0] == lROW) or (xl[1] == rROW)):
        ax.set_xlim((xl[0]*1.15, xl[1]*1.15))
    return(hROW, l)

def plot_Bmax(xc, **kwargs):
    """Plot the maximum magnetic field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
    kwargs:
        figure - matplotlib Figure object, target figure for plotting
        title - string, exact plot title
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""
    #get axes and x cutoff
    (fig, ax, xmax, linesym) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hB = ax.plot(xc.fields['Bmax'][-xmax:xmax], linesym,
            color = emf_plots_B_color, linewidth = emf_plots_fields_linewidth)
    lB = ['Magnetic Field (mG)']
    #plot wires
    hw, lw = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Bmax'])
    #adjust axis limits
    ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Bmax']))
    #plot ROW lines
    hROW, lROW = plot_ROW_edges(ax, xc.lROW, xc.rROW)
    #set axis text and legend
    ax.set_xlabel('Distance from Center of ROW (ft)')
    ax.set_ylabel('Maximum Magnetic Field (mG)')
    if('title' in kwargs):
        t = kwargs['title']
    else:
        t = 'Maximum Magnetic Field - %s' % xc.title
    ax.set_title(t)
    ax.legend(hB + hw + hROW, lB + lw + lROW, numpoints = 1)
    format_axes_legends(ax)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_Emax(xc, **kwargs):
    """Plot the maximum electric field along the ROW with conductor
    locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
    kwargs:
        figure - matplotlib Figure object, target figure for plotting
        title - string, exact plot title
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""
    #get axes and x cutoff
    (fig, ax, xmax, linesym) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hE = ax.plot(xc.fields['Emax'][-xmax:xmax], linesym,
            color = emf_plots_E_color, linewidth = emf_plots_fields_linewidth)
    lE = ['Electric Field (kV/m)']
    #plot wires
    hw, lw = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Emax'])
    #adjust axis limits
    ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Emax']))
    #plot ROW lines
    hROW, lROW = plot_ROW_edges(ax, xc.lROW, xc.rROW)
    #set axis text and legend
    ax.set_xlabel('Distance from Center of ROW (ft)')
    ax.set_ylabel('Maximum Electric Field (kV/m)')
    if('title' in kwargs):
        t = kwargs['title']
    else:
        t = 'Maximum Electric Field - %s' % xc.title
    ax.set_title(t)
    ax.legend(hE + hw + hROW, lE + lw + lROW, numpoints = 1)
    format_axes_legends(ax)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_max_fields(xc, **kwargs):
    """Plot the maximum magnetic and electric field on split vertical axes,
    along the ROW with conductor locations shown in artificial coordinates.
    args:
        xc - a CrossSection object
    kwargs:
        figure - matplotlib Figure object, target figure for plotting
        title - string, exact plot title
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""
    #get axes and x cutoff
    (fig, ax_B, xmax, linesym) = prepare_fig(xc, **kwargs)
    ax_E = ax_B.twinx()
    #plot the field curves
    hf = [ax_B.plot(xc.fields['Bmax'][-xmax:xmax], linesym,
        color = emf_plots_B_color, linewidth = emf_plots_fields_linewidth)[0],
        ax_E.plot(xc.fields['Emax'][-xmax:xmax], linesym,
        color = emf_plots_E_color, linewidth = emf_plots_fields_linewidth)[0]]
    lf = ['Magnetic Field (mG)', 'Electric Field (kV/m)']
    #plot wires
    hw, lw = plot_wires(ax_B, xc.hot, xc.gnd, xc.fields['Bmax'])
    #adjust axis limits
    ax_B.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Bmax']))
    ax_E.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Emax']))
    #plot ROW lines
    hROW, lROW = plot_ROW_edges(ax_B, xc.lROW, xc.rROW)
    #set axis text
    ax_B.set_xlabel('Distance from Center of ROW (ft)')
    ax_B.set_ylabel('Maximum Magnetic Field (mG)', color = emf_plots_B_color)
    ax_E.set_ylabel('Maximum Electric Field (kV/m)', color = emf_plots_E_color)
    if('title' in kwargs):
        t = kwargs['title']
    else:
        t = 'Maximum Magnetic and Electric Fields - %s' % xc.title
    ax_B.set_title(t)
    #set color of axis spines and ticklabels
    color_twin_axes(ax_B, emf_plots_B_color, ax_E, emf_plots_E_color)
    #legend
    ax_B.legend(hf + hw + hROW, lf + lw + lROW, numpoints = 1)
    format_axes_legends(ax_B, ax_E)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, field, hw, lw):
    """Handle plotting of DAT comparison features that don't require unique
    strings
    args:
        ax_abs - axis of absolute error plot
        ax_per - axis of percentage error plot
        ax_mag - axis of field magnitude plot
        pan - pandas.Panel object containing results and errors
        field - string, column label of field to be plotted (Bmax/Emax)
        hw - handles of the conductor symbol plots
        lw - labels for conductor symbols"""
    #plot absolute error
    h_abs = ax_abs.plot(pan['Absolute Difference'][field], 'k')
    ax_abs.set_ylabel('Absolute Difference (kV/m)')
    #plot percentage error
    h_per = ax_per.plot(pan['Percent Difference'][field], 'r')
    ax_per.set_ylabel('Percent Difference', color = 'r')
    #set error axes legend
    ax_abs.legend(h_abs + h_per, ['Absolute Difference','Percent Difference'])
    #plot results
    h_fld = ax_mag.plot(pan['FIELDS_output'][field], 'k')
    h_nm = ax_mag.plot(pan['New_model_output'][field], 'b')
    ax_mag.set_xlabel('Distance from ROW Center (ft)')
    #set results legend
    ax_mag.legend(h_fld + h_nm + hw, ['FIELDS','New Code'] + lw, numpoints = 1)

def plot_DAT_comparison(xc, pan, **kwargs):
    """Generate 2 subplots showing the FIELDS results (from a .DAT file)
    compared to the results of this code and the error.
    args:
        xc - CrossSection object
        pan - pandas.Panel object containing results and errors
    kwargs:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Bmax

    hw, lw = plot_wires(ax_mag, xc.hot, xc.gnd, pan['New_model_output']['Bmax'])
    plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Bmax', hw, lw)
    ax_abs.set_title('Absolute and Percent Difference, Max Magnetic Field')
    ax_mag.set_ylabel('Bmax (mG)')
    ax_mag.set_title('Model Results, Magnetic Field')
    color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'r')
    format_axes_legends(ax_abs)
    plt.tight_layout()
    format_axes_legends(ax_abs, ax_per, ax_mag)
    save_fig(xc.name + '-DAT_comparison_Bmax', fig, **kwargs)
    plt.close(fig)

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Emax
    hw, lw = plot_wires(ax_mag, xc.hot, xc.gnd, pan['New_model_output']['Emax'])
    plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Emax', hw, lw)
    ax_abs.set_title('Absolute and Percent Difference, Max Electric Field')
    ax_mag.set_ylabel('Emax (kV/m)')
    ax_mag.set_title('Model Results, Electric Field')
    color_twin_axes(ax_abs, mpl.rcParams['axes.labelcolor'], ax_per, 'r')
    plt.tight_layout()
    format_axes_legends(ax_abs, ax_per, ax_mag)
    save_fig(xc.name + '-DAT_comparison_Emax', fig, **kwargs)
    plt.close(fig)


#-------------------------------------------------------------------------------
#plotting routines working primarily with a SectionBook object

#useful globals for the section book plotting routine(s), unlikely to collide
#with other variables of the same name
emf_plots_sb_headspace = 0.4 #space at the top of plots for legend
emf_plots_sb_wireperc = 0.3 #percent of max field value to scale wire heights

def plot_group_fields(ax, xcs, field, **kwargs):
    """Plot the results of fields calculations
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from
        field - column label of the results to plot
    kwargs:
        xmax - cutoff distance from ROW center"""
    #check for an xmax keyword
    if('xmax' in kwargs):
        xmax = kwargs['xmax']
    else:
        xmax = False
    #plot the fields, keeping handles and finding max of all
    fields_list = [xc.fields[field] for xc in xcs]
    h = []
    l = []
    max_field = 0.
    for i in range(len(fields_list)):
        #plot
        if(xmax):
            h.append(ax.plot(fields_list[i][-xmax:xmax],
                    color = emf_plots_colormap[i%7],
                    linewidth = emf_plots_fields_linewidth)[0])
        else:
            h.append(ax.plot(fields_list[i], color = emf_plots_colormap[i%7],
                    linewidth = emf_plots_fields_linewidth)[0])
        l.append(field + ' - ' + xcs[i].title)
        #find max
        if(max(fields_list[i]) > max_field):
            max_field = max(fields_list[i])
    return(h, l, max_field)

def plot_group_wires(ax, xcs, max_field):
    """Plot the conductors of 1 or 2 CrossSections, using split color
    for conductor locations shared by 2 CrossSections
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from
        max_field - maximum value of all fields plotted in ax"""
    #conductor markers not available for more than 2 cross sections for now
    h, l = [], []
    if(len(xcs) == 2):
        #---only hot conductors are plotted
        #---use sets to group shared and unshared conductor locations
        #get x and y pairs of Conductors in each group
        xy_0 = [(c.x,c.y) for c in xcs[0].hot]
        xy_1 = [(c.x,c.y) for c in xcs[1].hot]
        #zero the y coordinate of any underground lines
        for i in range(len(xy_0)):
            if(xy_0[i][1] < 0):
                xy_0[i] = (xy_0[i][0], 0.0)
        for i in range(len(xy_1)):
            if(xy_1[i][1] < 0):
                xy_1[i] = (xy_1[i][0], 0.0)
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
        #plot shared conductors
        h = []
        l = []
        scale = emf_plots_xc_wireperc*max_field/np.max(np.absolute(all_y))
        #cross section 0 conductors only
        if(len(xy_0) > 0):
            x,y = zip(*xy_0)
            h.append(ax.plot(x, scale*np.array(y), 'd',
                    color = emf_plots_colormap[0])[0])
        else:
            #still need a handle for the legend
            h.append(lines.Line2D([], [], marker = 'd', linestyle = '',
                    color = emf_plots_colormap[0]))
        l.append('Conductors - ' + xcs[0].title)
        #cross section 1 conductors only
        if(len(xy_1) > 0):
            x,y = zip(*xy_1)
            h.append(ax.plot(x, scale*np.array(y), 'd',
                    color = emf_plots_colormap[1])[0])
        else:
            #still need a handle for the legend
            h.append(lines.Line2D([], [], marker = 'd', linestyle = '',
                    color = emf_plots_colormap[1]))
        l.append('Conductors - ' + xcs[1].title)
        #shared conductors
        if(len(shared) > 0):
            x,y = zip(*shared)
            ax.plot(x, scale*np.array(y), 'd', color = emf_plots_colormap[0],
                    fillstyle = 'left')
            ax.plot(x, scale*np.array(y), 'd', color = emf_plots_colormap[1],
                    fillstyle = 'right')
    #return handles and labels
    return(h, l)

def plot_group_ROW_edges(ax, xcs):
    """Plot lines delineating the Right-of-Way boundaries
    args:
        ax - target axis
        xcs - list of CrossSection objects to plot results from"""
    h, l = [], []
    if(len(xcs) == 2):
        yl = ax.get_ylim()
        #check if the left ROW edges are in the same place for both sections
        if(xcs[0].lROW == xcs[1].lROW):
            #plot overlapping dashed lines to create double colored dashes
            ax.plot([xcs[0].lROW]*2, yl, '--', color = emf_plots_colormap[0],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1,
                        dashes = [6,6,6,6])
            ax.plot([xcs[1].lROW]*2, yl, '--', color = emf_plots_colormap[1],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1,
                        dashes = [6,18,6,18])
        else:
            ax.plot([xcs[0].lROW]*2, yl, '--', color = emf_plots_colormap[0],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1)
            ax.plot([xcs[1].lROW]*2, yl, '--', color = emf_plots_colormap[1],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1)
        #check if the left ROW edges are in the same place for both sections
        if(xcs[0].rROW == xcs[1].rROW):
            #plot overlapping dashed lines to create double colored dashes
            ax.plot([xcs[0].rROW]*2, yl, '--', color = emf_plots_colormap[0],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1,
                        dashes = [6,6,6,6])
            ax.plot([xcs[1].rROW]*2, yl, '--', color = emf_plots_colormap[1],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1,
                        dashes = [6,18,6,18])
        else:
            ax.plot([xcs[0].rROW]*2, yl, '--', color = emf_plots_colormap[0],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1)
            ax.plot([xcs[1].rROW]*2, yl, '--', color = emf_plots_colormap[1],
                        linewidth = emf_plots_ROW_linewidth, zorder = -1)
        #create a line handle and label for each color of the dashed line
        h = [lines.Line2D([], [], linestyle = '--',color = emf_plots_colormap[0]),
            lines.Line2D([], [], linestyle = '--', color = emf_plots_colormap[1])]
        l = ['ROW Edges - ' + xcs[0].title, 'ROW Edges - ' + xcs[1].title]

    return(h, l)

def plot_groups(sb, **kwargs):
    """Plot the fields of grouped CrossSections in the same axis, a plot for
    both fields. Only plots groups of more than one CrossSection.
    args:
        sb - SectionBook object to pull plotting groups from
    kwargs:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')
        xmax - cutoff distance from ROW center"""
    #iterate over groups with more than 1 CrossSection
    for g in [group for group in sb.tag_groups if len(group) > 1]:
        xcs = [sb.xcs[i] for i in g]
        #BMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, frameon = emf_plots_ax_frameon)
        #plot the Bmax results for each xc in the group
        h, l, max_field = plot_group_fields(ax, xcs, 'Bmax', **kwargs)
        #plot wires
        hw, lw = plot_group_wires(ax, xcs, max_field)
        #adjust axis limits
        ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max_field)
        #plot ROW lines
        hROW, lROW = plot_group_ROW_edges(ax, xcs)
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)')
        ax.set_ylabel('Maximum Magnetic Field (mG)')
        t = 'Maximum Magnetic Field - '
        for i in range(len(xcs)):
            t += xcs[i].title + ', '
        t = t[:-2]
        ax.set_title(t)
        ax.legend(h+hw+hROW, l+lw+lROW, numpoints = 1)
        format_axes_legends(ax)
        #save the figure if keyword 'save' == True, and append fig to figs
        save_fig('group_%s-Bmax' % str(xcs[0].tag), fig, **kwargs)
        plt.close(fig)

        #EMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, frameon = emf_plots_ax_frameon)
        #plot the Bmax results for each xc in the group
        h, l, max_field = plot_group_fields(ax, xcs, 'Emax', **kwargs)
        #plot wires
        hw, lw = plot_group_wires(ax, xcs, max_field)
        #adjust axis limits
        ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max_field)
        #plot ROW lines
        hROW, lROW = plot_group_ROW_edges(ax, xcs)
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)')
        ax.set_ylabel('Maximum Electric Field (kV/m)')
        t = 'Maximum Electric Field - '
        for i in range(len(xcs)):
            t += xcs[i].title + ', '
        t = t[:-2]
        ax.set_title(t)
        ax.legend(h+hw+hROW, l+lw+lROW, numpoints = 1)
        format_axes_legends(ax)
        #save the figure if keyword 'save' == True, and append fig to figs
        save_fig('group_%s-Emax' % str(xcs[0].tag), fig, **kwargs)
        plt.close(fig)
