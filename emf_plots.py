import numpy as np
import matplotlib as mpl
plt = mpl.pyplot
lines = mpl.lines

import emf_funks

#useful globals used by any/all emf_plots functions
#rcparams for static formatting changes
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['font.family'] = 'calibri'
mpl.rcParams['text.color'] = (.2, .2, .2)
mpl.rcParams['axes.labelcolor'] = (.2, .2, .2)
mpl.rcParams['xtick.color'] = (.2, .2, .2)
mpl.rcParams['ytick.color'] = (.2, .2, .2)
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['figure.figsize'] = (10,6)
#other more specific and dynamic formatting variables
emf_plots_B_color = 'darkgreen'
emf_plots_E_color = 'midnightblue'
emf_plots_ROW_color = 'gray'
emf_plots_ax_frameon = False
emf_plots_ax_ticks_on = False
emf_plots_leg_edge_on = False
emf_plots_colormap = [(0, 0.4470, 0.7410),(0.8500, 0.3250,0.0980),
                    (0.9290, 0.6940, 0.1250),(0.4940, 0.1840, 0.5560),
                    (0.4660, 0.6740, 0.1880),(0.3010, 0.7450, 0.9330),
                    (0.6350, 0.0780, 0.1840)]
#-------------------------------------------------------------------------------
#general plotting routines

def format_axes_legends(*args):
    """Apply axis formatting commands to axes objects"""
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
                print('why cant I get the legend exactly in the corner?')

def save_fig(filename_if_needed, fig, **kwargs):
    """Snippet executed at the end of plotting methods to handle saving"""
    #force saving if a path is passed in
    if('path' in kwargs.keys()):
        kwargs['save'] = True
    #save the fig, or don't
    k = kwargs
    keys = k.keys()
    #condition filename and format strings for saving
    if('save' in keys):
        if(k['save']):
            #get filename
            fn = emf_funks.path_manage(filename_if_needed, '', **kwargs)
            #get format/extension
            if('format' in keys):
                fmt = k['format']
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
    object generation and some keywords initializing other params."""
    #prepare figure and axis
    k = kwargs
    keys = k.keys()
    if('figure' in keys):
        fig = k['figure']
    else:
        fig = plt.figure()
    ax = plt.gca()
    #get x cutoff, if any
    if('xmax' in keys):
        xmax = abs(k['xmax'])
    else:
        xmax = max(abs(xc.fields.index))
    #get appropriate linestyle based on number of sample points
    n = xc.fields[-xmax:xmax].shape[0]
    if(n > 200):
        linesym = '-'
    else:
        linesym = '.-'
    return(fig, ax, xmax, linesym)

def plot_wires(ax, hot, gnd, v):
    """Plot conductor symbols in ax, where hot and gnd are lists of
    Conductor objects, v is an iterable used to scale the Conductor heights,
    and perc is the percentage of the max value in v that the Conductor
    heights are scaled to. Returns handles for the hot and gnd conductors."""
    #get x and y coordinates
    x = np.array([c.x for c in hot + gnd])
    y = np.array([c.y for c in hot + gnd])
    #calculate the scaling factor
    scale = emf_plots_xc_wireperc*np.max(v)/np.max(np.absolute(y))
    #bring underground lines to zero
    y[y < 0.] = 0.
    #plot
    hhot, = ax.plot(x[:len(hot)], scale*y[:len(hot)], 'kd')
    hgnd, = ax.plot(x[len(hot):], scale*y[len(hot):], 'd', color = 'gray')
    return(hhot, hgnd)

def plot_ROW_edges(ax, lROW, rROW):
    """Plot dashed lines marking the left and right edges of the
    Right-of-Way in ax, the locations of which are given by lROW and rROW.
    The iterable v is used to scale the lines. Axis limits are also adjusted
    to allow extra space on the sides to make the ROW edge lines visible if needed.
    The ROW edge line handles are returned in a list."""
    yl = ax.get_ylim()
    hROW = ax.plot([lROW]*2, yl, '--', color = emf_plots_ROW_color)
    ax.plot([rROW]*2, yl, '--', color = emf_plots_ROW_color)
    xl = ax.get_xlim()
    if((xl[0] == lROW) or (xl[1] == rROW)):
        ax.set_xlim((xl[0]*1.15, xl[1]*1.15))
    return(hROW)

def plot_Bmax(xc, **kwargs):
    """Plot the maximum magnetic field along the ROW with conductor
    locations shown in artificial but to-scale locations. Pass in an
    existing figure with keyword argument 'figure' to recycle an object.
    Pass in a plot title with the keyword argument 'title' to specify an
    exact title, otherwise the title will be used. Use the kwarg
    'xmax' to cut the plotted fields at a certain distance from the ROW
    center. If the keyword argument 'save' is passed in, the plot
    will be saved. Use the keyword argument 'path' to specify the path
    or filename of the saved plot. The format can also be specified
    (usually 'png' or 'pdf') with the 'format' keyword."""
    #get axes and x cutoff
    (fig, ax, xmax, linesym) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hB, = ax.plot(xc.fields['Bmax'][-xmax:xmax], linesym, color = emf_plots_B_color)
    #plot wires
    hhot, hgnd = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Bmax'])
    #adjust axis limits
    ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Bmax']))
    #plot ROW lines
    hROW = plot_ROW_and_adjust(ax, xc.lROW, xc.rROW)
    #set axis text and legend
    ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
    ax.set_ylabel('Maximum Magnetic Field (mG)', fontsize = 14)
    if('title' in keys):
        t = k['title']
    else:
        t = 'Maximum Magnetic Field, %s' % xc.title
    ax.set_title(t)
    ax.legend(['Magnetic Field (mG)','Conductors','Grounded Conductors',
                'ROW Edge'], numpoints = 1, fontsize = 10)
    format_axes_legends(ax)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_Emax(xc, **kwargs):
    """Plot the maximum electric field along the ROW with conductor
    locations shown in artificial but to-scale locations. Pass in an
    existing figure with keyword argument 'figure' to recycle an object.
    Pass in a plot title with the keyword argument 'title' to specify an
    exact title, otherwise the title will be used. Use the kwarg
    'xmax' to cut the plotted fields at a certain distance from the ROW
    center. If the keyword argument 'save' is passed in True, the plot
    will be saved. Use the keyword argument 'path' to specify the path
    or filename of the saved plot. The format can also be specified
    (usually 'png' or 'pdf') with the 'format' keyword."""
    #get axes and x cutoff
    (fig, ax, xmax, linesym) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hE, = ax.plot(xc.fields['Emax'][-xmax:xmax], linesym, color = emf_plots_E_color)
    #plot wires
    hhot, hgnd = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Emax'])
    #adjust axis limits
    ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Emax']))
    #plot ROW lines
    hROW = plot_ROW_edges(ax, xc.lROW, xc.rROW)
    #set axis text and legend
    ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
    ax.set_ylabel('Maximum Electric Field (kV/m)', fontsize = 14)
    if('title' in keys):
        t = k['title']
    else:
        t = 'Maximum Electric Field, %s' % xc.title
    ax.set_title(t)
    ax.legend(['Electric Field (kV/m)','Conductors','Grounded Conductors',
                'ROW Edge'], numpoints = 1, fontsize = 10)
    format_axes_legends(ax)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_max_fields(xc, **kwargs):
    """Plot the maximum fields along the ROW with conductor
    locations shown in artificial but to-scale locations. Pass in an
    existing figure with keyword argument 'figure' to recycle an object.
    Pass in a plot title with the keyword argument 'title' to specify an
    exact title, otherwise the title will be used. Use the kwarg
    'xmax' to cut the plotted fields at a certain distance from the ROW
    center. If the keyword argument 'save' is passed in, the plot
    will be saved. Use the keyword argument 'path' to specify the path
    or filename of the saved plot. The format can also be specified
    (usually 'png' or 'pdf') with the 'format' keyword."""
    #get axes and x cutoff
    (fig, ax_B, xmax, linesym) = prepare_fig(xc, **kwargs)
    ax_E = ax_B.twinx()
    #plot the field curves
    hB, = ax_B.plot(xc.fields['Bmax'][-xmax:xmax], linesym, color = emf_plots_B_color)
    hE, = ax_E.plot(xc.fields['Emax'][-xmax:xmax], linesym, color = emf_plots_E_color)
    #plot wires
    hhot, hgnd = plot_wires(ax_B, xc.hot, xc.gnd, xc.fields['Bmax'])
    #adjust axis limits
    ax_B.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Bmax']))
    ax_E.set_ylim(0, (1 + emf_plots_xc_headspace)*max(xc.fields['Emax']))
    #plot ROW lines
    hROW = plot_ROW_edges(ax_B, xc.lROW, xc.rROW)
    #set axis text
    ax_B.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
    ax_B.set_ylabel('Maximum Magnetic Field (mG)',
                    fontsize = 14, color = emf_plots_B_color)
    ax_E.set_ylabel('Maximum Electric Field (kV/m)',
                    fontsize = 14, color = emf_plots_E_color)
    if('title' in kwargs.keys()):
        t = kwargs['title']
    else:
        t = '%s, Maximum Magnetic and Electric Fields' % xc.title
    ax_B.set_title(t)
    #set color of axis spines and ticklabels
    color_twin_axes(ax_B, emf_plots_B_color, ax_E, emf_plots_E_color)
    #legend
    ax_B.legend([hB, hE, hhot, hgnd, hROW[0]],
                ['Magnetic Field (mG)','Electric Field (kV/m)','Conductors',
                'Grounded Conductors','ROW Edge'],
                numpoints = 1, fontsize = 10)
    format_axes_legends(ax_B, ax_E)
    #save the fig or don't, depending on keywords
    save_fig(xc.name, fig, **kwargs)
    #return
    return(fig)

def plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, field):
    #plot absolute error
    h_abs, = ax_abs.plot(pan['Absolute Difference'][field], 'k')
    ax_abs.set_ylabel('Absolute Difference (kV/m)')
    #plot percentage error
    h_per, = ax_per.plot(pan['Percent Difference'][field], 'r')
    ax_per.set_ylabel('Percent Difference', color = 'r')
    #set error axes legend
    ax_abs.legend([h_abs,h_per], ['Absolute Difference','Percent Difference'],
                    fontsize = 10)
    #plot results
    h_fld, = ax_mag.plot(pan['FIELDS_output'][field], 'k')
    h_nm, = ax_mag.plot(pan['New_model_output'][field], 'b')
    ax_mag.set_xlabel('Distance from ROW Center (ft)')
    #set results legend
    ax_mag.legend([h_fld, h_nm], ['FIELDS','New Code'], fontsize = 10)

def plot_DAT_comparison(xc, pan, **kwargs):
    figs = []
    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Bmax
    plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Bmax')
    ax_abs.set_title('Absolute and Percent Difference, Max Magnetic Field')
    ax_mag.set_ylabel('Bmax (mG)')
    ax_mag.set_title('Model Results, Magnetic Field')
    color_twin_axes(ax_abs, 'k', ax_per, 'r')
    format_axes_legends(ax_abs)
    plt.tight_layout()
    format_axes_legends(ax_abs, ax_per, ax_mag)
    save_fig(xc.name + '-DAT_comparison_Emax', fig, **kwargs)
    figs.append(fig)

    #figure object and axes
    fig = plt.figure()
    ax_abs = fig.add_subplot(2,1,1)
    ax_per = ax_abs.twinx()
    ax_mag = fig.add_subplot(2,1,2)
    #Emax
    plot_DAT_repeatables(ax_abs, ax_per, ax_mag, pan, 'Emax')
    ax_abs.set_title('Absolute and Percent Difference, Max Electric Field')
    ax_mag.set_ylabel('Emax (kV/m)')
    ax_mag.set_title('Model Results, Electric Field')
    color_twin_axes(ax_abs, 'k', ax_per, 'r')
    plt.tight_layout()
    format_axes_legends(ax_abs, ax_per, ax_mag)
    save_fig(xc.name + '-DAT_comparison_Bmax', fig, **kwargs)
    figs.append(fig)

    return(figs)

#-------------------------------------------------------------------------------
#plotting routines working primarily with a SectionBook object

#useful globals for the section book plotting routine(s), unlikely to collide
#with other variables of the same name
emf_plots_sb_headspace = 0.4 #space at the top of plots for legend
emf_plots_sb_wireperc = 0.3 #percent of max field value to scale wire heights

def plot_group_fields(ax, xcs, field, **kwargs):
    #check for an xmax keyword
    if('xmax' in kwargs.keys()):
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
                    color = emf_plots_colormap[i%7])[0])
        else:
            h.append(ax.plot(fields_list[i], color = emf_plots_colormap[i%7])[0])
        l.append(xcs[i].title + ' ' + field)
        #find max
        if(max(fields_list[i]) > max_field):
            max_field = max(fields_list[i])
    return(h, l, max_field)

def plot_group_wires(ax, xcs, max_field, **kwargs):
    #conductor markers not available for more than 2 cross sections, the split
    #color symbols are too complex and the plot gets cluttered
    h, l = [], []
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
        l.append(xcs[0].title + ' Conductors')
        #cross section 1 conductors only
        if(len(xy_1) > 0):
            x,y = zip(*xy_1)
            h.append(ax.plot(x, scale*np.array(y), 'd',
                    color = emf_plots_colormap[1])[0])
        else:
            #still need a handle for the legend
            h.append(lines.Line2D([], [], marker = 'd', linestyle = '',
                    color = emf_plots_colormap[1]))
        l.append(xcs[1].title + ' Conductors')
        #shared conductors
        if(len(shared) > 0):
            x,y = zip(*shared)
            ax.plot(x, scale*np.array(y), 'd', color = emf_plots_colormap[0],
                    fillstyle = 'left')
            ax.plot(x, scale*np.array(y), 'd', color = emf_plots_colormap[1],
                    fillstyle = 'right')
    elif(len(xcs) == 1):
        #single set of conductors
        xc = xcs[0]
        x = np.array([c.x for c in xc.hot])
        y = np.array([c.y for c in xc.hot])
        scale = emf_plots_xc_wireperc*max_field/max(abs(y))
        h = ax.plot(x, scale*y, 'd', color = emf_plots_colormap[0])
        l = [xc.title + ' Conductors']
    #return handles and labels
    return(h, l)

def plot_group_ROW_edges(ax, xcs):
    yl = ax.get_ylim()
    left = [xc.lROW for xc in xcs]
    right = [xc.rROW for xc in xcs]
    #if all ROW edges are the same, just plot one set
    if((len(np.unique(left)) == 1) and (len(np.unique(right)) == 1)):
        l = [np.unique(left)]*2
        r = [np.unique(right)]*2
        h = ax.plot(l, yl, '--', color = emf_plots_ROW_color)
        ax.plot(r, yl, '--', color = emf_plots_ROW_color)
        l = ['ROW Edges']
    else:
        pass
        #unfinished

    return(h, l)

def plot_groups(sb, **kwargs):
    """Plot the fields of grouped CrossSections in the same axis, a plot for
    both fields. Use the 'save' keyword to save each plot and the 'path'
    keyword to specify the save location (don't inculde a file name or the
    plots will overwrite eachother). Use the kwarg 'xmax' to cut the plotted
    fields off at a certain distance from the ROW center. A list of figure
    objects is returned."""
    #figure list
    figs = []
    #iterate over the groups
    for g in sb.tag_groups:
        xcs = [sb.xcs[i] for i in g]

        #BMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, frameon = emf_plots_ax_frameon)
        #plot the Bmax results for each xc in the group
        h, l, max_field = plot_group_fields(ax, xcs, 'Bmax', **kwargs)
        #plot wires
        hw, lw = plot_group_wires(ax, xcs, max_field, **kwargs)
        #adjust axis limits
        ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max_field)
        #plot ROW lines
        hROW, lROW = plot_group_ROW_edges(ax, xcs)
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax.set_ylabel('Maximum Magnetic Field (mG)', fontsize = 14)
        t = 'Maximum Magnetic Field - '
        for i in range(len(xcs)):
            t += xcs[i].title + ', '
        t = t[:-2]
        ax.set_title(t)
        ax.legend(h+hw+hROW, l+lw+lROW, numpoints = 1, fontsize = 10)
        format_axes_legends(ax)
        #save the figure if keyword 'save' == True, and append fig to figs
        save_fig('B_field-group_%s' % str(xcs[0].tag), fig, **kwargs)
        figs.append(fig)

        #EMAX
        #get plotting objects
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, frameon = emf_plots_ax_frameon)
        #plot the Bmax results for each xc in the group
        h, l, max_field = plot_group_fields(ax, xcs, 'Emax', **kwargs)
        #plot wires
        hw, lw = plot_group_wires(ax, xcs, max_field, **kwargs)
        #adjust axis limits
        ax.set_ylim(0, (1 + emf_plots_xc_headspace)*max_field)
        #plot ROW lines
        hROW, lROW = plot_group_ROW_edges(ax, xcs)
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax.set_ylabel('Maximum Electric Field (kV/m)', fontsize = 14)
        t = 'Maximum Electric Field - '
        for i in range(len(xcs)):
            t += xcs[i].title + ', '
        t = t[:-2]
        ax.set_title(t)
        ax.legend(h+hw+hROW, l+lw+lROW, numpoints = 1, fontsize = 10)
        format_axes_legends(ax)
        #save the figure if keyword 'save' == True, and append fig to figs
        save_fig('E_field-group_%s' % str(xcs[0].tag), fig, **kwargs)
        figs.append(fig)

    return(figs)
