import numpy as np
import matplotlib.pyplot as plt

import emf_class
import emf_funks
import emf_calcs

#useful globals used by any/all emf_plots functions
emf_plots_figsize = (10,6) #figure inches (w,h), use None to revert to default

#-------------------------------------------------------------------------------
#plotting routines working with a CrossSection object

#useful globals for the cross section plotting routines, unlikely to collide
#with other variables of the same name
emf_plots_xc_headspace = 0.4 #space at the top of plots for legend
emf_plots_xc_wireperc = 0.3 #percent of max field value to scale wire heights

def prepare_fig(xc, **kwargs):
    """Snippet executed at the beginning of plotting methods to handle figure
    object generation and some keywords initializing other params."""
    #prepare figure and axis
    plt.rc('font', family = 'calibri')
    k = kwargs
    keys = k.keys()
    if('figure' in keys):
        fig = k['figure']
    else:
        fig = plt.figure(figsize = emf_plots_figsize)
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

def save_fig(xc, fig, **kwargs):
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
            fn = emf_funks.path_manage(xc.name, '', **kwargs)
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
    hROW = ax.plot([lROW]*2, yl, 'k--', [rROW]*2, yl, 'k--')
    xl = ax.get_xlim()
    if((xl[0] == lROW) or (xl[1] == rROW)):
        ax.set_xlim((xl[0]*1.15, xl[1]*1.15))
    return(hROW)

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
    hB, = ax.plot(xc.fields['Bmax'][-xmax:xmax], linesym, color = xc.B_color)
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
                'ROW Edge'], numpoints = 1, fontsize = 12)
    #save the fig or don't, depending on keywords
    save_fig(xc, fig, **kwargs)
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
    hE, = ax.plot(xc.fields['Emax'][-xmax:xmax], linesym, color = xc.E_color)
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
                'ROW Edge'], numpoints = 1, fontsize = 12)
    #save the fig or don't, depending on keywords
    save_fig(xc, fig, **kwargs)
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
    hB, = ax_B.plot(xc.fields['Bmax'][-xmax:xmax], linesym, color = xc.B_color)
    hE, = ax_E.plot(xc.fields['Emax'][-xmax:xmax], linesym, color = xc.E_color)
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
                    fontsize = 14, color = xc.B_color)
    ax_E.set_ylabel('Maximum Electric Field (kV/m)',
                    fontsize = 14, color = xc.E_color)
    if('title' in kwargs.keys()):
        t = kwargs['title']
    else:
        t = '%s, Maximum Magnetic and Electric Fields' % xc.title
    ax_B.set_title(t)
    #set color of axis spines and ticklabels
    color_twin_axes(ax_B, xc.B_color, ax_E, xc.E_color)
    #legend
    ax_B.legend([hB, hE, hhot, hgnd, hROW[0]],
                ['Magnetic Field (mG)','Electric Field (kV/m)','Conductors',
                'Grounded Conductors','ROW Edge'],
                numpoints = 1, fontsize = 12)
    #save the fig or don't, depending on keywords
    save_fig(xc, fig, **kwargs)
    #return
    return(fig)

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
            h.append(ax.plot(fields_list[i][-xmax:xmax])[0])
            l.append(xcs[i].title)
        else:
            h.append(ax.plot(fields_list[i])[0])
            l.append(xcs[i].title)
        #find max
        if(max(fields_list[i]) > max_field):
            max_field = max(fields_list[i])
    return(h, l, max_field)

def plot_group_wires(ax, xcs, max_field, **kwargs):
    #get x and y pairs of Conductors in each group
    pairs = []
    for xc in xcs:
        for c in xc.hot:
            pairs.append((c.x, c.y))
    #separate shared and unique Conductors
    shared = []
    unique = []
    for i in range(len(pairs)):
        xy = pairs[i]
        j = 0
        found = False
        while((j < len(pairs)) and (not found)):
            if(j != i):
                if(xc == pairs[j]):
                    shared.append(xy)
                    found = True
            j += 1
        if(not found):
            unique.append(xy)
    #plot shared conductors
    h = []
    l = []
    if(shared):
        x,y = zip(*shared)
        x,y = np.array(x),np.array(y)
        scale = emf_plots_sb_wireperc*max_field/np.max(np.absolute(y))
        h.append(ax.plot(x, scale*y, 'd')[0])
        l.append('Shared Conductors')
    #plot unique conductors
    if(unique):
        x,y = zip(*unique)
        x,y = np.array(x),np.array(y)
        scale = emf_plots_sb_wireperc*max_field/np.max(np.absolute(y))
        h.append(ax.plot(x, scale*y, 'kd')[0])
        l.append('Unique Conductors')
    #return handles
    return(h, l)

def plot_group_ROW_edges(ax, xcs):
    yl = ax.get_ylim()
    left = [xc.lROW for xc in xcs]
    right = [xc.rROW for xc in xcs]
    #if all ROW edges are the same, just plot one set
    if((len(np.unique(left)) == 1) and (len(np.unique(right)) == 1)):
        l = [np.unique(left)]*2
        r = [np.unique(right)]*2
        h = ax.plot(l, yl, 'k--', r, yl, 'k--')
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
        ax = fig.add_subplot(1,1,1)
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
        ax.legend(h+hw+hROW, l+lw+lROW, numpoints = 1, fontsize = 12)
        figs.append(fig)
        #save fig

        #EMAX

    return(figs)
