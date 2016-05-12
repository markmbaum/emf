import numpy as np
import matplotlib.pyplot as plt

import emf_class
import emf_funks
import emf_calcs

#-------------------------------------------------------------------------------
#plotting routines working primarily with a CrossSection object

def prepare_fig(xc, **kwargs):
    """Snippet executed at the beginning of plotting methods to handle figure
    object generation and some keywords"""
    #prepare figure and axis
    plt.rc('font', family = 'calibri')
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
    return(fig, ax, xmax)

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

def plot_wires(ax, hot, gnd, v, perc):
    """Plot conductor symbols in ax, where hot and gnd are lists of
    Conductor objects, v is an iterable used to scale the Conductor heights,
    and perc is the percentage of the max value in v that the Conductor
    heights are scaled to. Returns handles for the hot and gnd conductors."""
    #get x and y coordinates
    x = np.array([c.x for c in hot + gnd])
    y = np.array([c.y for c in hot + gnd])
    #calculate the scaling factor
    scale = perc*np.max(v)/np.max(np.absolute(y))
    #bring underground lines to zero
    y[y < 0.] = 0.
    #plot
    hhot, = ax.plot(x[:len(hot)], scale*y[:len(hot)], 'kd')
    hgnd, = ax.plot(x[len(hot):], scale*y[len(hot):], 'd', color = 'gray')
    return(hhot, hgnd)

def plot_ROW_and_adjust(ax, lROW, rROW, v, headroom):
    """Plot dashed lines marking the left and right edges of the
    Right-of-Way in ax, the locations of which are given by lROW and rROW.
    The iterable v is used to scale the lines. Axis limits are also adjusted
    to allow a percentage of empty space at the top, defined by headroom, and
    extra space on the sides to make the ROW edge lines visible if needed.
    The ROW edge line handles are returned in a list."""
    ax.set_ylim([0, max(v)*(1 + headroom)])
    yl = ax.get_ylim()
    hROW = ax.plot([lROW]*2, yl, 'k--', [rROW]*2, yl, 'k--')
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
    #keyword variables
    k = kwargs
    keys = k.keys()
    #get axes and x cutoff
    (fig, ax, xmax) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hB, = ax.plot(xc.fields['Bmax'][-xmax:xmax], '.-', color = xc.B_color)
    #plot wires
    hhot, hgnd = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Bmax'], .3)
    #plot ROW lines and adjust axis limits for legend and ROW lines
    hROW = plot_ROW_and_adjust(ax, xc.lROW, xc.rROW, xc.fields['Bmax'], .35)
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
    #keyword variables
    k = kwargs
    keys = k.keys()
    #get axes and x cutoff
    (fig, ax, xmax) = prepare_fig(xc, **kwargs)
    #plot the field curve
    hE, = ax.plot(xc.fields['Emax'][-xmax:xmax], '.-', color = xc.E_color)
    #plot wires
    hhot, hgnd = plot_wires(ax, xc.hot, xc.gnd, xc.fields['Emax'], .3)
    #plot ROW lines and adjust axis limits for legend and ROW lines
    hROW = plot_ROW_and_adjust(ax, xc.lROW, xc.rROW, xc.fields['Emax'], .35)
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
    #keyword variables
    k = kwargs
    keys = k.keys()
    #get axes and x cutoff
    (fig, ax_B, xmax) = prepare_fig(xc, **kwargs)
    ax_E = ax_B.twinx()
    #plot the field curves
    hB, = ax_B.plot(xc.fields['Bmax'][-xmax:xmax], '.-', color = xc.B_color)
    hE, = ax_E.plot(xc.fields['Emax'][-xmax:xmax], '.-', color = xc.E_color)
    #plot wires
    hhot, hgnd = plot_wires(ax_B, xc.hot, xc.gnd, xc.fields['Bmax'], .3)
    #plot ROW lines and adjust axis limits for legend and ROW lines
    hROW = plot_ROW_and_adjust(ax_B, xc.lROW, xc.rROW, xc.fields['Emax'], .35)
    #set axis text
    ax_B.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
    ax_B.set_ylabel('Maximum Magnetic Field (mG)',
                    fontsize = 14, color = xc.B_color)
    ax_E.set_ylabel('Maximum Electric Field (kV/m)',
                    fontsize = 14, color = xc.E_color)
    if('title' in keys):
        t = k['title']
    else:
        t = '%s, Maximum Magnetic and Electric Fields' % xc.title
    ax_B.set_title(t)
    #set color of axis spines and ticklabels
    ax_B.spines['left'].set_color(xc.B_color)
    ax_B.spines['right'].set_color(xc.E_color)
    ax_E.spines['left'].set_color(xc.B_color)
    ax_E.spines['right'].set_color(xc.E_color)
    ax_B.tick_params(axis = 'y', colors = xc.B_color)
    ax_E.tick_params(axis = 'y', colors = xc.E_color)
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

def plot_groups(sb, **kwargs):
    """Plot the fields of CrossSections in the same group in the same
    axis, one plot for both fields. Use the 'save' keyword to save
    each plot and the 'path' keyword to specify the save location (and
    filename if it's in the path string). Use the kwarg 'xmax' to cut the
    plotted fields off at a certain distance from the ROW center."""
    figs = []
    if('xmax' in kwargs.keys()):
        xmax = kwargs['xmax']
    else:
        xmax = False
    #iterate over the groups
    for g in sb.tag_groups:
        #BMAX
        #get plotting objects and x cutoff
        fig, ax = plt.subplots()
        #iterate over the indices in the group, plotting Bmax results
        h = []
        l = []
        maxBmax = 0.
        for idx in g:
            xc = sb.xcs[idx]
            #B field
            if(xmax):
                han, = ax.plot(xc.fields['Bmax'][-xmax:xmax], '.-')
            else:
                han, = ax.plot(xc.fields['Bmax'], '.-')
            if(max(xc.fields['Bmax']) > maxBmax):
                maxBmax = max(xc.fields['Bmax'])
            h.append(han) #store handles for legend
            l.append(xc.title) #store titles for legend labels
        #plot ROW lines and adjust axis limits for legend and ROW lines
        ax.set_ylim([0, maxBmax*1.35])
        yl = ax.get_ylim()
        hROW = ax.plot([xc.lROW]*2, yl, 'k--', [xc.rROW]*2, yl, 'k--')
        xl = ax.get_xlim()
        if((xl[0] == xc.lROW) or (xl[1] == xc.rROW)):
            ax.set_xlim((xl[0]*1.15, xl[1]*1.15))
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax.set_ylabel('Maximum Magnetic Field (mG)', fontsize = 14)
        t = 'Maximum Magnetic Field - '
        for i in range(len(l)):
            t += l[i] + ', '
        t = t[:-2]
        ax.set_title(t)
        ax.legend(h + [hROW[0]], l + ['ROW Edge'], numpoints = 1, fontsize = 12)
        save_fig(xc, fig, **kwargs)
        figs.append(fig)

        #EMAX
        #

    return(figs)
