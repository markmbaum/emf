from . import plt

from . import emf_funks

def _save_fig(filename_if_needed, fig, **kw):
    """Snippet executed at the end of plotting methods to handle saving
    args:
        filename_if_needed - string used for filename if it's not in 'path'
        fig - figure to save
    kw:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""

    #set current figure
    plt.figure(fig.number)

    #force saving if a path is passed in
    if('path' in kw):
        kw['save'] = True
    #condition filename and format strings for saving
    if('save' in kw):
        if(kw['save']):
            #get filename
            fn = emf_funks._path_manage(filename_if_needed, '', **kw)
            #get format/extension
            if('format' in kw):
                fmt = kw['format']
                if('.' in fmt):
                    fmt = fmt[fmt.index('.')+1:]
            else:
                fmt = 'png'
            #save the plot
            fn += '.' + fmt
            plt.savefig(fn, format=fmt)
            print('plot saved to: %s' % fn)

def _prepare_fig(**kw):
    """Snippet executed at the beginning of plotting methods to handle figure
    object generation
    kw:
        ax - Axes object to plot in
        fig - Figure object to plot in (first Axes are targeted, fig.axes[0])"""
    #use a preexisting figure or create a new one
    if('ax' in kw):
        ax = kw['ax']
        fig = ax.figure
    elif('fig' in kw):
        fig = kw['fig']
        if(fig.axes):
            ax = fig.axes[0]
        else:
            ax = fig.add_subplot(1,1,1)
    else:
        fig, ax = plt.subplots(1,1)

    return(fig, ax)
