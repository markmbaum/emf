from . import plt

import emf_funks

def _save_fig(filename_if_needed, fig, **kwargs):
    """Snippet executed at the end of plotting methods to handle saving
    args:
        filename_if_needed - string used for filename if it's not in 'path'
        fig - figure to save
    kwargs:
        save - bool, toggle plot saving
        path - string, destination/filename for saved figure
        format - string, saved plot format/extension (default 'png')"""

    #set current figure
    plt.figure(fig.number)

    #force saving if a path is passed in
    if('path' in kwargs):
        kwargs['save'] = True
    #condition filename and format strings for saving
    if('save' in kwargs):
        if(kwargs['save']):
            #get filename
            fn = emf_funks._path_manage(filename_if_needed, '', **kwargs)
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
            print('plot saved to: %s' % fn)
