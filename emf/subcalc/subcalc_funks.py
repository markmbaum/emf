from .. import os, np, pd, shutil

from ..emf_funks import (_path_manage, _check_extension, _is_number, _is_int,
                        _check_intable, _flatten, _sig_figs, _check_to_array,
                        _Levenshtein_group)

import subcalc_class
import SUBCALC_io

def drop_footprint_template(*args, **kw):
    """Copy the emf.subcalc footprint template in the current directory or a directory specified by an input string
    args:
        drop_path - string, path of copied template file"""
    #check inputs
    if(len(args) > 1):
        raise(fields_class.EMFError("""drop_footprint_template only accepts zero or one input argument. A string can be passed to specify the directory in which the template file is copied. With no arguments, the template file is copied into the current directory."""))
    elif(len(args) == 1):
        kw = {'path': args[0]}
    #get template file path
    fn_temp = os.path.dirname(os.path.dirname(__file__))
    fn_temp = os.path.join(fn_temp, 'templates')
    fn_temp = os.path.join(fn_temp, 'subcalc-footprint-template.xlsx')
    #get drop path
    fn_drop = _path_manage('subcalc-footprint-template', 'xlsx', **kw)
    #check for existing files
    if(os.path.isfile(fn_drop)):
        raise(subcalc_class.EMFError('A file with the path "%s" already exists. Move/delete it or pass a new path string to drop_footprint_template().' % fn_drop))
    #copy and notify
    shutil.copyfile(fn_temp, fn_drop)
    print('subcalc footprint template written to: %s' % fn_drop)

def drop_tower_template(*args, **kw):
    """Copy the emf.subcalc tower template in the current directory or a directory specified by an input string
    args:
        drop_path - string, path of copied template file"""
    #check inputs
    if(len(args) > 1):
        raise(fields_class.EMFError("""drop_tower_template only accepts zero or one input argument. A string can be passed to specify the directory in which the template file is copied. With no arguments, the template file is copied into the current directory."""))
    elif(len(args) == 1):
        kw = {'path': args[0]}
    #get template file path
    fn_temp = os.path.dirname(os.path.dirname(__file__))
    fn_temp = os.path.join(fn_temp, 'templates')
    fn_temp = os.path.join(fn_temp, 'subcalc-tower-template.xlsx')
    #get drop path
    fn_drop = _path_manage('subcalc-tower-template', 'xlsx', **kw)
    #check for existing files
    if(os.path.isfile(fn_drop)):
        raise(subcalc_class.EMFError('A file with the path "%s" already exists. Move/delete it or pass a new path string to drop_tower_template().' % fn_drop))
    #copy and notify
    shutil.copyfile(fn_temp, fn_drop)
    print('subcalc tower template written to: %s' % fn_drop)

def load_towers(fn, return_model=False, **kw):
    """Read an emf.subcalc tower template or a SUBCALC INP file and create a new list of Tower objects. The list of Tower objects is returned by default, but a new Model object containing the Tower objects can be returned by using return_model=True.
    args:
        fn - string, the path string to the emf.subcalc tower template file
             (excel or csv) or a SUBCALC INP file
    optional args:
        return_model - bool, if False, a list of Tower objects is returned.
                       If True, a new Model object containing the Towers is
                       returned.
    kw:
        sheet - string, if the template file is an excel workbook with multiple sheets,
                the target sheet must be specified
    returns:
        a list of Towers or a Model object, depending on return_model"""

    #check if the target file is an INP
    if('.INP' in fn):
        towers = SUBCALC_io.read_INP(fn)
        name = os.path.basename(fn)
        if('.' in name):
            name = name[:name.rfind('.')]
    else:
        #load the template file into a DataFrame
        if('.' in fn):
            ext = fn[fn.rfind('.')+1:]
            if(ext == 'xlsx'):
                dfs = pd.read_excel(fn, sheetname=None)
                if(len(dfs.keys()) > 1):
                    if('sheet' in kw):
                        name = kw['sheet']
                        df = dfs[name]
                    else:
                        raise(subcalc_class.EMFError("""If an excel file with multiple sheets is passed to load_towers, the target sheet must be specified with the keyword argument 'sheet'."""))
                else:
                    df = dfs[dfs.keys()[0]]
            elif(ext == 'csv'):
                df = pd.read_csv(fn)
                name = os.path.basename(fn)
                if('.' in name):
                    name = name[:name.rfind('.')]
            else:
                raise(subcalc_class.EMFError("Only csv and xlsx files can be passed to load_towers."))
        else:
            raise(subcalc_class.EMFError("""No extension was detected at the end of file name "%s" passed to load_towers. File names passed to load_towers must have .csv or .xlsx extensions.""" % fn))
        #condition data a little
        df = df.dropna(how='all').fillna(method='ffill')
        #match the columns with string distance method
        cols = ['group', 'sequence', 'tower x', 'tower y', 'rotation', 'h', 'v', 'I', 'phase']
        df.columns = _Levenshtein_group(df.columns, cols)
        #parse to Towers
        towers = []
        for (group, seq), df in df.groupby(['group','sequence']):
            row = df.iloc[0]
            t = subcalc_class.Tower(group, seq, row['tower x'], row['tower y'], row['rotation'],
                    df['h'].values, df['v'].values, df['I'].values, df['phase'].values)
            towers.append(t)

    #return Towers directly or in a Model
    if(return_model):
        return(subcalc_class.Model(name=name.replace('-towers',''), towers=towers))
    else:
        return(towers)

def load_results(*args, **kw):
    """Read a .REF output file and load the data into a Results object
    args:
        results_path - string, path to the output .REF file of field results or to
                the excel file exported by a Results object
        footprint_path - string, optional, path to the csv file of
                         footprint data
    kw:
        Bkey - string, sets 'component' of magnetic field results that the
               returned Results object accesses by default
                     - can be 'Bx', 'By', 'Bz', 'Bmax', or 'Bres'
                     - default is 'Bmax'
                     - all components are stored, none are lost
    returns
        res - Results object containing results"""

    #check for a Bkey kwarg
    if('Bkey' in kw):
        Bkey = kw['Bkey']
    else:
        Bkey = 'Bmax'

    #check extensions
    try:
        fn = _check_extension(args[0], '.REF', '')
    except(subcalc_class.EMFError):
        fn = _check_extension(args[0], '.xlsx', """
        Can only load Results from .REF or .xlsx files""")

    if(fn[-3:] == 'REF'):

        #pull data from the REF file
        data, info = SUBCALC_io.read_REF(args[0])
        #get the gridded arrays
        data = mesh_dict_grids(data)
        #initialize Results object
        res = subcalc_class.Results(data, info, Bkey=Bkey)
        #check for footprint file path and load if present
        if(len(args) > 1):
            res.load_footprints(args[1])

    elif(fn[-4:] == 'xlsx'):
        #get a dict of all sheets in excel file
        dfs = pd.read_excel(args[0], sheetname=None)
        bkeys = dfs.keys()
        if('info' in bkeys):
            bkeys.remove('info')
        if('footprints' in bkeys):
            bkeys.remove('footprints')
        #slice out grid data
        x = [float(i) for i in dfs[bkeys[0]].columns]
        y = [float(i) for i in dfs[bkeys[0]].index]
        X, Y = np.meshgrid(x, y)
        data = {'X': X, 'Y': Y}
        for k in bkeys:
            data[str(k)] = dfs[k].values
        #slice out info dictionary
        if('info' in dfs):
            info = dfs['info']
            params = info[info.columns[0]].values
            values = info[info.columns[1]].values
            info = dict(zip(params, values))
            #initialize Results with metadata
            res = subcalc_class.Results(data, info, Bkey=Bkey)
        else:
            #initialize Results object without metadata dict
            res = subcalc_class.Results(data, Bkey=Bkey)

        #check for footprints
        if(len(args) > 1):
            #check for footprint file path and load if present
            res.load_footprints(args[1])
        elif('footprints' in dfs):
            res.load_footprints(dfs['footprints'])

    else:
        raise(subcalc_class.EMFError("""
        Results must be loaded from .REF file or excel files"""))

    #set the Results object's name with the input filename
    n = os.path.basename(args[0])
    if('.' in n):
        n = n[:n.rfind('.')]
    res.name = n

    #return
    return(res)

def mesh_dict_grids(flat_data):
    """Convert raw grid data read from a SubCalc output file (by the read_REF function) into meshed grids of X, Y coordinates and their corresponding B field values. This function would generally only be used on the results of subcalc_funks.read_REF to reshape 1D arrays of values pulled from REF files. However, the load_results function calls these functions to read a REF file, reshape the results, and return them in a Results object. Thus, the load_results function is almost always the more convenient way to load REF data into Python objects.
    args:
        flat_data - dict, keyed by 'x', 'y', 'bx', 'by', 'bz', 'bmax', 'bres'
    returns:
        grid_data - dict with 2D arrays keyed by
                'X', 'Y', 'Bx', 'By', 'Bz', 'Bmax', 'Bres'
    """
    #find the number of points in a row
    x = flat_data['x']
    y = flat_data['y']
    count = 0
    v = y[0]
    while(y[count] == v):
        count += 1
    #get ncols and nrows
    L = len(x)
    ncols = count
    nrows = L/ncols
    #map old to new keys
    mapk = dict(zip(['x','y','bx','by','bz','bmax','bres'],
                    ['X','Y','Bx','By','Bz','Bmax','Bres']))
    #replace with 2D arrays
    k = [mapk[i] for i in flat_data]
    v = [np.reshape(flat_data[i], (nrows, ncols)) for i in flat_data]
    grid_data = dict(zip(k, v))

    return(grid_data)

def _bilinear_interp(res, x, y):
    """Use Results to interpolate linearly in two dimensions for an estimate of any x,y coordinate inside the grid.
    args:
        res - Results object
        x - float, x coordinate to interpolate at
        y - float, y coordinate to interpolate at
    returns:
        B_interp - float, interpolated field value"""
    #first find the 4 point grid cell containing x,y
    #   (the point is assumed to lie inside the grid)
    _, xidx = _double_min(np.abs(res.x - x))
    _, yidx = _double_min(np.abs(res.y - y))
    #get coordinates and values
    x1, x2 = res.x[xidx]
    y1, y2 = res.y[yidx]
    B11 = res.B[yidx[0], xidx[0]]
    B12 = res.B[yidx[0], xidx[1]]
    B21 = res.B[yidx[1], xidx[0]]
    B22 = res.B[yidx[1], xidx[1]]
    #interpolate
    ym1 = y - y1
    xm1 = x - x1
    y2m = y2 - y
    x2m = x2 - x
    B_interp = ((1.0/((x2 - x1)*(y2 - y1)))
                *(x2m*(B11*y2m + B12*ym1) + xm1*(B21*y2m + B22*ym1)))

    return(B_interp)

def _2Dmax(G):
    """Find the indices and value of the maximum value in a 2 dimensional array
    args:
        G - 2D numpy array
    returns:
        m - the maximum value
        i - index of max along 0th axis
        j - index of max along 1st axis"""
    imax, jmax = 0, 0
    m = np.min(G)
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            if(G[i,j] > m):
                m = G[i,j]
                imax = i
                jmax = j
    return(m, imax, jmax)

def _2Dmin(G):
    """Find the indices and value of the minimum value in a 2 dimensional array
    args:
        G - 2D numpy array
    returns:
        m - the minimum value
        i - index of max along 0th axis
        j - index of max along 1st axis"""
    imin, jmin = 0, 0
    m = np.max(G)
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            if(G[i,j] < m):
                m = G[i,j]
                imin = i
                jmin = j
    return(m, imin, jmin)

def _double_min(v):
    """Find the lowest two values in an array and their indices
    args:
        v - iterable
    returns:
        mins - array of minima, the first one being the smallest
        idxs - array of indices of minima"""
    if(len(v) < 2):
        raise(subcalc_class.EMFError("""
        Cannot find lowest two values in an array of length less than 2."""))
    m = max(v) #store the max for initialization
    mins = np.array([m, m], dtype=float)
    idxs = np.zeros((2,), dtype=int)
    for i in range(len(v)):
        if(v[i] < mins[0]):
            #swap first minimum to second
            mins[1] = mins[0]
            idxs[1] = idxs[0]
            #store new first minimum
            mins[0] = v[i]
            idxs[0] = i
        elif(v[i] < mins[1]):
            #store new second minimum
            mins[1] = v[i]
            idxs[1] = i

    return(mins, idxs)

_dist = lambda x1, x2, y1, y2: np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
_dist3 = lambda x1, x2, y1, y2, z1, z2: np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

def _interp_path_2D(x, y, n):
    """Interpolate points along linear segments of a path in 2d space so
    that the total number of points is approximately n by distributing the
    points along the segments (according to their lenghts). The end points
    of each segment along the path are preserved.
    args:
        x - iterable of x coordinates
        y - iterable of y coordinates
        n - target number of samples
    returns:
        x_interp - array of x coordinates
        y_interp - array of y coordinates"""
    L = len(x) - 1
    rL = range(L)
    #calculate distances of each segment
    d = np.array([_dist(x[i], x[i+1], y[i], y[i+1]) for i in rL])
    #convert distances to fractions
    d = d/sum(d)
    #approximately distribute the total number of points to each segment
    n = np.ceil(n*d)
    #make sure there are at least two samples in each segment
    n[n < 1] = 1
    #interpolate over each segment
    x_interp = _flatten([np.linspace(x[i], x[i+1], n[i], False) for i in rL])
    y_interp = _flatten([np.linspace(y[i], y[i+1], n[i], False) for i in rL])
    x_interp.append(float(x[-1]))
    y_interp.append(float(y[-1]))
    return(np.array(x_interp), np.array(y_interp))

def _interp_path_3D(x, y, z, n):
    """Interpolate points along linear segments of a path in 2d space so
    that the total number of points is approximately n by distributing the
    points along the segments (according to their lenghts). The end points
    of each segment along the path are preserved.
    args:
        x - iterable of x coordinates
        y - iterable of y coordinates
        z - iterable of z coordinates
        n - target number of samples
    returns:
        x_interp - array of x coordinates
        y_interp - array of y coordinates
        z_interp - array of z coordinates"""
    L = len(x) - 1
    rL = range(L)
    #calculate distances of each segment
    d = np.array([_dist3(x[i], x[i+1], y[i], y[i+1], z[i], z[i+1]) for i in rL])
    #convert distances to fractions
    d = d/sum(d)
    #approximately distribute the total number of points to each segment
    n = np.ceil(n*d)
    #make sure there are at least two samples in each segment
    n[n < 1] = 1
    #interpolate over each segment
    x_interp = _flatten([np.linspace(x[i], x[i+1], n[i], False) for i in rL])
    y_interp = _flatten([np.linspace(y[i], y[i+1], n[i], False) for i in rL])
    z_interp = _flatten([np.linspace(z[i], z[i+1], n[i], False) for i in rL])
    x_interp.append(float(x[-1]))
    y_interp.append(float(y[-1]))
    z_interp.append(float(z[-1]))
    return(np.array(x_interp), np.array(y_interp), np.array(z_interp))

def cumulative_distance(*args):
    """calculate the cumlative linear distances along a path of x,y coordinates, including the first point (the returned array always starts with zero)
    args:
        either a single iterable of points or separate iterables containing
        the coordinates. Can accept two or three dimensional input.
    returns:
        dist - array, cumulative distance along the points"""

    largs = len(args)
    if(largs == 1):
        args = zip(*args[0])
        if(len(args) == 2):
            return(_cum_dist_2D(args[0], args[1]))
        elif(len(args) == 3):
            return(_cum_dist_3D(args[0], args[1], args[2]))
        else:
            raise(subcalc_class.EMFError("Input error. cumulative_distance accepts either a single iterable of points or separate iterables containing the coordinates. Can accept two or three dimensional input."))
    elif(largs == 2):
        if(len(args[0]) != len(args[1])):
            raise(subcalc_class.EMFError("Inputs must have the same length."))
        return(_cum_dist_2D(args[0], args[1]))
    elif(largs == 3):
        if(not (len(args[0]) == len(args[1]) == len(args[2]))):
            raise(subcalc_class.EMFError("Inputs must have the same length."))
        return(_cum_dist_3D(args[0], args[1], args[2]))
    else:
        raise(subcalc_class.EMFError("Cannot accept more than three inputs."))

def _cum_dist_2D(x, y):
    d = np.zeros((len(x),), dtype=float)
    for i in range(1, len(x)):
        d[i] = _dist(x[i], x[i-1], y[i], y[i-1]) + d[i-1]
    return(d)

def _cum_dist_3D(x, y, z):
    d = np.zeros((len(x),), dtype=float)
    for i in range(1, len(x)):
        d[i] = _dist3(x[i], x[i-1], y[i], y[i-1], z[i], z[i-1]) + d[i-1]
    return(d)
