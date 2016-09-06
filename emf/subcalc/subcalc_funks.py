from .. import np

from ..emf_funks import (_path_manage, _check_extension, _is_number, _check_intable,
                _flatten, _sig_figs)

import subcalc_class

def load_model(*args):
    """Read a .REF output file and load the data into a Model object
    args:
        REF_path - string, path to the output .REF file of field results
        footprint_path - string, optional, path to the csv file of
                         footprint data
    returns
        mod - Model object containing results"""

    #load field results into a model
    data, grid_info = read_REF(args[0])
    mod = subcalc_class.Model(data, grid_info)

    if(len(args) > 1):
        #load footprints into the model
        mod.load_footprints(args[1])

    return(mod)

def read_REF(file_path):
    """Reads a .REF output file generated by the SUBCALC program and pulls
    out information about the reference grid of the model with the Res and
    Max magnetic fields.
    args:
        file_path - string, path to saved .REF output file
    returns:
        data - dict, keys are 'x', 'y', 'Bmax', 'Bres', 'Bx', and 'By'
        info - dict, reference grid and other information"""

    #check the extension
    file_path = _check_extension(file_path, 'REF', """
        SubCalc results are saved to text files with .REF extensions.
        The input path:
            "%s"
        does not have the correct extension.""" % file_path)

    #allocate dictionaries
    info = {} #dictionary storing reference grid information
    keys = ['X Coord', 'Y Coord', 'X Mag', 'Y Mag', 'Max', 'Res']
    return_keys = ['x', 'y', 'Bx', 'By', 'Bmax', 'Bres']
    data = dict(zip(keys, [[] for i in range(len(keys))]))

    #pull data out
    with open(file_path, 'r') as ifile:
        #store information about the grid
        for i in range(24):
            line = ifile.readline().strip()
            if(':' in line):
                idx = line.find(':')
                line = [line[:idx], line[idx+1:]]
                if(_is_number(line[1])):
                    info[line[0]] = float(line[1])
                else:
                    info[line[0]] = line[1].strip()
        #read through the rest of the data
        for line in ifile:
            for k in keys:
                if(k == line[:len(k)]):
                    L = line[line.index(':')+1:]
                    data[k].append([float(i) for i in L.split()])

    #flatten the lists in data
    for k in data:
        data[k] = np.array(_flatten(data[k]))

    #switch the keys
    data = dict(zip(return_keys, [data[k] for k in keys]))

    return(data, info)

def _bilinear_interp(mod, x, y):
    """Use Model results to interpolate linearly in two dimensions for an
    estimate of any x,y coordinate inside the grid.
    args:
        mod - Model object
        x - float, x coordinate to interpolate at
        y - float, y coordinate to interpolate at
    returns:
        B_interp - float, interpolated field value"""
    #first find the 4 point grid cell containing x,y
    #   (the point is assumed to lie inside the grid)
    _, xidx = _double_min(np.abs(mod.x - x))
    _, yidx = _double_min(np.abs(mod.y - y))
    #get coordinates and values
    x1, x2 = mod.x[xidx]
    y1, y2 = mod.y[yidx]
    B11 = mod.B[yidx[0], xidx[0]]
    B12 = mod.B[yidx[0], xidx[1]]
    B21 = mod.B[yidx[1], xidx[0]]
    B22 = mod.B[yidx[1], xidx[1]]
    #interpolate
    B_interp = (1.0/((x2 - x1)*(y2 - y1)))*(
        B11*(x2 - x)*(y2 - y) + B21*(x - x1)*(y2 - y)
        + B12*(x2 - x)*(y - y1) + B22*(x - x1)*(y - y1))

    return(B_interp)

def _2Dmax(G):
    """Find the indices of the maximum value in a 2 dimensional array
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
    mins = np.array([m, m], dtype = float)
    idxs = np.zeros((2,), dtype = int)
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
