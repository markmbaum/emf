import numpy as np
import pandas as pd

import subcalc_funks

from ..emf_class import EMFError

class Footprint(object):

    def get_x(self):
        if(self.draw_as_loop):
            return(self._x + [self._x[0]])
        else:
            return(self._x)
    x = property(get_x)

    def get_y(self):
        if(self.draw_as_loop):
            return(self._y + [self._y[0]])
        else:
            return(self._y)
    y = property(get_y)

    def __init__(self, name, x, y, power_line, of_concern, draw_as_loop, group):
        """
        args:
            name - string, the name of the Footprint, i.e. "Substation"
            x - iterable, x coordinates of Footprint
            y - iterable, y coordinates of Footprint
            of_concern - bool, True if the Footprint represents an area
                        that is potentially concerned about EMF (homes)
            draw_as_loop - bool, True if the footprint should be plotted
                            as a closed loop
            tag - anything, identifier used to group footprints"""
        #check x and y are the same length
        if(len(x) != len(y)):
            raise(EMFError("""
            Footprints must have the same number of x and y values to form
            spatial coordinates"""))
        #set attributes
        self.name = name #string
        self._x = list(x)
        self._y = list(y)
        self.power_line = power_line #bool
        self.of_concern = of_concern #bool
        self.draw_as_loop = draw_as_loop #bool
        self.group = group #string

class Model(object):

    def __init__(self, data, grid_info):
        #dictionary of reference grid information from the SubCalc model
        self.grid_info = grid_info
        #2D reference grid arrays
        self.X, self.Y, self.B = self._meshgrid(data, 'Bmax')
        #    Sets which component of the results are used ^^^
        #1D grid row and column coordinates
        self.x = self.X[0,:]
        self.y = self.Y[:,0]
        #other reference objects in the model, like substation boundaries,
        #stored in a list of Footprint objects
        self.footprints = []
        self.footprint_groups = [] #list of lists of integers for grouping

    def load_footprints(self, footprint_path):
        """Read footprint data from a csv file and organize it in
        Footprint objects stored in self.footprints
        args:
            footprint_path - string, path to the footprint csv data"""
        #check extension
        footprint_path = subcalc_funks._check_extension(footprint_path, 'csv',
            """Footprint files must be csv files.""")
        #load data
        df = pd.read_csv(footprint_path)
        message = """
            The column:
                "%s"
            must contain only one repeated value for each footprint.
            It contained multiple values for footprint name:
                "%s" """
        fields = ['Power Line?', 'Of Concern?', 'Draw as Loop?', 'Group']
        #create a footprint for unique entries in the 'Name' field
        for n in df['Name'].unique():
            s = df[df['Name'] == n]
            #check that certain fields only contain a single entry
            for f in fields:
                if(len(s[f].unique()) > 1):
                    print s[f].unique()
                    raise(EMFError(message % (f,n)))
            fp = Footprint(n, s['X'].values, s['Y'].values,
                    bool(s['Power Line?'].unique()[0]),
                    bool(s['Of Concern?'].unique()[0]),
                    bool(s['Draw as Loop?'].unique()[0]),
                    s['Group'].unique()[0])
            self.footprints.append(fp)
        #update things
        self.update()

    def cross_section(self, p1, p2):
        """Interpolate the field along a line between two points
        args:
            p1 - iterable, an x-y pair
            p2 - iterable, an x-y pair
        returns:
            x - array, x coordinates of interpolated values
            y - array, y coordinates of interpolated values
            B_interp - array, interpolated field values"""
        #check point lengths
        if((len(p1) != 2) or (len(p2) != 2)):
            raise(EMFError('Points must consist of two values, an xy pair.'))
        #create x and y vectors
        x, y = np.linspace(p1[0], p2[0], 1000), np.linspace(p1[1], p2[1], 1000)
        B_interp = self.interp(x, y)
        return(x, y, B_interp)


    def interp(self, x, y):
        """Interpolate in the x and y directions to find an estimated B
        value at an x,y location within the model
        args:
            x - iterable or scalar, x coordinate(s) to interpolate at
            y - iterable or scalar, y coordinate(s) to interpolate at
        return:
            B_interp - array or float, the interpolated field value"""
        #make x,y iterable if scalars are passed in
        if(not (hasattr(x, '__len__') and hasattr(y, '__len__'))):
            scalar = True
            x = np.array([x], dtype = float)
            y = np.array([y], dtype = float)
        else:
            scalar = False
        #check that all points are in the grid
        if(not all([self.in_grid(x[i],y[i]) for i in range(len(x))])):
            raise(EMFError("""
            x,y coordinates must fall inside the reference grid"""))
        #interpolate
        B_interp = np.array([subcalc_funks._bilinear_interp(self, x[i], y[i])
                                for i in range(len(x))])
        #return
        if(scalar):
            return(B_interp[0])
        else:
            return(B_interp)

    def in_grid(self, x, y):
        """Check if an x,y coordinate pair is inside the Model grid
        args:
            x - float, x coordinate
            y - float, y coordinate
        returns:
            b - bool, True if x,y is in the grid, False if it's not"""
        if((x > self.grid_info['Maximum X Coordinate']) or
            (x < self.grid_info['Minimum X Coordinate']) or
            (y > self.grid_info['Maximum Y Coordinate']) or
            (y < self.grid_info['Minimum Y Coordinate'])):
            return(False)
        else:
            return(True)

    def _meshgrid(self, data, B_key):
        """Convert raw grid data read from a SubCalc output file
        (by subcalc_funks.read_REF) into meshed grids of X, Y coordinates
        and their corresponding B field values
        args:
            data - dict, keyed by 'x', 'y', and 'B' with 1D arrays in
                    each, of equal length
            B_key - string, can be 'Bx', 'By', 'B', or 'Bres' and is
                    used to select which component of the magnetic field
                    results are stored by the model
        returns:
            X - 2D array, x coordinates
            Y - 2D array, y coordinates
            B - 2D array, magnetic field values"""

        #find the number of points in a row
        x = data['x']
        y = data['y']
        b = data[B_key]
        count = 0
        v = y[count]
        while(y[count+1] == y[count]):
            count += 1
        count += 1
        #get ncols and nrows
        L = len(x)
        ncols = count
        nrows = L/ncols
        #fill 2D arrays
        X, Y, B = [np.empty((nrows,ncols)) for i in range(3)]
        count = 0
        for i in range(nrows):
            X[i,:] = x[count:count+ncols]
            Y[i,:] = y[count:count+ncols]
            B[i,:] = b[count:count+ncols]
            count += ncols

        return(X, Y, B)

    def update(self):
        """Execute all update methods"""
        self._update_tag_groups()

    def _update_tag_groups(self):
        """Generate a list of lists of Footprint indices with identical tags"""
        u = list(set([f.group for f in self.footprints]))
        self.footprint_groups = [[] for i in range(len(u))]
        for i in range(len(self.footprints)):
            self.footprint_groups[u.index(self.footprints[i].group)].append(i)
