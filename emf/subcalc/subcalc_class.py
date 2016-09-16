from .. import np
from .. import pd

from ..emf_class import EMFError

import subcalc_funks

class Model(object):

    def __init__(self, data, info):
        #dictionary of reference grid information from the SubCalc model
        self.info = info
        #sets the 'component' of the magnetic field results stored by the model
        #   can be 'Bx', 'By', 'Bmax', or 'Bres'
        self.B_key = 'Bmax'
        #2D reference grid arrays
        self.X, self.Y, self.B = self._meshgrid(data)
        #grid dimensions for reference
        self.grid_limits = {'xmax': np.max(self.x), 'xmin': np.min(self.x),
                                'ymax': np.max(self.y), 'ymin': np.min(self.y)}
        #other reference objects in the model, like substation boundaries,
        #stored in a list of Footprint objects
        self.footprints = []
        self.footprint_groups = [] #list of lists of integers for grouping
        #angle of the northern direction with respect to the grid
        #   where 0 degrees is the positive y axis and clockwise is increasing
        self._north_angle = None

    #north angle property
    #properties
    def _get_x(self):
        return(self.X[0,:])
    x = property(_get_x, None, None, 'Unique x values in model grid')

    def _get_y(self):
        return(self.Y[:,0])
    y = property(_get_y, None, None, 'Unique y values in model grid')

    def _get_north_angle(self):
        return(self._north_angle)
    def _set_north_angle(self, angle):
        if(not subcalc_funks._is_number(angle)):
            raise(EMFError("""
            The 'north_angle' attribute of a Model object must be a number."""))
        else:
            self._north_angle = float(angle)
    def _del_north_angle(self):
        del(self._north_angle)
    north_angle = property(_get_north_angle, _set_north_angle, _del_north_angle,
            """Angle of the Northern direction in degrees, where 0 represents
            the vertical or Y direction and clockwise represents increasing
            angle""")

    def load_footprints(self, footprint_path, **kwargs):
        """Read footprint data from a csv file and organize it in
        Footprint objects stored in self.footprints
        args:
            footprint_path - string, path to the footprint csv/excel data.
                            If footprint data is in an excel workbook with,
                            multiple sheets, the sheet name must be passed
                            to the kwarg 'sheet'"""
        #check extension
        footprint_path = subcalc_funks._check_extension(footprint_path, 'csv',
            """Footprint files must be csv files.""")
        #load data
        df = pd.read_csv(footprint_path)
        #check all columns are present
        cols = ['Group', 'Name', 'X', 'Y', 'Power Line?', 'Of Concern?',
                'Draw as Loop?', 'Group']
        if(set(cols) - set(df.columns)):
            raise(EMFError("""
            Footprint csv files must have the following column names:
            %s
            The footprint csv file with path:
            %s
            is missing or misspells these columns:
            %s""" % (str(cols), footprint_path,
                    str(list(set(cols) - set(df.columns))))))
        message = """
            The column:
                "%s"
            must contain only one repeated value for each footprint.
            It contained multiple values for footprint name:
                "%s" """
        fields = cols[4:]
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
        #update things (footprint_groups)
        self.update()

    def cross_section(self, p1, p2, **kwargs):
        """Interpolate the field along a line between two points
        args:
            p1 - iterable, an x-y pair
            p2 - iterable, an x-y pair
        kwargs:
            n - integer, number of points sampled (default 1000)
        returns:
            x - array, x coordinates of interpolated values
            y - array, y coordinates of interpolated values
            B_interp - array, interpolated field values"""
        #check point lengths
        if((len(p1) != 2) or (len(p2) != 2)):
            raise(EMFError('Points must consist of two values, xy pairs.'))
        #check kwargs
        if('n' in kwargs):
            n = kwargs['n']
        else:
            n = 1000
        #create x and y vectors
        x, y = np.linspace(p1[0], p2[0], 1000), np.linspace(p1[1], p2[1], n)
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
        if((x > self.grid_limits['xmax']) or
            (x < self.grid_limits['xmin']) or
                (y > self.grid_limits['ymax']) or
                    (y < self.grid_limits['ymin'])):
            return(False)
        else:
            return(True)

    def resample(self, **kwargs):
        """Resample the model grid along a new number of x,y values, a new
        selection of x,y values, or a new number of total values
        kwargs:
            x - int or iterable, new number of x samples or new selection of
                x samples
            y - int or iterable, new number of y samples or new selection of
                y samples
            N - int, new approximate number of total samples
                    - overrides x and y kwargs
                    - preserves approx ratio of number of x and y values
                    - rounds up to nearest possible whole number of points
        returns:
            X - 2D numpy array with resampled x coordinates
            Y - 2D numpy array with resampled y coordinates
            B_resample - 2D numpy array with resampled field values"""
        #get x and y coordinates to resample from arguments
        if('N' in kwargs):
            N = kwargs['N']
            if(type(N) is not int):
                raise(EMFError('Keyword argument "N" must be an integer'))
            N = float(N)
            aspect = float(self.B.shape[0])/self.B.shape[1]
            N_y = np.ceil(np.sqrt(N/aspect))
            N_x = np.ceil(N/N_y)
            x = np.linspace(self.grid_limits['xmin'],
                            self.grid_limits['xmax'], N_x)
            y = np.linspace(self.grid_limits['ymin'],
                            self.grid_limits['ymax'], N_y)
        else:
            if('x' in kwargs):
                x = kwargs['x']
                if(type(x) is int):
                    x = np.linspace(self.grid_limits['xmin'],
                                    self.grid_limits['xmax'], x)
                else:
                    x = self.x
            if('y' in kwargs):
                y = kwargs['y']
                if(type(y) is int):
                    y = np.linspace(self.grid_limits['ymin'],
                                    self.grid_limits['ymax'], y)
                else:
                    y = self.y
        #flip y coordinates so that Y prints intuitively
        y = y[::-1]
        #resample the grid
        X, Y = np.meshgrid(x, y)
        if(False):
            B_resample = np.empty((len(y), len(x)))
            for i in range(len(y)):
                for j in range(len(x)):
                    B_resample[i,j] = subcalc_funks._bilinear_interp(self,
                            x[j], y[i])
        return(X, Y, B_resample)

    def _meshgrid(self, data):
        """Convert raw grid data read from a SubCalc output file
        (by subcalc_funks.read_REF) into meshed grids of X, Y coordinates
        and their corresponding B field values
        args:
            data - dict, keyed by 'x', 'y', and self.B_key with 1D arrays
                    in each, of equal length
        returns:
            X - 2D array, x coordinates
            Y - 2D array, y coordinates
            B - 2D array, magnetic field values"""

        #find the number of points in a row
        x = data['x']
        y = data['y']
        b = data[self.B_key]
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

class Footprint(object):

    def _get_x(self):
        if(self.draw_as_loop):
            return(self._x + [self._x[0]])
        else:
            return(self._x)
    def _set_x(self, value):
        self._x = value
    def _del_x(self):
        del(self._x)
    x = property(_get_x, _set_x, _del_x, 'x coordinates of Footprint vertices')

    def _get_y(self):
        if(self.draw_as_loop):
            return(self._y + [self._y[0]])
        else:
            return(self._y)
    def _set_y(self, value):
        self._y = value
    def _del_y(self):
        del(self._y)
    y = property(_get_y, _set_y, _del_y, 'y coordinates of Footprint vertices')

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
            tag - anything, identifier used to group footprints
            label_ha - string
            label_va - string"""
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

    def __str__(self):
        """quick and dirty printing"""
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            s += str(k) + ': ' + str(v[k]) + '\n'
        return(s)
