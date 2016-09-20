from .. import np, pd, os, interpn

from ..emf_class import EMFError

import subcalc_funks

class Model(object):

    def __init__(self, X, Y, B, *args):
        """Grid data must be passed in
        args:
            X - 2D array of x coordinates
            Y - 2D array of y coordinats
            B - 2D array of magnetic field magnitudes
            info - dictionary of model information"""
        #2D reference grid arrays
        self.X, self.Y, self.B = X, Y, B
        #grid dimensions for reference
        self.grid_limits = {'xmax': np.max(self.x), 'xmin': np.min(self.x),
                                'ymax': np.max(self.y), 'ymin': np.min(self.y)}
        #dictionary of reference grid information from the SubCalc model
        if(len(args) > 0):
            info = args[0]
            if(type(info) is not dict):
                raise(EMFError("""
                The fourth argument to Model.__init__ must be a dictionary
                of model result information"""))
            self.info = info
        #other reference objects in the model, like substation boundaries,
        #stored in a list of Footprint objects
        self.footprint_df = None
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
    north_angle = property(_get_north_angle, _set_north_angle, None,
            """Angle of the Northern direction in degrees, where 0 represents
            the vertical or Y direction and clockwise represents increasing
            angle""")

    def load_footprints(self, footprint_info, **kwargs):
        """Read footprint data from a csv file and organize it in
        Footprint objects stored in self.footprints
        args:
            footprint_info - string, path to the footprint csv/excel data.
                            If footprint data is in an excel workbook with,
                            multiple sheets, the sheet name must be passed
                            to the kwarg 'sheet'

                                or

                            an existing DataFrame with footprint data"""
        #load file if footprint_info is not a DataFrame
        if(not (type(footprint_info) is pd.DataFrame)):
            #check extension
            footprint_info = subcalc_funks._check_extension(footprint_info,
                'csv', 'Footprint files must be csv files.')
            #load data
            df = pd.read_csv(footprint_info)
        else:
            df = footprint_info
        #store the DataFrame
        self.footprint_df = df
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
            %s""" % (str(cols), footprint_info,
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
        #store grid x,y extents
        xmin, xmax = self.grid_limits['xmin'], self.grid_limits['xmax']
        ymin, ymax = self.grid_limits['ymin'], self.grid_limits['ymax']
        #get 1D vectors of x and y coordinates to resample at, from kwargs
        if('N' in kwargs):
            N = kwargs['N']
            if(not subcalc_funks._is_int(N)):
                raise(EMFError('Keyword argument "N" must be a whole number'))
            N = float(N)
            aspect = float(self.B.shape[0])/self.B.shape[1]
            N_y = np.ceil(np.sqrt(N/aspect))
            N_x = np.ceil(N/N_y)
            x = np.linspace(xmin, xmax, N_x)
            y = np.linspace(ymin, ymax, N_y)
        else:
            if('x' in kwargs):
                x = kwargs['x']
                if(subcalc_funks._is_int(x)):
                    x = np.linspace(xmin, xmax, x)
            else:
                x = self.x
            if('y' in kwargs):
                y = kwargs['y']
                if(subcalc_funks._is_int(y)):
                    y = np.linspace(ymin, ymax, y)
            else:
                y = self.y
        #flip y coordinates so that Y prints intuitively
        y = y[::-1]
        #resample the grid
        X, Y = np.meshgrid(x, y)
        #some arrays have to be flipped to conform to conventions of interpn
        B_resample = interpn((self.y[::-1], self.x),
                                self.B[::-1,:],
                                (Y[::-1,:], X))

        #return with re-flipped results
        return(X, Y, B_resample[::-1,:])

    def export(self, **kwargs):
        """Export the grid data and accompanying info to an excel file with
        two tabs, one for the grid and a second for the info
        kwargs:
            path - string, output destination/filename for workbook"""
        #get appropriate export filename
        fn = os.path.basename(self.info['REF_path'])
        fn = subcalc_funks._path_manage(fn, '.xlsx', **kwargs)
        #create excel writing object
        xl = pd.ExcelWriter(fn, engine='xlsxwriter')
        #write grid data
        pd.DataFrame(self.B, columns=self.x, index=self.y).to_excel(
                xl, sheet_name='grid')
        #write model information
        pd.DataFrame([self.info[k] for k in self.info], index=self.info.keys(),
                columns=['Parameter Value']).sort_index().to_excel(
                        xl, sheet_name='info', index_label='Parameter Name')
        #write footprint DataFrame if present
        if(self.footprint_df is not None):
            self.footprint_df.to_excel(xl, sheet_name='footprints', index=False)
        #save and print
        xl.save()
        print('model saved to: %s' % fn)

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
