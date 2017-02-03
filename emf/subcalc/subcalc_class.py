from .. import np, pd, os, copy, datetime, textwrap, _interpn

from ..emf_class import EMFError
from ..emf_funks import _print_str_list

import subcalc_funks
import subcalc_calcs
import subcalc_print

class Model(object):
    """Model objects store Tower and/or Conductor objects, information about the desired model grid, and provide the means of computing model results and returning them in a Results object. Model objects are designed to make it easy to compute the fields in a uniform height 2 dimensional grid. Model objects store x limits, y limits, a z coordinate/height, and a grid spacing to automatically construct a 2d grid of sample points and calculate magnetic fields at those points. The calculate() method returns a grid of results in the form of a Results object, which can be used to create plots and analyze the grid. Alternatively, an arbitrary set of x, y, z coordinates can be passed to the sample() method to calculate fields at any set of points in 3d space and bypass the Results object entirely."""

    def __init__(self, **kw):
        """
        kw:
            name - string, a name identifying the Model, used for filenames
            towers - an iterable of Tower objects to include in the Model
            conductors - an iterable of Conductor objects to include in the Model
            xmin - float, the minimum x value in the Model area, default is zero (ft)
            xmax - float, the maximum x value in the Model area, will use
                   Tower and Conductor information to pick value automatically
                   if left unset and return None if left unset with no Towers or
                   Conductors in the Model
            ymin - float, the minimum y value in the Model area, default is zero (ft)
            ymax - float, the maximum y value in the Model area, will use
                   Tower and Conductor information to pick value automatically
                   if left unset and return None if left unset with no Towers or
                   Conductors in the Model
            z - float, the vertical height to calculate magnetic fields at (ft),
                defaults to 3.28 ft (the SUBCALC default)
            spacing - float, the spacing between grid points, defaults to 1 (ft)"""

        self._name = None
        self._towers = []
        self._conductors = []
        self._xmin = None
        self._xmax = None
        self._ymin = None
        self._ymax = None
        self._z = 3.28
        self._spacing = None

        self._auto_lim_frac = 0.1

        if('name' in kw):
            self.name = kw['name']
        if('towers' in kw):
            for i in kw['towers']:
                self.add_tower(i)
        if('conductors' in kw):
            for i in kw['conductors']:
                self.add_conductor(i)
        if('xmin' in kw):
            self.xmin = kw['xmin']
        if('xmax' in kw):
            self.xmax = kw['xmax']
        if('ymin' in kw):
            self.ymin = kw['ymin']
        if('ymax' in kw):
            self.ymax = kw['ymax']
        if('z' in kw):
            self.z = kw['z']
        if('spacing' in kw):
            self.spacing = kw['spacing']

    #---------------------------------------------------------------------------
    #properties

    def _get_name(self):
        if(not self._name):
            return('unnamed-model')
        else:
            return(self._name)
    def _set_name(self, value):
        if(not value):
            raise(EMFError("Cannot set the 'name' property of a Model object to an implicitly false object. Use a string."))
        self._name = str(value)
    name = property(_get_name, _set_name, None, 'A string identifying the Model, used for filenames')

    def _get_towers(self): return(self._towers)
    towers = property(_get_towers, None, None, 'List of Tower objects in the Model')

    def _get_tower_groups(self):
        """Generate a list of lists of Tower objects with identical group
        strings, with the sublists sorted according to the Tower seq property"""
        u = list(self.tower_group_names)
        groups = [[] for i in range(len(u))]
        #create the groups
        for i in range(len(self.towers)):
            t = self.towers[i]
            groups[u.index(t.group)].append(t)
        #sort the groups by seq
        k = lambda t: t.seq
        for i in range(len(groups)):
            groups[i] = sorted(groups[i], key=k)
        return(groups)
    tower_groups = property(_get_tower_groups, None, None, 'A list of lists of Tower objects with identical group strings, with the sublists sorted according to the Tower seq property')

    def _get_tower_group_names(self):
        return(set([t.group for t in self.towers]))
    tower_group_names = property(_get_tower_group_names, None, None, 'Set of unique group strings representing the Towers in Model.towers')

    def _get_conductors(self): return(self._conductors)
    conductors = property(_get_conductors, None, None, 'List of Conductor objects in the Model')

    def _get_xmin(self):
        if(self._xmin is None):
            xr, yr = self._xy_ranges()
            if((xr is not None) and (yr is not None)):
                m = max([xr[1] - xr[0], yr[1] - yr[0]])
                pad = m*self.auto_lim_frac
                return(xr[0] - pad)
            else:
                return(None)
        else:
            return(self._xmin)
    def _set_xmin(self, value):
        if(value is None):
            self._xmin = None
        elif(float(value) >= self.xmax):
            raise(EMFError('xmin must be less than xmax.'))
        else:
            self._xmin = float(value)
    xmin = property(_get_xmin, _set_xmin, None, 'The minimum x value in the Model area in feet, will use Tower and Conductor information to pick a value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_xmax(self):
        if(self._xmax is None):
            xr, yr = self._xy_ranges()
            if((xr is not None) and (yr is not None)):
                m = max([xr[1] - xr[0], yr[1] - yr[0]])
                pad = m*self.auto_lim_frac
                return(xr[1] + pad)
            else:
                return(None)
        else:
            return(self._xmax)
    def _set_xmax(self, value):
        if(value is None):
            self._xmax = None
        elif(float(value) <= self.xmin):
            raise(EMFError('xmax must be greater than xmin.'))
        else:
            self._xmax = float(value)
    xmax = property(_get_xmax, _set_xmax, None, 'The maximum x value in the Model area in feet, will use Tower and Conductor information to pick a value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_xlim(self): return(self.xmin, self.xmax)
    def _set_xlim(self, value): self.xmin, self.xmax = value
    xlim = property(_get_xlim, _set_xlim, None, 'The range of x coordinates in the Model: (xmin, xmax)')

    def _get_ymin(self):
        if(self._ymin is None):
            xr, yr = self._xy_ranges()
            if((xr is not None) and (yr is not None)):
                m = max([xr[1] - xr[0], yr[1] - yr[0]])
                pad = m*self.auto_lim_frac
                return(yr[0] - pad)
            else:
                return(None)
        else:
            return(self._ymin)
    def _set_ymin(self, value):
        if(value is None):
            self._ymin = None
        elif(float(value) >= self.ymax):
            raise(EMFError('ymin must be less than ymax.'))
        else:
            self._ymin = float(value)
    ymin = property(_get_ymin, _set_ymin, None, 'The minimum y value in the Model area in feet, will use Tower and Conductor information to pick a value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_ymax(self):
        if(self._ymax is None):
            xr, yr = self._xy_ranges()
            if((xr is not None) and (yr is not None)):
                m = max([xr[1] - xr[0], yr[1] - yr[0]])
                pad = m*self.auto_lim_frac
                return(yr[1] + pad)
            else:
                return(None)
        else:
            return(self._ymax)
    def _set_ymax(self, value):
        if(value is None):
            self._ymax = None
        elif(float(value) <= self.xmin):
            raise(EMFError('ymax must be greater than ymin.'))
        else:
            self._ymax = float(value)
    ymax = property(_get_ymax, _set_ymax, None, 'The maximum y value in the Model area in feet, will use Tower and Conductor information to pick a value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_ylim(self): return(self.ymin, self.ymax)
    def _set_ylim(self, value): self.ymin, self.ymax = value
    ylim = property(_get_ylim, _set_ylim, None, 'The range of y coordinates in the Model: (ymin, ymax)')

    def _get_auto_lim_frac(self): return(self._auto_lim_frac)
    def _set_auto_lim_frac(self, value):
        if(not subcalc_funks._is_number(value)):
            raise(EMFError("auto_lim_frac must be a nubmer greater than or equal to zero."))
        if(value < 0):
            raise(EMFError("auto_lim_frac must be a nubmer greater than or equal to zero."))
        self._auto_lim_frac = float(value)
    auto_lim_frac = property(_get_auto_lim_frac, _set_auto_lim_frac, None, 'When xmin, xmax, ymin, and ymax are left unset, they are automatically determined by the ranges of x and y coordinates of Tower and Conductor objects in the Model. They are set to the edges of those ranges with some padding. auto_lim_frac determines how much padding is automatically used between the limits of the Model domain and the wires in the Model.')

    def _get_z(self): return(self._z)
    def _set_z(self, value):
        self._z = float(value)
    z = property(_get_z, _set_z, None, 'Vertical (z) coordinate to calculate magnetic fields at (ft)')

    def _get_spacing(self):
        if(self._spacing is None):
            return(1.0)
        else:
            return(self._spacing)
    def _set_spacing(self, value):
        s = float(value)
        if(s <= 0):
            raise(EMFError)
        self._spacing = float(value)
    spacing = property(_get_spacing, _set_spacing, None, 'Spacing between sample points along x and y axes of the Model grid (ft)')

    def _get_N(self): return(len(self.x)*len(self.y))
    N = property(_get_N, None, None, 'The total number of sample points in the Model. Cannot be set directly (use Model.spacing to adjust).')

    def _get_x(self):
        if((self.xmax is None) or (self.xmin is None)):
            return(None)
        inc = self.spacing/2.0
        x = np.arange(self.xmin, self.xmax + inc, self.spacing, dtype=float)
        return(x)
    x = property(_get_x, None, None, 'Array of unique x coordinates in the Model grid')

    def _get_y(self):
        if((self.ymax is None) or (self.ymin is None)):
            return(None)
        inc = self.spacing/2.0
        y = np.arange(self.ymin, self.ymax + inc, self.spacing, dtype=float)
        return(y[::-1])
    y = property(_get_y, None, None, 'Array of unique y coordinates in the Model grid')

    def _get_X(self):
        x, y = self.x, self.y
        if((x is None) or (y is None)):
            return(None)
        X = np.tile(x, (len(y),1))
        return(X)
    X = property(_get_X, None, None, '2D array of x coordinates in the Model Grid')

    def _get_Y(self):
        x, y = self.x, self.y
        if((x is None) or (y is None)):
            return(None)
        Y = np.tile(np.reshape(y, (len(y),1)), (1,len(x)))
        return(Y)
    Y = property(_get_Y, None, None, '2D array of y coordinates in the Model Grid')

    def _get_segments(self):
        """Parse lists of Conductor and Tower objects in the Model into individual wire segments, for emf calculations. A list of tuples is returned, each representing a single wire and containing the wire's start and end points (a and b), the current amplitude, and the phase:

            ((xa, ya, za), (xb, yb, zb), I, phase)"""
        segs = []
        #create segments from the list of Towers
        tg = self.tower_groups
        if(any([(len(g) < 2) for g in tg])):
            raise(EMFError('Tower groups must have more than one tower in them. A single Tower cannot be parsed into wire segments.'))
        for g in tg:
            #parse to segments
            for i in range(len(g) - 1):
                t1, t2 = g[i], g[i+1]
                x1, y1, z1 = t1.conductor_x, t1.conductor_y, t1.conductor_z
                x2, y2, z2 = t2.conductor_x, t2.conductor_y, t2.conductor_z
                for j in range(len(t1)):
                    segs.append(
                        (
                            np.array([x1[j], y1[j], z1[j]]),
                            np.array([x2[j], y2[j], z2[j]]),
                            t1.I[j],
                            t1.phase[j]
                        )
                    )
        #create segments from Conductors
        for c in self.conductors:
            for i in range(len(c) - 1):
                segs.append(
                    (
                        np.array([c.x[i], c.y[i], c.z[i]]),
                        np.array([c.x[i+1], c.y[i+1], c.z[i+1]]),
                        c.I,
                        c.phase
                    )
                )

        return(segs)
    segments = property(_get_segments, None, None, "Parse lists of Conductor and Tower objects in the Model into individual wire segments, for emf calculations. A list of tuples is returned, each representing a single wire and containing the wire's start and end points (a and b), the current amplitude, and the phase: ((xa, ya, za), (xb, yb, zb), I, phase)")

    def _get_footprints(self):
        fps = []
        #create Footprints from Tower objects
        for g in self.tower_groups:
            fps.append(Footprint(g[0].group, g[0].group,
                    [t.tower_x for t in g], [t.tower_y for t in g],
                    True, False, False))
        #create Footprints from Conductor objects
        for c in self.conductors:
            fps.append(Footprint(c.name, c.name, c.x, c.y, True, False, False))
        return(fps)
    footprints = property(_get_footprints, None, None, 'List of Footprint objects create from Tower and Conductor objects in the Model.')

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(subcalc_print._str_Model(self))

    def add_tower(self, tnew):
        """Add a Tower object to the Model
        args:
            tnew - Tower object to add (a copy is added)"""
        #check type
        if(type(tnew) is not Tower):
            raise(EMFError('The add_tower method of Model objects can only add Tower objects.'))
        #check that the Tower's number of conductors, configuration of currents,
        #and configuration of phases matches those of conductors in the same
        #group
        group_names = self.tower_group_names
        if(tnew.group in group_names):
            #collect Towers in the same group
            group = [t for t in self.towers if (t.group == tnew.group)]
            #get an existing Tower in the same group to compare against
            tcomp = group[0]
            #check lengths
            if(len(tnew) != len(tcomp)):
                raise(EMFError('All Towers with the same group must have the same length. The Tower cannot be added because it has a different number of wires than other Towers in the "%s" group.' % tnew.group))
            #check currents
            if(np.any(tnew.I != tcomp.I)):
                raise(EMFError('All Towers with the same group must have the same currents. The Tower cannot be added because it has a different current configuration than other Towers in the "%s" group.' % tnew.group))
            #check phasing
            if(np.any(tnew.phase != tcomp.phase)):
                raise(EMFError('All Towers with the same group must have the same phasing. The Tower cannot be added because it has a different phase configuration than other Towers in the "%s" group.' % tnew.group))
            #check sequencing
            if(tnew.seq in [t.seq for t in group]):
                raise(EMFError('All Towers with the same group must have the different "seq" numbers. The Tower cannot be added because it has the same "seq" as another Tower in the "%s" group.' % tnew.group))

        self._towers.append(copy.deepcopy(tnew))

    def load_towers(self, fn, **kw):
        """Read a subcalc tower template and create a new list of Tower objects, adding each of the new towers to the Model's existing list of Towers.
        args:
            fn - string, the path string to the subcalc tower template file (excel or csv)
        kw:
            sheet - string, if the template file is an excel workbook with multiple sheets,
                    the target sheet must be specified"""
        #use the funks function
        for t in subcalc_funks.load_towers(fn, return_model=False, **kw):
            self.add_tower(t)

    def add_conductor(self, c):
        """Add a Conductor object to the Model
        args:
            c - Conductor object to add (a copy is added)"""
        if(type(c) is not Conductor):
            raise(EMFError('The add_conductor method of Model objects can only add Conductor objects.'))
        if(c.name in [i.name for i in self.conductors]):
            raise(EMFError('Cannot add conductors with identical names. The name "%s" is already used.' % c.name))
        self._conductors.append(copy.deepcopy(c))

    def calculate(self, components='all', footprints=None, clear=False, **kw):
        """Parse the Tower and Conductor objects into individual wire segements and compute the magnetic fields produced by the collection of wire segments in the model domain, returning a Results object. The Tower and Conductor objects will also be converted to Footprint objects in the newly created Results object.
        optional args:
            components - list of components of the field to be included in
                         the returned Results object. Can be any combination
                         of 'Bmax', 'Bres', 'Bx', 'By', 'Bz' or it can be
                         'all' to store all of those options. The default is
                         'all'. Choosing fewer components will not speed up
                         calculations, but it will reduce memory demand. The
                         components chosen will become the available 'Bkeys'
                         in the returned Results object.
            footprints - string, path to the footprint csv/excel data.
                        If footprint data is in an excel workbook with
                        multiple sheets, the sheet name must be passed
                        to the kwarg 'sheet'

                            or

                        an existing DataFrame with footprint data
            clear - bool, if False, loaded Footprints will be appended to the
                    list of Footprints created from Tower and Conductor objects.
                    The default is False. If True, only the loaded Footprints
                    are used.
        kw:
            sheet - str, specifies the target sheet if an excel workbook with
                    multiple sheets is passed in to the footprints arg
                    (not needed for a csv or an excel book with 1 sheet)
        returns:
            res - a Results object"""
        #check properties
        if((not self.conductors) and (not self.towers)):
            raise(EMFError('Cannot calculate fields because there are no Tower or Conductor objects in the Model. All fields would be zero.'))
        #get components
        all_components = ['Bmax', 'Bres', 'Bx', 'By', 'Bz']
        if(components == 'all'):
            components = all_components
        else:
            for c in components:
                if(c not in all_components):
                    raise(EMFError("The only valid components are 'Bmax', 'Bres', 'Bx', 'By', and 'Bz'. The input component %s is not valid." % str(c)))
        #start time
        t_start = datetime.datetime.now()
        #get wire segments
        S = self.segments
        #get sample points
        x, y, z = self.x, np.sort(self.y), self.z
        #compute one set of phasors to get arrays started
        Ph_x, Ph_y, Ph_z = subcalc_calcs.B_field_grid(
                S[0][0], S[0][1], S[0][2], S[0][3], x, y, z)
        #updates
        print('segment calculations complete: 1/%d' % len(S)),
        #compute the fields of all other segments and add the phasors
        for i in range(1,len(S)):
            Ph = subcalc_calcs.B_field_grid(
                    S[i][0], S[i][1], S[i][2], S[i][3], x, y, z)
            Ph_x += Ph[0]
            Ph_y += Ph[1]
            Ph_z += Ph[2]
            print('\rsegment calculations complete: %d/%d' % (i+1, len(S))),
        #convert the results into a 2D grid
        Ph_x,Ph_y,Ph_z,X,Y = subcalc_calcs.grid_segment_results(Ph_x,Ph_y,Ph_z,x,y)
        #get the real components
        Bx,By,Bz,Bres,Bmax = subcalc_calcs.phasors_to_magnitudes(Ph_x,Ph_y,Ph_z)
        #create a Results object with info about the original calculations
        data = dict(X=X, Y=Y)
        for c in components:
            exec("data['%s'] = %s" % (c, c))
        info = {'Created on': str(datetime.datetime.now()),
                'Z Height (ft)': z,
                'Minimum X Coordinate': np.min(x),
                'Minimum Y Coordinate': np.min(y),
                'Maximum X Coordinate': np.max(x),
                'Maximum Y Coordinate': np.max(y),
                'X Grid Increment': self.spacing,
                'Y Grid Increment': self.spacing,
                'Number X Grid Points': len(x),
                'Number Y Grid Points': len(y),
                'Total Grid Points': len(x)*len(y),
                'Distance Units': 'feet',
                'B-Field Units': 'mG'}
        res = Results(data, info)
        #load Footprints
        if(footprints is not None):
            if(clear):
                res.load_footprints(footprints, True, **kw)
            else:
                res._footprints = self.footprints
                res.load_footprints(footprints, False, **kw)
        else:
            res._footprints = self.footprints
        #copy the name over if it has been set
        if(self._name):
            res.name = self.name
        #print elapsed time
        t_end = datetime.datetime.now()
        print('\ntotal calculation time: %g seconds' % (t_end - t_start).total_seconds())
        #return
        return(res)

    def sample(self, *args, **kw):
        """Parse the Tower and Conductor objects into individual wire segements and compute the magnetic fields produced by the collection of wire segments at an arbitrary set of points in 3d space, defined by the x, y, and z inputs.
        args:
            x - iterable of numbers or a single number, the x coordinate(s)
                of the sample point(s)
            y - iterable of numbers or a single number, the y coordinate(s)
                of the sample point(s)
            z - iterable of numbers or a single number, the z coordinate(s)
                of the sample point(s)

                note: if any of the inputs are iterables, they must have
                      the same length as other iterable inputs.
        kw:
            n - integer, if provided, the input coordinates will be treated
                as a path and new points will be generated along it with
                linear interpolation to achieve a total number of points
                close to "n"
        returns:
            df - pandas DataFrame containing all components of the computed
                 fields, the sampled coordinates, and the cumulative distance
                 along the sampling path. The samples form a MultiIndex in
                 the returned DataFrame, which allows the data frame to be
                 indexed by sampled points
                    e.g. df[x,y,z]"""

        #start time
        t_start = datetime.datetime.now()
        #deal with lengths and make sure to work with np arrays
        args = [subcalc_funks._check_to_array(i) for i in args]
        L = max([len(i) for i in args])
        if(L != 1):
            for i in range(len(args)):
                if(len(args[i]) == 1):
                    args[i] = np.array([float(args[i][0])]*L)
                elif(len(args[i]) != L):
                    raise(EMFError("Cannot parse inputs if multiple inputs are iterables (vectors) and they have different lenghts. Inputs must be numbers or iterables with the same length as other iterable inputs."))
            #check for interpolation
            if('n' in kw):
                n = kw['n']
                args = subcalc_funks._interp_path_3D(args[0], args[1], args[2], n)
        x, y, z = args[0], args[1], args[2]
        #get wire segments
        S = self.segments
        #compute one set of phasors to get arrays started
        Ph_x, Ph_y, Ph_z = subcalc_calcs.B_field_general(
                S[0][0], S[0][1], S[0][2], S[0][3], x, y, z)
        #updates
        print('segment calculations complete: 1/%d' % len(S)),
        #compute the fields of all other segments and add the phasors
        for i in range(1,len(S)):
            Ph = subcalc_calcs.B_field_general(
                    S[i][0], S[i][1], S[i][2], S[i][3], x, y, z)
            Ph_x += Ph[0]
            Ph_y += Ph[1]
            Ph_z += Ph[2]
            print('\rsegment calculations complete: %d/%d' % (i+1, len(S))),
        #get the real components
        Bx,By,Bz,Bres,Bmax = subcalc_calcs.phasors_to_magnitudes(Ph_x,Ph_y,Ph_z)
        #create a dataframe
        dist = subcalc_funks.cumulative_distance(x, y, z)
        df = pd.DataFrame(
                data=dict(Bmax=Bmax, Bres=Bres, Bx=Bx, By=By, Bz=Bz, dist=dist),
                index=pd.MultiIndex.from_arrays([x,y,z], names=['x', 'y', 'z']))
        #print elapsed time
        t_end = datetime.datetime.now()
        print('\ntotal calculation time: %g seconds' % (t_end - t_start).total_seconds())
        #return
        return(df)

    def _xy_ranges(self):
        'Find the range of x and y coordinaes occupied by Tower and Conductor objects in the Model'
        x_ma = -np.inf
        x_mi = np.inf
        y_ma = -np.inf
        y_mi = np.inf
        x = subcalc_funks._flatten([t.conductor_x for t in self.towers]
                                    + [c.x for c in self.conductors])
        y = subcalc_funks._flatten([t.conductor_y for t in self.towers]
                                    + [c.y for c in self.conductors])
        for i in range(len(x)):
            if(x[i] > x_ma):
                x_ma = x[i]
            if(x[i] < x_mi):
                x_mi = x[i]
            if(y[i] > y_ma):
                y_ma = y[i]
            if(y[i] < y_mi):
                y_mi = y[i]
        if(x_mi != np.inf):
            xr = (x_mi, x_ma)
            yr = (y_mi, y_ma)
        else:
            xr, yr = None, None
        return(xr, yr)

class Tower(object):
    """Tower objects are used to define the locations of power lines in a Model.

    Towers are connected by their 'group' strings. Towers with the same 'group' are assumed to carry the same wires through a model domain. Thus, because Tower objects in the same group are carrying the same conductors, the order of their 'h', 'v', 'I', and 'phase' properties matters. For example, for a single Tower, all of the 0th elements of those properties correspond to a the same individual wire. For the path of that wire to be accurately calculated, it must also be in the 0th position in Towers of the same group. Towers in the same group must also have the same number of wires, meaning their 'h', 'v', 'I', and 'phase' arrays must have the same lengths.

    The sequence of Towers in a group is defined by thier 'seq' property. Current is assumed to flow from the Tower with Tower with the lowest 'seq' value to the Tower with the next lowest 'seq' value, and so on."""

    def __init__(self, group, seq, tower_x, tower_y, rot, h, v, I, phase):
        """
        args:
            group - str, a string used to group towers together into continuous
                    circuits
            seq - int, a number describing the position of this tower among
                  other towers in the same group, where current is assumed to
                  run from the 0th tower to the 1st, 2nd, ...
            tower_x - the x coordinate of the tower/pole in the model grid (ft)
            tower_y - the y coordinate of the tower/pole in the model grid (ft)
            tower_rot - the rotation angle of the tower in degrees, where
                        zero degrees points along the positive x axis and the
                        rotation increases clockwise
            h - iterable of numbers, the horizontal locations of conductors
                on the tower (ft)
            v - iterable of numbers, the vertical locations of conductors
                on the tower (ft)
            I - iterable of numbers, the current amplitudes of conductors
                on the tower (Amps)
            phase - iterable of numbers, the phase angles of conductors on the
                    tower (degrees)"""

        self._group = None
        self._seq = None
        self._tower_x = None
        self._tower_y = None
        self._rot = None
        self._h = None
        self._v = None
        self._I = None
        self._phase = None

        self.group = group
        self.seq = seq
        self.tower_x = tower_x
        self.tower_y = tower_y
        self.rot = rot
        self.h = h
        self.v = v
        self.I = I
        self.phase = phase

    #-----------------------------------------------------------------------
    #properties

    def _get_group(self): return(self._group)
    def _set_group(self, value): self._group = str(value)
    group = property(_get_group, _set_group, None, 'String used to group towers together into continuous circuits. All Towers in the same circuit have the same group string.')

    def _get_seq(self): return(self._seq)
    def _set_seq(self, value): self._seq = int(value)
    seq = property(_get_seq, _set_seq, None, 'Integer describing the position of this tower among other towers in the same group, where current is assumed to run from the 0th tower to the 1st, 2nd, ...')

    def _get_tower_x(self): return(self._tower_x)
    def _set_tower_x(self, value): self._tower_x = float(value)
    tower_x = property(_get_tower_x, _set_tower_x, None, 'Tower x coordinate within the model domain (ft)')

    def _get_tower_y(self): return(self._tower_y)
    def _set_tower_y(self, value): self._tower_y = float(value)
    tower_y = property(_get_tower_y, _set_tower_y, None, 'Tower y coordinate within the model domain (ft)')

    def _get_conductor_x(self): return(self.tower_x + np.cos(self.rot_rad)*self.h)
    conductor_x = property(_get_conductor_x, None, None, 'x-axis coordinates of the conductors on the tower in the model domain, calculated from Tower.x, Tower.h, and Tower.rot (ft)')

    def _get_conductor_y(self): return(self.tower_y + np.sin(self.rot_rad)*self.h)
    conductor_y = property(_get_conductor_y, None, None, 'y-axis coordinates of the conductors on the tower in the model domain, calculated from Tower.y, Tower.h, and Tower.rot (ft)')

    def _get_conductor_z(self): return(self.v)
    conductor_z = property(_get_conductor_z, None, None, 'z-axis coordinates of the conductors on the tower in the model domain, identical to Tower.h (ft)')

    def _get_rot(self): return(self._rot)
    def _set_rot(self, value): self._rot = float(value)
    rot = property(_get_rot, _set_rot, None, 'Tower rotation (degrees), where zero degrees points along the positive x axis and the rotation increases clockwise')

    def _get_rot_rad(self): return((self._rot/360.0)*2*np.pi)
    def _set_rot_rad(self, value): self._rot = (float(value)/(2*np.pi))*360.0
    rot_rad = property(_get_rot_rad, _set_rot_rad, None, 'Tower rotation (radians), where zero degrees points along the positive x axis and the rotation increases clockwise')

    def _get_h(self): return(self._h)
    def _set_h(self, value): self._h = np.array(value, dtype=float)
    h = property(_get_h, _set_h, None, 'Horizontal locations of conductors on the tower (ft)')

    def _get_v(self): return(self._v)
    def _set_v(self, value): self._v = np.array(value, dtype=float)
    v = property(_get_v, _set_v, None, 'Vertical locations of conductors on the tower (ft)')

    def _get_I(self): return(self._I)
    def _set_I(self, value): self._I = np.array(value, dtype=float)
    I = property(_get_I, _set_I, None, 'Current amplitudes of conductors on the tower (Amps)')

    def _get_phase(self): return(self._phase)
    def _set_phase(self, value): self._phase = np.array(value, dtype=float)
    phase = property(_get_phase, _set_phase, None, 'Phase angles of conductors on the tower (degrees)')

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(subcalc_print._str_Tower(self))

    def __len__(self):
        if(not (len(self.h) == len(self.v) == len(self.I) == len(self.phase))):
            raise(EMFError("""The 'h', 'v', 'I', and 'phase' arrays of a Tower object must have the same length. Tower with group '%s' and seq '%s' has nonuniform array lengths for those properties.""" % (self.group, str(self.seq))))
        return(len(self.h))


class Conductor(object):
    """Conductor objects represent single wires running through a Model domain. They can be used to specify the coordinates of conductor paths directly, instead of through the properties of a Tower object. Tower and Conductor objects together are the means of creating wires in a model domain."""

    def __init__(self, name, x, y, z, I, phase):
        """
        args:
            name - string, name of the conductor (i.e. 'Wire 1')
            x - iterable of at least two numbers, the x coordinates of the
                Conductor path in the model domain (ft)
            y - iterable of at least two numbers, the y coordinates of the
                Conductor path in the model domain (ft)
            z - iterable of at least two numbers, the z coordinates of the
                Conductor path in the model domain (ft)
            I - float, amplitude of the current in the Conductor, assumed to
                flow from (x[0], y[0], z[0]) to (x[-1], y[-1], z[-1]) (Amps)
            phase - float, the phase angle of the current flow (degrees)"""

        self._name = None
        self._x = None
        self._y = None
        self._z = None
        self._I = None
        self._phase = None

        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.I = I
        self.phase = phase

    #---------------------------------------------------------------------------
    #properties

    def _get_name(self): return(self._name)
    def _set_name(self, value): self._name = str(value)
    name = property(_get_name, _set_name, None, 'Name of the Conductor')

    def _get_x(self): return(self._x)
    def _set_x(self, value):
        if(not hasattr(value, '__len__')):
            raise(EMFError('The x property of a Conductor must be an iterable of at least two numbers.'))
        elif(len(value) < 2):
            raise(EMFError('The x property of a Conductor must have at least two numeric elements.'))
        self._x = np.array(value, dtype=float)
    x = property(_get_x, _set_x, None, 'x coordinates of the Conductor path in the model domain (ft)')

    def _get_y(self): return(self._y)
    def _set_y(self, value):
        if(not hasattr(value, '__len__')):
            raise(EMFError('The y property of a Conductor must be an iterable of at least two numbers.'))
        elif(len(value) < 2):
            raise(EMFError('The y property of a Conductor must have at least two numeric elements.'))
        self._y = np.array(value, dtype=float)
    y = property(_get_y, _set_y, None, 'y coordinates of the Conductor path in the model domain (ft)')

    def _get_z(self): return(self._z)
    def _set_z(self, value):
        if(not hasattr(value, '__len__')):
            raise(EMFError('The z property of a Conductor must be an iterable of at least two numbers.'))
        elif(len(value) < 2):
            raise(EMFError('The z property of a Conductor must have at least two numeric elements.'))
        self._z = np.array(value, dtype=float)
    z = property(_get_z, _set_z, None, 'z coordinates of the Conductor path in the model domain (ft)')

    def _get_I(self): return(self._I)
    def _set_I(self, value):
        if(not subcalc_funks._is_number(value)):
            raise(EMFError('The I property of Conductor objects must be a number defining the amplitude of the current in the Conductor (Amps).'))
        self._I = float(value)
    I = property(_get_I, _set_I, None, 'amplitude of the current in the Conductor, assumed to flow from (x[0], y[0], z[0]) to (x[-1], y[-1], z[-1]) (Amps)')

    def _get_phase(self): return(self._phase)
    def _set_phase(self, value):
        if(not subcalc_funks._is_number(value)):
            raise(EMFError('The phase property of Conductor objects must be a number defining the phase angle of the current (degrees).'))
        self._phase = float(value)
    phase = property(_get_phase, _set_phase, None, 'the phase angle of the current flow (degrees)')

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(subcalc_print._str_Conductor(self))

    def __len__(self):
        if(not (len(self.x) == len(self.y) == len(self.z))):
            raise(EMFError("Conductor arrays in the 'x', 'y', and 'z' properties must have the same lengths."))
        return(len(self.x))

class Results(object):
    """Results objects store the results of magnetic field modeling generated by the SubCalc program. A Results can be created by passing grid results directly or by passing a dictionary of grid results. The function subcalc.load_results is the best way to generate a Results object from a .REF file containing SubCalc results. The Results object can be saved to a more flexible (and smaller) excel file with Results.export(). Then the excel file can be read back into a Results object using subcalc.load_results.

    Results objects have a 'Bkey' property that determines which component of the magnetic field results is accessed by the Results.B property. For example, when Results.Bkey == 'Bmax' all operations involving the grid of magnetic field results accessed by Results.B will operate on the 'max' field component. When Bmax == 'Bz' all operations deal with the vertical 'z' component of the field, and so on. This includes plotting.

    Footprints of objects in the Results domain like buildings, the power lines, fences, etc. can also be stored in Results objects. Data for these objects can be saved in csv template files. The path of the footprint csv files can be passed to subcalc.load_results() for automatic inclusion in a newly generated Results object or it can be passed to an existing Results with Results.load_footprints(). The footprint data is stored in Footprint objects that have very little functionality and are mostly just organizational objects.

    Several methods are available for interpolating new values from the grid of field results: Results.interp(), Results.segment(), Results.path(), and Results.resample().

    There are also methods for selecting a subset of a Results object's domain and for shifting its x,y coordinates. These are Results.zoom() and Results.rereference() respectively.

    Contour plots of the fields (Results.B) can be automatically generated with subcalc.plot_contour() and colormesh plots can be automatically generated with subcalc.plot_pcolormesh(), both of which accept Results objects. The fields along a path through the Results domain (essentially a cross section) can be plotted with subcalc.plot_path. A contour or pcolormesh can be combined with cross sections using the subcalc.plot_cross_sections() function."""

    def __init__(self, *args, **kw):
        """Grid data must be passed in
        args:
            either:
                X - 2D array of x coordinates
                Y - 2D array of y coordinats
                B - 2D array of magnetic field magnitudes
                info - dict, optional, dictionary of results metadata
            or:
                data - dict, dictionary with X, Y, and B grids, should be
                            keyed by 'X','Y','Bmax','Bres','Bx','By','Bz'
                info - dict, optional, dictionary of results metadata
        kw:
            name - str, name for the object, used for filenames and such
            Bkey - str, selects which component of field results to use
                    or defines which component is passed as grid arrays
                    ('Bmax', 'Bres', 'Bx', 'By', 'Bz')."""

        largs = len(args)

        self._name = None
        if('name' in kw):
            self.name = kw['name']

        #inputs must correspond to the second case, dictionary of data and
        #optional dict of metadata
        if(largs <= 2):
            #check args[0]
            if(type(args[0]) is not dict):
                raise(EMFError("""The first argument to Results() must be a dictionary of result information when passing 1 or 2 arguments to initialize the Results, not %s""" % type(args[0])))
            #check keys
            s = set(['X', 'Y', 'Bmax', 'Bres', 'Bx', 'By', 'Bz'])
            k = set(args[0].keys())
            if(any([(i not in s) for i in k])):
                raise(EMFError("""If passing a dictionary to initialize a Results object, the dict can have the following keys only:
                    %s""" % str(s)))
            if(('X' not in k) or ('Y' not in k)):
                raise(EMFError("""If passing a dictionary to initialize a Results object, the dictionary must have 'X' and 'Y' keys, which lead to grid coordinate arrays."""))
            if(len(k) < 3):
                raise(EMFError("""If passing a dictionary to initialize a Results object, the dictionary must be keyed by 'X', 'Y', and any number of the following keys:
                    %s""" % str(['Bmax','Bres','Bx','By','Bz'])))
            #store data
            self._grid = args[0]
            #deal with Bkey
            if('Bkey' in kw):
                self.Bkey = kw['Bkey']
            elif('Bmax' in k):
                self.Bkey = 'Bmax'
            else:
                k = list(k)
                k.remove('X')
                k.remove('Y')
                self.Bkey = k[0]
            #store the info dict if present
            if(largs == 2):
                if(type(args[1]) is not dict):
                    raise(EMFError("""The fourth argument to Results() must be a dictionary of  results information, not %s""" % type(args[1])))
                self._info = args[1]
            else:
                self._info = {}

        #inputs must correspond to the first case, grids passed in directly
        elif(largs <= 4):
            #check input types
            msg = """If passing three or four arguments to initialize a Results, the first three arguments must be 2D numpy arrays representing X, Y, and B grids respectively, each with the same shape."""
            for i in range(3):
                if(type(args[i]) is not np.ndarray):
                    raise(EMFError(msg))
                elif(len(args[i].shape) != 2):
                    raise(EMFError(msg))
                elif(args[i].shape != args[1].shape):
                    raise(EMFError(msg))
            #2D reference grid arrays
            if('Bkey' in kw):
                self._Bkey = kw['Bkey']
            else:
                self._Bkey = 'unknown'
            self._grid = {'X': args[0], 'Y': args[1], self._Bkey: args[2]}
            #store the info dict if present
            if(largs == 4):
                if(type(args[3]) is not dict):
                    raise(EMFError("""The fourth argument to Results() must be a dictionary of results information, not type %s""" % type(args[3])))
                self._info = args[3]
            else:
                self._info = None

        #other reference objects in the results domain, like substatio
        #boundaries, stored in a list of Footprint objects
        self._footprints = []
        #angle of the northern direction with respect to the grid
        #   where 0 degrees is the positive y axis and clockwise is increasing
        self._north_angle = None
        #private list of column names for parsing and creating footprint dfs
        self._footprint_cols = ['Name', 'Group', 'X', 'Y', 'Power Line?',
                                    'Of Concern?', 'Draw as Loop?']

    #---------------------------------------------------------------------------
    #properties

    def _get_name(self):
        if(self._name is None):
            return('unnamed-results')
        else:
            return(self._name)
    def _set_name(self, value):
        if(not value):
            raise(EMFError('Results.name should be a string.'))
        self._name = str(value)
    name = property(_get_name, _set_name, None, 'Name of results object, used for filenames when exporting')

    def _get_B(self):
        return(self._grid[self._Bkey])
    B = property(_get_B, None, None, '2D grid of magnetic field results with y-coordinates decreasing down the rows and x-coordinates inreasing along the columns. For example, B[0,0] retrieves the result at the lowest x value and highest y value and B[-1,-1] retrieves the result at the higest x value and lowest y value. The grid is arranged so that the coordinate arrays print in the same way they are arranged on a cartesian plane.')

    def _get_Bkeys(self):
        k = self._grid.keys()
        k.remove('X')
        k.remove('Y')
        return(set(k))
    Bkeys = property(_get_Bkeys, None, None, 'A set of available Bkey values in the Results')

    def _get_Bkey(self):
        return(self._Bkey)
    def _set_Bkey(self, value):
        if(value not in self.Bkeys):
            raise(EMFError('Bkey must be set to one of the following elements:\n    %s' % ', '.join(sorted([repr(i) for i in self.Bkeys]))))
        else:
            self._Bkey = value
    Bkey = property(_get_Bkey, _set_Bkey, None, 'Component of magnetic field accessed by the B property')

    def _get_info(self): return(self._info)
    info = property(_get_info, None, None, 'Dictionary of results metadata')

    def _get_footprints(self): return(self._footprints)
    footprints = property(_get_footprints, None, None, 'List of Footprint objects')

    def _get_X(self):
        return(self._grid['X'])
    X = property(_get_X, None, None, '2D grid of reference grid x coordinates')

    def _get_Y(self):
        return(self._grid['Y'])
    Y = property(_get_Y, None, None, '2D grid of reference grid y coordinates')

    def _get_x(self):
        return(self.X[0,:])
    x = property(_get_x, None, None, 'Unique x values in results grid (column positions)')

    def _get_y(self):
        return(self.Y[:,0])
    y = property(_get_y, None, None, 'Unique y values in results grid (row positions)')

    def _get_xmax(self):
        return(np.max(self.x))
    xmax = property(_get_xmax, None, None, 'Maximum horizontal coordinate in results')

    def _get_xmin(self):
        return(np.min(self.x))
    xmin = property(_get_xmin, None, None, 'Minimum horizontal coordinate in results')

    def _get_xlim(self): return(self.xmin, self.xmax)
    xlim = property(_get_xlim, None, None, '(xmin, xmax)')

    def _get_ymax(self):
        return(np.max(self.y))
    ymax = property(_get_ymax, None, None, 'Maximum vertical coordinate in results')

    def _get_ymin(self):
        return(np.min(self.y))
    ymin = property(_get_ymin, None, None, 'Minimum vertical coordinate in results')

    def _get_ylim(self): return(self.ymin, self.ymax)
    ylim = property(_get_ylim, None, None, '(ymin, ymax)')

    def _get_Bmax(self): return(subcalc_funks._2Dmax(self.B)[0])
    Bmax = property(_get_Bmax, None, None, 'Maximum value of Results.B')

    def _get_loc_Bmax(self):
        m,i,j = subcalc_funks._2Dmax(self.B)
        return(self.x[j], self.y[i])
    loc_Bmax = property(_get_loc_Bmax, None, None, 'x,y coordinates of maximum of Results.B')

    def _get_idx_Bmax(self):
        m,i,j = subcalc_funks._2Dmax(self.B)
        return(i, j)
    idx_Bmax = property(_get_idx_Bmax, None, None, 'indices of maximum of Results.B, along the 0th axis then the 1st axis')

    def _get_Bmin(self): return(subcalc_funks._2Dmin(self.B)[0])
    Bmin = property(_get_Bmin, None, None, 'Minimum value of Results.B')

    def _get_spacing(self):
        if(np.any(np.diff(np.diff(self.x)) > (np.mean(self.x)*1e-6))):
            sx = 'non-uniform'
        else:
            sx = abs(self.x[1] - self.x[0])
        if(np.any(np.diff(np.diff(self.y)) > (np.mean(self.x)*1e-6))):
            sy = 'non-uniform'
        else:
            sy = abs(self.y[0] - self.y[1])
        if(sx == sy):
            return(sx)
        else:
            return(sx, sy)
    spacing = property(_get_spacing, None, None, 'Distance between grid points along the x and y axes (ft). If the spacing is different along the axes, a tuple of two floats is returned. Otherwise only a single float is returned. If the spacing along an axis is non-uniform (can be done with resampling), a string is returned.')

    def _get_N(self): return(self.B.shape[0]*self.B.shape[1])
    N = property(_get_N, None, None, 'The total number of sample points in the Results grid')

    def _get_north_angle(self):
        return(self._north_angle)
    def _set_north_angle(self, angle):
        if(not subcalc_funks._is_number(angle)):
            raise(EMFError("""The 'north_angle' attribute of a Results object must be a number."""))
        else:
            self._north_angle = float(angle)
    north_angle = property(_get_north_angle, _set_north_angle, None, """Angle of the Northern direction in degrees, where 0 represents the vertical or Y direction and clockwise represents increasing angle""")

    def _get_footprint_groups(self):
        """Generate a list of lists of Footprints with identical tags"""
        u = list(self.footprint_group_names)
        groups = [[] for i in range(len(u))]
        for i in range(len(self.footprints)):
            fp = self.footprints[i]
            groups[u.index(fp.group)].append(fp)
        return(groups)
    footprint_groups = property(_get_footprint_groups)

    def _get_footprint_group_names(self): return(set([fp.group for fp in self.footprints]))
    footprint_group_names = property(_get_footprint_group_names, None, None, 'A set of the unique footprint group names in the Results object')

    def _get_footprint_df(self):
        #create a df for each fp
        cols = self._footprint_cols
        dfs = []
        for fp in self.footprints:
            x, y = fp.x, fp.y
            r = range(len(x))
            dfs.append(pd.DataFrame({
                cols[0]: [fp.name for i in r],
                cols[1]: [fp.group for i in r],
                cols[2]: x,
                cols[3]: y,
                cols[4]: [int(fp.power_line) for i in r],
                cols[5]: [int(fp.of_concern) for i in r],
                cols[6]: [int(fp.draw_as_loop) for i in r]
            }))
        #concatenate the dfs and return
        df = pd.concat(dfs, axis=0, ignore_index=True)
        return(df)
    footprint_df = property(_get_footprint_df, None, None, 'A DataFrame of Footprint information in the same format as the emf.subcalc footprint template')

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(subcalc_print._str_Results(self))

    def load_footprints(self, footprints, clear=True, **kw):
        """Read footprint data from a csv/excel file and organize it in Footprint objects stored in self.footprints
        args:
            footprints - string, path to the footprint csv/excel data.
                        If footprint data is in an excel workbook with,
                        multiple sheets, the sheet name must be passed
                        to the keyword argumment 'sheet'

                            or

                        an existing DataFrame with footprint data
        optional args:
            clear - bool, True will replace all existing footprints with the
                    loaded ones, False will simply add loaded footprints to
                    the list of existing ones.
        kw:
            sheet - str, specifies the target sheet if an excel workbook with
                    multiple sheets is passed in"""
        #load file if footprints is not a DataFrame
        t = type(footprints)
        if(t is pd.DataFrame):
            df = footprints
        elif((t is str) or (t is unicode)):
            if('.' in footprints):
                ext = footprints[footprints.rfind('.')+1:]
                if(ext == 'xlsx'):
                    dfs = pd.read_excel(footprints, sheetname=None)
                    if(len(dfs.keys()) > 1):
                        if('sheet' in kw):
                            df = dfs[kw['sheet']]
                        else:
                            raise(EMFError("""If an excel file with multiple sheets is passed to Results.load_footprints, the target sheet must be specified with the keyword argument 'sheet'."""))
                    else:
                        df = dfs[dfs.keys()[0]]
                elif(ext == 'csv'):
                    df = pd.read_csv(footprints)
                else:
                    raise(EMFError("Only csv and xlsx files can be passed to Results.load_footprints."))
            else:
                raise(EMFError("""No extension was detected at the end of file name "%s" passed to Results.load_footprints. File names passed to load_footprints must have .csv or .xlsx extensions.""" % str(footprints)))
        else:
            raise(EMFError('Only pandas DataFrames and strings can be passed to Results.load_footprints, not type "%s"' % str(t)))

        #do a little data conditioning
        df = df.dropna(how='all').fillna(method='ffill')
        #match columns with string distance method
        cols = self._footprint_cols
        df.columns = subcalc_funks._Levenshtein_group(df.columns.values, cols)
        #store an error message
        msg = """
            The column:
                "%s"
            must contain only a single value for each unique footprint.
            The value should simply be repeated to fill all cells.
            The column contains multiple values for footprint name:
                "%s" """
        #pick out some columns
        fields = cols[4:]
        #clear the footprints list
        if(clear):
            self._footprints = []
        #create a footprint out of rows with the same "Name"
        for name, df in df.groupby(cols[0]):
            #check that certain fields only contain a single entry
            for f in fields:
                if(len(df[f].unique()) > 1):
                    raise(EMFError(msg % (f,n)))
            #create a footprint object
            row = df.iloc[0]
            fp = Footprint(row[cols[0]], row[cols[1]],
                    df[cols[2]].values, df[cols[3]].values,
                    row[cols[4]], row[cols[5]], row[cols[6]])
            #append the Footprint to the Results object's list
            self._footprints.append(fp)

    def concern_points(self, n=101):
        """Find the maximum fields at all footprint objects in the Results with
        their "of_concern" properties set to True.
        optional args:
            n - int, number of samples to interpolate along each footprint
        returns:
            x_concern - array, x coordinates of the concern points
            y_concern - array, y coordinates of the concern points
            B_concern - array, the magnetic field magnitudes at the points"""
        #get footprints of concern
        fps = [fp for fp in self.footprints if (fp.of_concern)]
        L = len(fps)
        #find the maximum field along each footprints path
        x_concern = np.empty((L,), dtype=float)
        y_concern = np.empty((L,), dtype=float)
        B_concern = np.empty((L,), dtype=float)
        for i in range(L):
            fp = fps[i]
            x, y, B = self.path(zip(fp.x, fp.y), n=n)
            idx = np.argmax(B)
            x_concern[i], y_concern[i], B_concern[i] = x[i], y[i], B[i]
        return(x_concern, y_concern, B_concern)

    def path(self, points, n=101, close_path=False):
        """Interpolate the field along a path defined by lists of x and y
        coordinates
        args:
            points - an iterable of x,y pairs representing a path through the
                     results grid, for example: [(1,2), (1,3), (2,4)]
        optional args:
            n - integer, approximate total number of points sampled. Each
                segment will have at least two samples (the beginning and end
                of the segment). The default is 101.
            close_path - bool, if True, append a segment to the end of the
                         path connecting the first and last points
        returns:
            x - array of x coordinates sampled along the path, will
                include all input coordinates
            y - array of y coordinates sampled along the path, will
                include all input coordinates
            B_interp - interpolated values corresponding to the input
                       coordinates"""
        #check closing path kw
        if(close_path):
            points = list(points)
            points.append(points[-1])
        #get sample points
        x, y = subcalc_funks._interp_path_2D(x, y, n)
        #interpolate
        B_interp = self.interp(x, y)

        return(np.array(x), np.array(y), np.array(B_interp))

    def segment(self, p1, p2, n=101):
        """Interpolate the field along a line between two points
        args:
            p1 - iterable, an x-y pair
            p2 - iterable, an x-y pair
        optional args:
            n - integer, number of points sampled (default 101)
        returns:
            x - array, x coordinates of interpolated values
            y - array, y coordinates of interpolated values
            B_interp - array, interpolated field values"""
        #check point lengths
        if((len(p1) != 2) or (len(p2) != 2)):
            raise(EMFError('Points must consist of two values (xy pairs).'))
        #create x and y vectors
        x, y = np.linspace(p1[0], p2[0], n), np.linspace(p1[1], p2[1], n)
        B_interp = self.interp(x, y)
        return(x, y, B_interp)

    def interp(self, x, y):
        """Interpolate in the x and y directions to find an estimated B value at an x,y location within the results grid
        args:
            x - iterable or scalar, x coordinate(s) to interpolate at
            y - iterable or scalar, y coordinate(s) to interpolate at
        returns:
            B_interp - array or float, the interpolated field value"""
        #make x,y iterable if scalars are passed in
        if(not (hasattr(x, '__iter__') and hasattr(y, '__iter__'))):
            scalar = True
            x = np.array([x], dtype=float)
            y = np.array([y], dtype=float)
        else:
            scalar = False
        #check that all points are in the grid
        if(not all([self.in_grid(x[i],y[i]) for i in range(len(x))])):
            raise(EMFError("""
            x,y coordinates must fall inside the reference grid:
                range of x coordinates: %g to %g
                range of y coordinates: %g to %g""" %
                (self.xmin, self.xmax, self.ymin, self.ymax)))
        #interpolate
        B_interp = np.array([subcalc_funks._bilinear_interp(self, x[i], y[i])
                                for i in range(len(x))])
        #return
        if(scalar):
            return(B_interp[0])
        else:
            return(B_interp)

    def in_grid(self, x, y):
        """Check if an x,y coordinate pair is inside the Results grid
        args:
            x - float, x coordinate
            y - float, y coordinate
        returns:
            b - bool, True if x,y is in the grid, False if it's not"""
        if((x > self.xmax) or
            (x < self.xmin) or
                (y > self.ymax) or
                    (y < self.ymin)):
            return(False)
        else:
            return(True)

    def resample(self, **kw):
        """Resample the results grid along a new number of x,y values, a new selection of x,y values, or a new number of total values
        kw:
            x - int or iterable, new number of x samples or new selection of
                x samples
            y - int or iterable, new number of y samples or new selection of
                y samples
            N - int, new approximate number of total samples
                    - overrides x and y kw
                    - preserves approx ratio of number of x and y values
                    - rounds up to nearest possible whole number of points
            inplace - resample in place or return a new object, default is False
        returns:
            res_resample - a new Results object containing the resampled grid"""

        #check inputs
        if(('x' not in kw) and ('y' not in kw) and ('N' not in kw)):
            raise(EMFError('Keyword arguments x, y, or N are required by Results.resample'))
        if('inplace' in kw):
            inplace = bool(kw['inplace'])
        else:
            inplace = False
        #in place or not
        if(inplace):
            res = self
        else:
            res = copy.deepcopy(self)
        #store grid x,y extents
        xmin, xmax = res.xmin, res.xmax
        ymin, ymax = res.ymin, res.ymax
        #get 1D vectors of x and y coordinates to resample at, from kw
        if('N' in kw):
            N = kw['N']
            if(not subcalc_funks._is_int(N)):
                raise(EMFError('Keyword argument "N" must be a whole number'))
            N = float(N)
            aspect = float(ymax - ymin)/(xmax - xmin)
            N_y = np.ceil(np.sqrt(N*aspect))
            N_x = np.ceil(N/N_y)
            x = np.linspace(xmin, xmax, N_x)
            y = np.linspace(ymin, ymax, N_y)
        else:
            if('x' in kw):
                x = kw['x']
                if(subcalc_funks._is_int(x)):
                    x = np.linspace(xmin, xmax, x)
                else:
                    x = np.sort(np.array(x, dtype=float))
            else:
                x = res.x
            if('y' in kw):
                y = kw['y']
                if(subcalc_funks._is_int(y)):
                    y = np.linspace(ymin, ymax, y)
                else:
                    y = np.sort(np.array(y, dtype=float))
            else:
                y = res.y
        #flip y coordinates so that Y prints intuitively
        y = y[::-1]
        X, Y = np.meshgrid(x, y)
        #resample each grid
        for k in res._grid:
            if((k != 'X') and (k != 'Y')):
                #some arrays are flipped to conform to _interpn conventions
                res._grid[k] = _interpn((res.y[::-1], res.x), res._grid[k],
                                            (Y[::-1,:], X))
        #spatial grids
        res._grid['X'] = X
        res._grid['Y'] = Y

        #return
        if(not inplace):
            return(res)

    def rereference(self, x_ref=0, y_ref=0, inplace=False):
        """Redefine the coordinates of the bottom left corner of the results grid (the lowest values along each axis) and increment values in the spatial grids accordingly. If no value is provided, spatial grids are adjusted start at (0, 0).
        optional args:
            x_ref - float, new starting value for the x axis, default is zero
            y_ref - float, new starting value for the y axis, default is zero
            inplace - bool, if True, the Results is rereferenced in place,
                      otherwise a rereferenced copy is returned."""
        #get variables
        if(inplace):
            res = self
        else:
            res = self.copy()
        X, Y = res.X, res.Y
        #rereference the grid and child footprint objects
        xdif, ydif = np.min(X) - x_ref, np.min(Y) - y_ref
        X -= xdif
        Y -= ydif
        for i in range(len(res.footprints)):
            x, y = res._footprints[i]._x, res._footprints[i]._y
            res._footprints[i]._x = [(j - xdif) for j in x]
            res._footprints[i]._y = [(j - ydif) for j in y]
        if(not inplace):
            return(res)

    def zoom(self, x_range, y_range, inplace=False, rereference=False):
        """Select a sub-area of the Results grid
        args:
            x_range - iterable of floats, two values representing the range of
                      x values to include in the zoomed grid, if input is
                      implicitly False all values are included (no x zooming)
            y_range - iterable of floats, two values representing the range of
                      y values to include in the zoomed grid, if input is
                      implicitly False all values are included (no y zooming)
        optional args:
            inplace - bool, if True the Results grid is changed in place,
                      otherwise a new Results is returned with the zoomed grid
            rereference - bool or number, if a number, the spatial grids will
                          be rereferenced to begin at that value. If True, the
                          spatial grids are rereferenced to begin at zero."""
        #check the ranges
        if(x_range is False):
            x_range = (self.xmin, self.xmax)
        if(y_range is False):
            y_range = (self.ymin, self.ymax)
        x_range, y_range = sorted(x_range), sorted(y_range)
        p1, p2 = (x_range[0], y_range[0]), (x_range[1], y_range[1])
        if((not self.in_grid(p1[0], p1[1]) or (not self.in_grid(p2[0], p2[1])))):
            raise(EMFError('Cannot zoom outside the results grid limits. The x limits are [%g, %g] and the y limits are [%g, %g].' % (self.xmin, self.xmax, self.ymin, self.ymax)))
        #get variables
        if(inplace):
            res = self
        else:
            res = self.copy()
        x, y = res.x, res.y
        #zoom
        bx = (x <= max(x_range)) & (x >= min(x_range))
        by = (y <= max(y_range)) & (y >= min(y_range))
        for k in res._grid:
            res._grid[k] = res._grid[k][by,:][:,bx]
        #check rereferencing
        if((type(rereference) is int) or (type(rereference) is float)):
            res.rereference(rereference, rereference, inplace=True)
        elif(rereference):
            res.rereference(inplace=True)
        if(not inplace):
            return(res)

    def flatten(self):
        """Create 1 dimensional versions of the gridded X, Y, and B arrays
        returns:
            x - 1D numpy array with x coordinates
            y - 1D numpy array with y coordinates
            b - 1D numpy array with magnetic field values"""
        nrows, ncols = self.B.shape
        n = nrows*ncols
        x = np.reshape(self.X, n)
        y = np.reshape(self.Y, n)
        b = np.reshape(self.B, n)
        return(x, y, b)

    def export(self, **kw):
        """Export the grid data and accompanying info to an excel file with tabs for each Bfield component, another for the info dict, and a final one for footprints if they're present
        kw:
            Bkeys - iterable of strings, selects which components to export.
            path - string, output destination/filename for workbook"""
        #get appropriate export filename
        fn = subcalc_funks._path_manage(self.name, '.xlsx', **kw)
        #get components to export
        if('Bkeys' in kw):
            Bkeys = set(kw['Bkeys'])
            s = self.Bkeys
            for k in Bkeys:
                if(k not in s):
                    raise(EMFError('The only acceptable Bkeys are %s. The input "%s" is not valid.' % (_print_str_list(s), k)))
        else:
            Bkeys = self.Bkeys
        #create excel writing object
        xl = pd.ExcelWriter(fn, engine='xlsxwriter')
        #write grid data
        for k in self._grid:
            if((k != 'X') and (k != 'Y') and (k in Bkeys)):
                pd.DataFrame(self._grid[k], columns=self.x, index=self.y
                        ).to_excel(xl, sheet_name=k)
        #write metadata if present
        if(self.info):
            pd.DataFrame([self.info[k] for k in self.info], index=self.info.keys(),
                    columns=['Parameter Value']).sort_index().to_excel(
                            xl, sheet_name='info', index_label='Parameter Name')
        #write footprints if present
        if(self.footprints):
            self.footprint_df.to_excel(xl, sheet_name='footprints', index=False,
                    columns=self._footprint_cols)
        #save and print
        xl.save()
        print('results saved to: %s' % fn)

    def copy(self):
        'Return a deep copy of the Results object'
        return(copy.deepcopy(self))

class Footprint(object):

    def __init__(self, name, group, x, y, power_line, of_concern, draw_as_loop):
        """
        args:
            name - string, the name of the Footprint, i.e. "Substation"
            group - string, group strings that are identical between
                    Footprint objects designate them as part of the same
                    group for plotting, i.e. "Nearby Houses"
            x - iterable, x coordinates of Footprint
            y - iterable, y coordinates of Footprint
            power_line - bool, True indicates the footprint corresponds to
                         a modeled power line or circuit
            of_concern - bool, True if the Footprint represents an area
                         that is potentially concerned about EMF (homes)
            draw_as_loop - bool, True if the footprint should be plotted
                           as a closed loop"""
        #check x and y are the same length
        if(len(x) != len(y)):
            raise(EMFError("""
            Footprints must have the same number of x and y values to form
            spatial coordinates"""))

        #init variables
        self._name = None
        self._group = None
        self._x = None
        self._y = None
        self._power_line = None
        self._of_concern = None
        self._draw_as_loop = None

        #set attributes
        self.name = name
        self.group = group
        self.x = x
        self.y = y
        self.power_line = power_line
        self.of_concern = of_concern
        self.draw_as_loop = draw_as_loop

    def _get_name(self): return(self._name)
    def _set_name(self, value):
        if(not value):
            raise(EMFError("""Cannot set Footprint 'name' property to an implicitly False object. Use a string"""))
        self._name = str(value)
    name = property(_get_name, _set_name, None, 'Footprint name')

    def _get_group(self):
        if(pd.isnull(self._group) or (self._group is None)):
            return('')
        else:
            return(self._group)
    def _set_group(self, value):
        if(not value):
            raise(EMFError("""Cannot set Footprint 'group' property to an implicitly False object. Use a string."""))
        self._group = str(value)
    group = property(_get_group, _set_group, None, 'Footprint group')

    def _get_x(self):
        if(self.draw_as_loop):
            return(self._x + [self._x[0]])
        else:
            return(list(self._x))
    def _set_x(self, value):
        self._x = [float(i) for i in value]
    x = property(_get_x, _set_x, None, 'x coordinates of Footprint vertices')

    def _get_y(self):
        if(self.draw_as_loop):
            return(self._y + [self._y[0]])
        else:
            return(list(self._y))
    def _set_y(self, value):
        self._y = [float(i) for i in value]
    y = property(_get_y, _set_y, None, 'y coordinates of Footprint vertices')

    def _get_power_line(self): return(self._power_line)
    def _set_power_line(self, value): self._power_line = bool(value)
    power_line = property(_get_power_line, _set_power_line, None, 'bool, flag indicating whether the Footprint represents the path of a power line or circuit (1 = yes)')

    def _get_of_concern(self): return(self._of_concern)
    def _set_of_concern(self, value): self._of_concern = bool(value)
    of_concern = property(_get_of_concern, _set_of_concern, None, 'bool, flag indicating whether the Footprint is "of concern," which means that the maximum fields along its border will be annoted on plots')

    def _get_draw_as_loop(self): return(self._draw_as_loop)
    def _set_draw_as_loop(self, value): self._draw_as_loop = bool(value)
    draw_as_loop = property(_get_draw_as_loop, _set_draw_as_loop, None, 'bool, flag indicating whether the Footprint should be drawn as a closed loop')

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(subcalc_print._str_Footprint(self))
