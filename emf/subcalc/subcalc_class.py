from .. import np, pd, os, copy, _interpn

from ..emf_class import EMFError

import subcalc_funks
import subcalc_calcs

class Model(object):
    """Model objects store Tower and/or Conductor objects, information about the desired model grid, and provide the means of computing model results and returning them in a Results object."""

    def __init__(self, name='unnamed-model', **kw):
        """
        optional args:
            name - string, a name identifying the Model, used for filenames
        kw:
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

        self.name = name
        if('towers' in kw):
            for t in kw['towers']:
                self.add_tower(t)
        if('conductors' in kw):
            for c in kw['conductors']:
                self.add_conductor(c)
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

    def _get_name(self): return(self._name)
    def _set_name(self, value): self._name = str(value)
    name = property(_get_name, _set_name, None, 'A string identifying the Model, used for filenames')

    def _get_towers(self): return(self._towers)
    towers = property(_get_towers, None, None, 'List of Tower objects in the Model')

    def _get_tower_groups(self):
        """Generate a list of lists of Tower objects with identical group
        strings, with the sublists sorted according to the Tower seq property"""
        u = list(set([t.group for t in self.towers]))
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

    def _get_conductors(self): return(self._conductors)
    conductors = property(_get_conductors, None, None, 'List of Conductor objects in the Model')

    def _get_xmin(self):
        if(self._xmin is None):
            return(0.0)
        else:
            return(self._xmin)
    def _set_xmin(self, value):
        if(value >= self.xmax):
            raise(EMFError('xmin must be less than xmax.'))
        self._xmin = float(value)
    xmin = property(_get_xmin, _set_xmin, None, 'The minimum x value in the Model area in feet, defualts to zero')

    def _get_xmax(self):
        if(self._xmax is None):
            if((not self.conductors) and (not self.towers)):
                return(None)
            else:
                ma = -np.inf
                mi = np.inf
                for t in self.towers:
                    o = np.max(t.conductor_x)
                    if(o > ma):
                        ma = o
                    o = np.min(t.conductor_x)
                    if(o < mi):
                        mi = o
                for c in self.conductors:
                    o = np.max(c.x)
                    if(o > ma):
                        ma = o
                    o = np.min(c.x)
                    if(o < mi):
                        mi = o

                return(ma + abs(mi - self.xmin))
    def _set_xmax(self, value):
        if(value <= self.xmin):
            raise(EMFError('xmax must be greater than xmin'))
        self._xmax = float(value)
    xmax = property(_get_xmax, _set_xmax, None, 'The maximum x value in the Model area in feet, will use Tower and Conductor information to pick value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_ymin(self):
        if(self._ymin is None):
            return(0.0)
        else:
            return(self._ymin)
    def _set_ymin(self, value):
        if(value >= self.ymax):
            raise(EMFError('ymin must be less than ymax.'))
        self._ymin = float(value)
    ymin = property(_get_ymin, _set_ymin, None, 'The minimum y value in the Model area in feet, defualts to zero')

    def _get_ymax(self):
        if(self._ymax is None):
            if((not self.conductors) and (not self.towers)):
                return(None)
            else:
                ma = -np.inf
                mi = np.inf
                for t in self.towers:
                    o = np.max(t.conductor_y)
                    if(o > ma):
                        ma = o
                    o = np.min(t.conductor_y)
                    if(o < mi):
                        mi = o
                for c in self.conductors:
                    o = np.max(c.y)
                    if(o > ma):
                        ma = o
                    o = np.min(c.y)
                    if(o < mi):
                        mi = o

                return(ma + abs(mi - self.ymin))
        else:
            return(self._ymax)
    def _set_ymax(self, value):
        if(value <= self.ymin):
            raise(EMFError('ymax must be greater than ymin.'))
        self._ymax = float(value)
    ymax = property(_get_ymax, _set_ymax, None, 'The maximum y value in the Model area in feet, will use Tower and Conductor information to pick value automatically if left unset and return None if left unset with no Towers or Conductors in the Model')

    def _get_z(self): return(self._z)
    def _set_z(self, value):
        self._z = float(z)
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
            #check lengths
            if(any([(len(g[i]) != len(g[i+1])) for i in range(len(g) - 1)])):
                raise(EMFError('All Towers with the same group must have the same length. Group "%s" contains Towers with unequal lengths.' % g[0].group))
            #check currents and phases
            for i in range(len(g[0])):
                if(any([(g[j].I[i] != g[j+1].I[i]) for j in range(len(g) - 1)])):
                    raise(EMFError("""Wires with different currents were found between Towers in group "%s". Check the currents.""" % g[0].group))
                if(any([(g[j].phase[i] != g[j+1].phase[i]) for j in range(len(g) - 1)])):
                    raise(EMFError("""Wires with different phases were found between Towers in group "%s". Check the phasing.""" % g[0].group))
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
    segments = property(_get_segments, None, None, '')

    #---------------------------------------------------------------------------
    #methods

    def add_tower(self, t):
        """Add a Tower object to the Model
        args:
            t - Tower object to add (a copy is added)"""
        if(type(t) is not Tower):
            raise(EMFError('The add_tower method of Model objects can only add Tower objects.'))
        self._towers.append(copy.deepcopy(t))

    def add_conductor(self, c):
        """Add a Conductor object to the Model
        args:
            c - Conductor object to add (a copy is added)"""
        if(type(c) is not Conductor):
            raise(EMFError('The add_conductor method of Model objects can only add Conductor objects.'))
        self._conductors.append(copy.deepcopy(c))

    def calculate(self, *args):
        """Parse the Tower and Conductor objects into individual wire segements and compute the magnetic fields produced by the collection of wire segments in the model domain
        optional args:
            footprints - string, path to the footprint csv/excel data.
                        If footprint data is in an excel workbook with,
                        multiple sheets, the sheet name must be passed
                        to the kwarg 'sheet'

                            or

                        an existing DataFrame with footprint data
        returns:
            res - a Results object"""
        #check properties
        if((not self.conductors) and (not self.towers)):
            raise(EMFError('Cannot calculate fields because there are no Tower or Conductor objects in the Model. All fields would be zero.'))
        scc = subcalc_calcs
        #get wire segments
        S = self.segments
        #get sample points
        x, y, z = self.x, np.sort(self.y), self.z
        #compute one set of phasors to get arrays started
        Ph_x, Ph_y, Ph_z = scc.B_field_segment(S[0][0], S[0][1], S[0][2], S[0][3], x, y, z)
        #compute the fields of all other segments and add the phasors
        for i in range(1,len(S)):
            Ph = scc.B_field_segment(S[i][0], S[i][1], S[i][2], S[i][3], x, y, z)
            Ph_x += Ph[0]
            Ph_y += Ph[1]
            Ph_z += Ph[2]
        #convert the results into a 2D grid
        Ph_x, Ph_y, Ph_z, X, Y = scc.grid_segment_results(Ph_x, Ph_y, Ph_z, x, y)
        #get the real components
        Bx, By, Bz, Bres, Bmax = scc.phasors_to_magnitudes(Ph_x, Ph_y, Ph_z)
        #return a Results object
        data = dict(X=X, Y=Y, Bx=Bx, By=By, Bz=Bz, Bres=Bres, Bmax=Bmax)
        info = {'Z Height (ft)': z,
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
        if(len(args) > 0):
            res.load_footprints(args[0])
        return(res)

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

    def __len__(self):
        if(not (len(self.h) == len(self.v) == len(self.I) == len(self.phase))):
            raise(EMFError("""The 'h', 'v', 'I', and 'phase' arrays of a Tower object must have the same length. Tower with group '%s' and seq '%s' has nonuniform array lengths for those properties.""" % (self.group, str(self.seq))))
        return(len(self.h))


class Conductor(object):
    """Conductor objects represent single wires running through a Model domain. They can be used to specify the coordinates of conductor paths directly, instead of through the properties of a Tower object. Tower and Conductor objects together are the means of creating wires in a model domain."""

    def __init__(self, x, y, z, I, phase):
        """
        args:
            x - iterable of at least two numbers, the x coordinates of the
                Conductor path in the model domain (ft)
            y - iterable of at least two numbers, the y coordinates of the
                Conductor path in the model domain (ft)
            z - iterable of at least two numbers, the z coordinates of the
                Conductor path in the model domain (ft)
            I - float, amplitude of the current in the Conductor, assumed to
                flow from (x[0], y[0], z[0]) to (x[-1], y[-1], z[-1]) (Amps)
            phase - float, the phase angle of the current flow (degrees)"""

        self._x = None
        self._y = None
        self._z = None
        self._I = None
        self._phase = None

        self.x = x
        self.y = y
        self.z = z
        self.I = I
        self.phase = phase

    #---------------------------------------------------------------------------
    #properties

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

    def __len__(self):
        if(not (len(self.x) == len(self.y) == len(self.z))):
            raise(EMFError("Conductor arrays in the 'x', 'y', and 'z' properties must have the same lengths."))
        return(len(self.x))

class Results(object):
    """Results objects store the results of magnetic field modeling generated by the SubCalc program. A Results can be created by passing grid results directly or by passing a dictionary of grid results. The function subcalc.load_results is the best way to generate a Results object from a .REF file containing SubCalc results. The Results object can be saved to a more flexible (and smaller) excel file with Results.export(). Then the excel file can be read back into a Results object using subcalc.load_results.

    Results objects have a 'Bkey' property that determines which component of the magnetic field results is accessed by the Results.B property. For example, when Results.Bkey == 'Bmax' all operations involving the grid of magnetic field results accessed by Results.B will operate on the 'max' field component. When Bmax == 'Bz' all operations deal with the vertical 'z' component of the field, and so on. This includes plotting.

    Footprints of objects in the Results domain like buildings, the power lines, fences, etc. can also be stored in Results objects. Data for these objects can be saved in csv template files. The path of the footprint csv files can be passed to subcalc.load_results for automatic inclusion in a newly generated Results object or it can be passed to an existing Results with Results.load_footprints. The footprint data is stored in Footprint objects that have very little functionality and are mostly just organizational objects.

    Several methods are available for interpolating new values from the grid of field results: Results.interp, Results.segment, Results.path, and Results.resample.

    There are also methods for selecting a subset of a Results object's domain and for shifting its x,y coordinates. These are Results.zoom and Results.rereference respectively.

    Contour plots of the results (Results.B) can be automatically generated with subcalc.plot_contour(Results) and colormesh plots can be automatically generated with subcalc.plot_pcolormesh(Results). The fields along a path through the Results domain (essentially a cross section) can be plotted with subcalc.plot_path. A contour or pcolormesh can be combined with cross sections using the subcalc.plot_cross_sections function."""

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
            s = set(['X','Y','Bmax','Bres','Bx','By','Bz'])
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
                self.Bkey = list(k)[0]
            #store the info dict if present
            if(largs == 2):
                if(type(args[1]) is not dict):
                    raise(EMFError("""The fourth argument to Results() must be a dictionary of  results information, not %s""" % type(args[1])))
                self._info = args[1]
            else:
                self._info = None

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
        self.footprint_df = None #DataFrame of Footprint information
        self._footprints = []
        #angle of the northern direction with respect to the grid
        #   where 0 degrees is the positive y axis and clockwise is increasing
        self._north_angle = None

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
            raise(EMFError('Bkey must be set to one of the following elements:\n%s' % str(self.Bkeys)))
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

    def _get_ymax(self):
        return(np.max(self.y))
    ymax = property(_get_ymax, None, None, 'Maximum vertical coordinate in results')

    def _get_ymin(self):
        return(np.min(self.y))
    ymin = property(_get_ymin, None, None, 'Minimum vertical coordinate in results')

    def _get_Bmax(self): return(subcalc_funks._2Dmax(self.B)[0])
    Bmax = property(_get_Bmax, None, None, 'Maximum value of Results.B')

    def _get_loc_Bmax(self):
        m,i,j = subcalc_funks._2Dmax(self.B)
        return(self.x[j], self.y[j])
    loc_Bmax = property(_get_loc_Bmax, None, None, 'x,y coordinates of maximum of Results.B')

    def _get_idx_Bmax(self):
        m,i,j = subcalc_funks._2Dmax(self.B)
        return(i, j)
    idx_Bmax = property(_get_idx_Bmax, None, None, 'indices of maximum of Results.B, along the 0th axis then the 1st axis')

    def _get_Bmin(self): return(subcalc_funks._2Dmin(self.B)[0])
    Bmin = property(_get_Bmin, None, None, 'Minimum value of Results.B')

    def _get_spacing(self):
        sx = abs(self.x[1] - self.x[0])
        sy = abs(self.y[1] - self.y[0])
        if(sx != sy):
            raise(EMFError('Sample spacing is different in the x and y dimensions.'))
        return(sx)
    spacing = property(_get_spacing, None, None, 'Distance between grid points along the x and y axes (ft)')

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
        u = list(set([f.group for f in self.footprints]))
        groups = [[] for i in range(len(u))]
        for i in range(len(self.footprints)):
            fp = self.footprints[i]
            groups[u.index(fp.group)].append(fp)
        return(groups)
    footprint_groups = property(_get_footprint_groups)

    #---------------------------------------------------------------------------
    #methods

    def __str__(self):
        return(
        'Results object\n    name: %s\n    components/Bkeys: %s\n    B field range (%s): %g to %g mG\n    x limits: [%g, %g] ft\n    y limits: [%g, %g] ft\n    total samples: %d' %
        (repr(self.name), str(self.Bkeys), self.Bkey, self.Bmin, self.Bmax, self.xmin, self.xmax, self.ymin, self.ymax, len(self.x)*len(self.y))
        )

    def load_footprints(self, footprints, **kw):
        """Read footprint data from a csv/excel file and organize it in Footprint objects stored in self.footprints
        args:
            footprints - string, path to the footprint csv/excel data.
                        If footprint data is in an excel workbook with,
                        multiple sheets, the sheet name must be passed
                        to the keyword argumment 'sheet'

                            or

                        an existing DataFrame with footprint data
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
        df.fillna(method='ffill', inplace=True)
        #store the DataFrame
        self.footprint_df = df
        #match columns with string distance method
        cols = ['Group', 'Name', 'X', 'Y', 'Power Line?',
                'Of Concern?', 'Draw as Loop?']
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
        self._footprints = []
        #create a footprint out of rows with the same "Name"
        for name, df in df.groupby('Name'):
            #check that certain fields only contain a single entry
            for f in fields:
                if(len(df[f].unique()) > 1):
                    raise(EMFError(msg % (f,n)))
            #create a footprint object
            row = df.iloc[0]
            fp = Footprint(name, row['Group'], df['X'].values, df['Y'].values,
                    row['Power Line?'], row['Of Concern?'], row['Draw as Loop?'])
            #append the Footprint to the Results object's list
            self._footprints.append(fp)

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
        L = len(points) - 1
        rL = range(L)
        #distribute point count for multiple segements, or don't
        if(L > 1):
            #get x and y
            x, y = zip(*points)
            #calculate distances of each segment
            d = np.array([np.sqrt((x[i]-x[i+1])**2+(y[i]-y[i+1])**2) for i in rL])
            #convert distances to fractions
            d = d/sum(d)
            #approximately distribute the total number of points to each segment
            n = np.ceil(n*d)
            #make sure there are at least two samples in each segment
            n[n < 2] = 2
        else:
            n = [n]
        #interpolate over each segment
        segs = [self.segment(points[i], points[i+1], n[i]) for i in rL]
        x, y, B_interp = (subcalc_funks._flatten(i) for i in zip(*segs))

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
        if(not (hasattr(x, '__len__') and hasattr(y, '__len__'))):
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
        returns:
            res_resample - a new Results object containing the resampled grid"""

        #check inputs
        if(not kw):
            raise(EMFError('Keyword arguments are required by Results.resample'))
        #store grid x,y extents
        xmin, xmax = self.xmin, self.xmax
        ymin, ymax = self.ymin, self.ymax
        #get 1D vectors of x and y coordinates to resample at, from kw
        if('N' in kw):
            N = kw['N']
            if(not subcalc_funks._is_int(N)):
                raise(EMFError('Keyword argument "N" must be a whole number'))
            N = float(N)
            aspect = float(self.B.shape[0])/self.B.shape[1]
            N_y = np.ceil(np.sqrt(N/aspect))
            N_x = np.ceil(N/N_y)
            x = np.linspace(xmin, xmax, N_x)
            y = np.linspace(ymin, ymax, N_y)
        else:
            if('x' in kw):
                x = kw['x']
                if(subcalc_funks._is_int(x)):
                    x = np.linspace(xmin, xmax, x)
            else:
                x = self.x
            if('y' in kw):
                y = kw['y']
                if(subcalc_funks._is_int(y)):
                    y = np.linspace(ymin, ymax, y)
            else:
                y = self.y
        #flip y coordinates so that Y prints intuitively
        y = y[::-1]
        #resample the grid
        X, Y = np.meshgrid(x, y)
        #some arrays have to be flipped to conform to conventions of _interpn
        B_resample = _interpn((self.y[::-1], self.x),
                                self.B[::-1,:],
                                (Y[::-1,:], X))

        #return with re-flipped results
        res_resample = Results(X, Y, B_resample[::-1,:],
                copy.deepcopy(self.info),
                Bkey=self.Bkey)
        res_resample.load_footprints(self.footprint_df)
        return(res_resample)

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
            res._footprints[i]._x = [j - xdif for j in x]
            res._footprints[i]._y = [j - ydif for j in y]
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
            path - string, output destination/filename for workbook"""
        #get appropriate export filename
        fn = subcalc_funks._path_manage(self.name, '.xlsx', **kw)
        #create excel writing object
        xl = pd.ExcelWriter(fn, engine='xlsxwriter')
        #write grid data
        for k in self._grid:
            if((k != 'X') and (k != 'Y')):
                pd.DataFrame(self._grid[k], columns=self.x, index=self.y
                        ).to_excel(xl, sheet_name=k)
        #write results information if present
        if(self.info is not None):
            pd.DataFrame([self.info[k] for k in self.info], index=self.info.keys(),
                    columns=['Parameter Value']).sort_index().to_excel(
                            xl, sheet_name='info', index_label='Parameter Name')
        #write footprint DataFrame if present
        if(self.footprint_df is not None):
            self.footprint_df.to_excel(xl, sheet_name='footprints', index=False)
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

        self._res = None #parent Results object

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

    def __str__(self):
        """quick and dirty printing"""
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            s += str(k) + ': ' + str(v[k]) + '\n'
        return(s)
