from .. import np, pd, copy

from ..emf_class import EMFError

from . import fields_funks
from . import fields_plots
from . import fields_calcs
from . import fields_print
from . import FIELDS_io

class Conductor(object):
    """A single power line object representing an infintely long conductor and storing basic physical parameters like the line's 2D coordinates, size, voltage, current, phase angle, etc. The Conductor class is the functional unit of emf.fields. Groups of Conductors are organized by higher level CrossSection objects to form complete 2D EMF models of parellel sets of power lines. CrossSection objects use the physical parameters stored in Conductor objects to compute predicted EMF values at a fixed height across a preset distance perpendicular to the power lines."""

    def __init__(self, name, params=None, copy_cond=None):
        """
        args:
            name - scalar hashable, mandatory, essentially just a name
        optional args:
            params - dict or list, available Conductor parameters are:

                    x, y, subconds, d_cond, d_bund, V, I, phase

                A list can be passed with the values of any number of the above
                paremeters (in the order listed above) or a dict can be passed
                with any number of those parameters as keys.

            copy_cond - an existing Conductor object, parameters that are not
                        set with the second arg will be borrowed from the
                        Conductor"""

        self._name = None #conductor label (string)
        self._xs = None #parent CrossSection object
        self._freq = 60.0 #phase frequency (Hz)
        self._x = None #x coordinate (ft)
        self._y = None #y coordinate (ft)
        self._subconds = 1 #number of subconductors per bundle (integer)
        self._d_cond = None #conductor diameter (inches)
        self._d_bund = None #bundle diameter (inches)
        self._V = None #line voltage (kV)
        self._I = None #line current (Amps)
        self._phase = None #phase angle (degrees)

        #set properties property
        self.name = name
        #parameters
        if(params is not None):

            a = params
            p = self.params

            if(type(a) is list):
                for i in range(len(a)):
                    exec('self.%s = %s' % (p[i], repr(str(a[i]))))
                p = p[i+1:]

            elif(type(a) is dict):
                for k in a:
                    if(k not in p):
                        raise(EMFError("""Unexpected dictionary key "%s" encountered in Conductor initialization.""" % str(k)))
                    exec('self.%s = %s' % (k, repr(str(a[k]))))
                    p.remove(k)

            else:
                raise(EMFError("""The second argument to Conductor initialization must be a list or dict. See the Conductor class doc string for more info on the second arg."""))

        #conductor to copy from
        if(copy_cond is not None):

            if(type(copy_cond) is not type(self)):
                raise(EMFError("""The third argument must be an existing Conductor object from which all unset parameters in the newly created Conductor are copied."""))
            for q in p:
                if(eval('copy_cond.%s' % q) is not None):
                    exec('self.%s = copy_cond.%s' % (q, q))

    #---------------------------------------------------------------------------
    #PROPERTIES

    def _get_params(self): return(['x', 'y', 'subconds', 'd_cond', 'd_bund', 'V', 'I', 'phase'])
    params = property(_get_params, None, None, 'Get a list of the parameter/property names for Conductor objects')

    def _check_complete(self):
        """Check if all Conductor variables have been set
        returns:
            b - bool, True if all physical parameters are not None
            s - sentence indicating which parameter is empty"""
        for p in self.params:
            val = eval('self.%s' % p)
            if((val is None) or pd.isnull(val)):
                if(self._xs is not None):
                    return(False, 'The "%s" parameter in Conductor "%s" in CrossSection "%s" is unset.' % (p, self.name, self._xs.sheet))
                else:
                    return(False, 'The "%s" parameter in Conductor "%s" is unset.' % (p, self.name))
        return(True, None)
    complete = property(_check_complete)

    def _reset_xs_fields(self, old_value, new_value):
        """if the Conductor has a parent CrossSection, set the parent object's fields property to None"""
        if((self._xs is not None) and (old_value != new_value)):
            self._xs._fields = None

    def _check_to_float(self, value, prop):
        """check that an incoming value can be converted to a float and return the float version of it or raise an error"""
        if((not fields_funks._is_number(value)) or pd.isnull(value)):
            raise(EMFError("""Conductor property '%s' must have a non-null numeric value. It cannot be set to: %s""" % (prop, repr(value))))
        return(float(value))

    def _check_to_int(self, x, prop):
        """check that an incoming value can be converted to an integer and return the int version of it or raise an error"""
        if(not fields_funks._is_number(x)):
            raise(EMFError("""Conductor property '%s' must be numeric. It cannot be set to: %s""" % (prop, repr(x))))
        elif(int(x) != float(x)):
            raise(EMFError("""Conductor property '%s' must be an integer. It cannot be set to: %s""" % (prop, repr(x))))
        return(int(x))

    def _get_name(self): return(self._name)
    def _set_name(self, new_name):
        if(not new_name):
            raise(EMFError("""Conductor names cannot be implicitly False objects."""))
        if(self._xs is not None):
            if(new_name in self._xs.names):
                raise(EMFError("""Cannot change Conductor name from "%s" to "%s" because the latter is alread present in the Conductor's parent CrossSection object. CrossSections must contain Conductors with unique names.""" % (str(self._name), str(new_name))))
        old_name = self._name
        self._name = new_name
        if((new_name != old_name) and (self._xs is not None)):
            self._xs._update_name2idx()
    name = property(_get_name, _set_name, None, """Conductor object identification, usually a string but can be anything hashable""")

    def _get_freq(self): return(self._freq)
    freq = property(_get_freq, None, None, """Conductor frequency (Hz), cannot be set because Conductors in the same CrossSection are assumed to have equal frequencies, making the frequency irrelevant. This is a stand in value that does not currently affect EMF calculations in any way.""")

    def _get_x(self): return(self._x)
    def _set_x(self, new_value):
        old_value = self._x
        self._x = self._check_to_float(new_value, 'x')
        self._reset_xs_fields(old_value, new_value)
    x = property(_get_x, _set_x, None, """Conductor x coordinate (ft)""")

    def _get_y(self): return(self._y)
    def _set_y(self, new_value):
        old_value = self._y
        self._y = self._check_to_float(new_value, 'y')
        self._reset_xs_fields(old_value, new_value)
    y = property(_get_y, _set_y, None, """Conductor y coordinate (ft)""")

    def _get_subconds(self): return(self._subconds)
    def _set_subconds(self, new_value):
        old_value = self._subconds
        self._subconds = self._check_to_int(new_value, 'subconds')
        self._reset_xs_fields(old_value, new_value)
    subconds = property(_get_subconds, _set_subconds, None, """The number of subconductors in the conductor "bundle," defaults to 1 if left unset""")

    def _get_d_cond(self): return(self._d_cond)
    def _set_d_cond(self, new_value):
        old_value = self._d_cond
        self._d_cond = self._check_to_float(new_value, 'd_cond')
        self._reset_xs_fields(old_value, new_value)
    d_cond = property(_get_d_cond, _set_d_cond, None, """Diameter of the conductor or a single subconductor if subconds > 1 (inches)""")

    def _get_d_bund(self):
        if(self.subconds == 1):
            if(self._d_bund is None):
                return(self.d_cond)
        return(self._d_bund)
    def _set_d_bund(self, new_value):
        old_value = self._d_bund
        self._d_bund = self._check_to_float(new_value, 'd_bund')
        self._reset_xs_fields(old_value, new_value)
    d_bund = property(_get_d_bund, _set_d_bund, None, """Diameter of the conductor bundle, defaults to d_cond if subconds == 1 (inches)""")

    def _get_V(self): return(self._V)
    def _set_V(self, new_value):
        old_value = self._V
        self._V = self._check_to_float(new_value, 'V')
        self._reset_xs_fields(old_value, new_value)
    V = property(_get_V, _set_V, None, """Conductor phase voltage, not line voltage (kilovolts, kV)""")

    def _get_I(self): return(self._I)
    def _set_I(self, new_value):
        old_value = self._I
        self._I = self._check_to_float(new_value, 'I')
        self._reset_xs_fields(old_value, new_value)
    I = property(_get_I, _set_I, None, """Conductor phase current, not line current (Amps)""")

    def _get_phase(self): return(self._phase)
    def _set_phase(self, new_value):
        old_value = self._phase
        self._phase = self._check_to_float(new_value, 'phase')
        self._reset_xs_fields(old_value, new_value)
    phase = property(_get_phase, _set_phase, None, """Conductor phase angle (degrees)""")

    #---------------------------------------------------------------------------
    #METHODS

    def __str__(self):
        return(fields_print._str_Conductor(self))

    def __eq__(self, other):
        if(type(self) is not type(other)):
            return(False)
        for p in self.params:
            if(eval('self.%s' % p) != eval('other.%s' % p)):
                return(False)
        return(True)

    def __ne__(self, other):
        if(self == other):
            return(False)
        else:
            return(True)

    def copy(self):
        'Return a deep copy of the Conductor object'
        return(copy.deepcopy(self))

class CrossSection(object):
    """An object representing a single model of a set of parellel power lines. CrossSections are containers for Conductor objects, providing dict-like access to the Conductors by their 'name' properties, which also store modeling information like the sampling locations and the horizontal locations of right-of-way (ROW) edges, which usually mark the edge of the utility company's property and the location of regulatory interest. CrossSections use the physical parameters stored in Conductor objects to compute predicted EMF values."""

    def __init__(self, sheet, conds=[], title='', group=None):
        """CrossSection must be initialized with a sheet string, and a list of Conductor objects can optionally be passed in in the second arg
        args:
            sheet - string, sets CrossSection.sheet (basically a name)
        optional args:
            conds - optional list of Conductor objects to copy into the xs
            title - string, sets CrossSection.title
            group - group name, sets CrossSection.group"""

        self._sheet = None #mandatory, short, template worksheet name
        self._group = None #identifier linking multiple CrossSection objects
        self._title = None #longer form, used for plotting text
        self._sb = None #parent SectionBook object
        self.soil_resistivity = 100.0 #?
        self._max_dist = None #maximum simulated distance from the ROW center
        self._step = 1.0 #step size for calculations
        self._sample_height = 3.0 #uniform sample height
        self._lROW = None #exact coordinate of the left ROW edge
        self._rROW = None #exact coordinate of the left ROW edge
        self._conds = [] #list of Conductor objects
        #dictionary mapping Conductor names to Conductor objects
        self._name2idx = dict()
        #integer indexer
        self._i = _IntegerIndexer(self._conds)
        #DataFrame storing results, populated with _calculate_fields()
        self._fields = None

        #set properties
        self.sheet = sheet
        for c in conds:
            self.add_conductor(c)
        self.title = title
        self.group = group

    #---------------------------------------------------------------------------
    #PROPERTIES

    def _get_params(self): return(['max_dist', 'step', 'lROW', 'rROW'])
    params = property(_get_params, None, None, 'Get a list of the parameter/property names which must be set for the object to be "complete".')

    def _check_complete(self):
        """Check that all variables in each conductor in the CrossSection are set, and not left with the initial None values. Also check that CrossSection properties needed for calculations (max_dist, step, lROW, rROW) are set.
        returns:
            b - True if all Conductors are complete, False if not
            s - sentence indicating which conductor is unset and which
                parameter in it is the cause"""
        if(not self.conds):
            return(False, 'There are no Conductors in CrossSection "%s".' % self.sheet)
        for p in self.params:
            if(eval('self.%s' % p) is None):
                return(False, 'The "%s" property of CrossSection "%s" is unset.' % (p, self.sheet))
        for c in self.conds:
            b, s = c.complete
            if(b is False):
                return(b, s)
        return(True, None)
    complete = property(_check_complete)

    def _check_to_float(self, value, prop):
        """check that an incoming value can be converted to a float and return the float version of it, or raise an error"""
        if(not fields_funks._is_number(value)):
            raise(EMFError("""CrossSection property '%s' must be numeric. It cannot be set to: %s""" % (prop, repr(value))))
        return(float(value))

    def _reset_fields(self, old_value, new_value):
        """set the fields property to None"""
        if(old_value != new_value):
            self._fields = None

    def _update_parent_sb_sheets(self, old_sheet, new_sheet):
        """if a parent SectionBook is present, update its indexing dictionary"""
        if(self._sb is not None):
            if(old_sheet != new_sheet):
                self._sb._update_sheet2idx()

    def _get_conds(self): return(self._conds)
    conds = property(_get_conds, None, None, """List of Conductor objects added to the CrossSection. Do not modify it directly, but use the CrossSection methods add_conductor() and remove_conductor() to manage which Conductor objects are in a CrossSection.""")

    def _get_i(self): return(self._i)
    i = property(_get_i, None, None, """Use CrossSection.i[] to retrieve Conductor objects by the integer order in which they were added to the CrossSection. For example, CrossSection.i[0] gets the first Conductor added, CrossSection.i[4] gets the fifth one added, etc.""")

    def _get_names(self): return([c.name for c in self.conds])
    names = property(_get_names, None, None, 'Return a list of the names of Conductors in the CrossSection, in the order in which the Conductors were added')

    def _get_sheet(self): return(self._sheet)
    def _set_sheet(self, new_value):
        if(not new_value):
            raise(EMFError('The CrossSection sheet property cannot be implicitly False.'))
        old_value = self._sheet
        self._sheet = new_value
        self._update_parent_sb_sheets(old_value, new_value)
    sheet = property(_get_sheet, _set_sheet, None, 'The principle identification variable of a CrossSection, should be a string briefly identifying the object. CrossSections are retrieved from SectionBook objects by their sheet strings.')

    def _get_group(self): return(self._group)
    def _set_group(self, value): self._group = value
    group = property(_get_group, _set_group, None, """Very short variable used to group CrossSections in SectionBooks. CrossSection objects in SectionBooks that have the same group properties are grouped for functions that involve comparisons.""")

    def _get_title(self):
        if(not self._title):
            return(self.sheet)
        else:
            return(self._title)
    def _set_title(self, value): self._title = value
    title = property(_get_title, _set_title, None, """Longer form identification string, used in plot titles and such. If unset, will fall back to CrossSection.sheet""")

    def _get_max_dist(self): return(self._max_dist)
    def _set_max_dist(self, new_value):
        #check input
        new_value = self._check_to_float(new_value, 'max_dist')
        if(new_value < 0):
            raise(EMFError('max_dist must be a positive number.'))
        #make sure max_dist will not cut out the ROW edge locations
        if((self.lROW is not None) and (abs(self.lROW) > new_value)):
            raise(EMFError('A max_dist value of %s would cut off the left ROW edge location of %s. The max_dist value must be greater than or equal to the magnitudes of each ROW edge location.' % (new_value, self.lROW)))
        if((self.rROW is not None) and (abs(self.rROW) > new_value)):
            raise(EMFError('A max_dist value of %s would cut off the right ROW edge location of %s. The max_dist value must be greater than or equal to the magnitudes of each ROW edge location.' % (new_value, self.rROW)))
        #set
        old_value = self._max_dist
        self._max_dist = new_value
        self._reset_fields(old_value, new_value)
    max_dist = property(_get_max_dist, _set_max_dist, None, """Maximum horizontal distance of sample points for EMF calculations (ft), default is 100. Resolution or step size is controlled by the 'step' property.""")

    def _get_step(self): return(self._step)
    def _set_step(self, new_value):
        old_value = self._step
        self._step = self._check_to_float(new_value, 'step')
        self._reset_fields(old_value, new_value)
    step = property(_get_step, _set_step, None, 'Horizontal distance between sample points for EMF calculations (ft), default is 1')

    def _get_sample_height(self): return(self._sample_height)
    def _set_sample_height(self, new_value):
        new_value = self._check_to_float(new_value, 'sample_height')
        if(new_value < 0):
            raise(EMFError('CrossSection sample_height must be greater than or equal to zero.'))
        old_value = self._sample_height
        self._sample_height = new_value
        self._reset_fields(old_value, new_value)
    sample_height = property(_get_sample_height, _set_sample_height, None, 'Height of all sample points (ft), default is 3')

    def _get_lROW(self): return(self._lROW)
    def _set_lROW(self, new_value):
        #check input
        new_value = self._check_to_float(new_value, 'lROW')
        if((self.rROW is not None) and (new_value >= self.rROW)):
            raise(EMFError('lROW must be less than rROW.'))
        if((self.max_dist is not None) and (abs(new_value) > self.max_dist)):
            raise(EMFError('The magnitude of lROW cannot exceed max_dist.'))
        #set
        old_value = self._lROW
        self._lROW = new_value
        self._reset_fields(old_value, new_value)
    lROW = property(_get_lROW, _set_lROW, None, """Horizontal location of the left (negative x) edge of the right-of-way (ROW), a point of interest on the left side of the model (ft)""")

    def _get_rROW(self): return(self._rROW)
    def _set_rROW(self, new_value):
        #check input
        new_value = self._check_to_float(new_value, 'rROW')
        if((self.lROW is not None) and (new_value <= self.lROW)):
            raise(EMFError('rROW must be greater than lROW.'))
        if((self.max_dist is not None) and (abs(new_value) > self.max_dist)):
            raise(EMFError('The magnitude of rROW cannot exceed max_dist.'))
        #set
        old_value = self._rROW
        self._rROW = new_value
        self._reset_fields(old_value, new_value)
    rROW = property(_get_rROW, _set_rROW, None, """Horizontal location of the right (positive x) edge of the right-of-way (ROW), a point of interest on the left side of the model (ft)""")

    def _get_hot(self): return([c for c in self.conds if c.V != 0])
    hot = property(_get_hot, None, None, """Return a list of the hot (nonzero voltage) Conductor objects in the CrossSection""")

    def _get_gnd(self): return([c for c in self.conds if c.V == 0])
    gnd = property(_get_gnd, None, None, """Return a list of the grounded (zero voltage) Conductor objects in the CrossSection""")

    def _get_x_sample(self):
        u = np.floor(self.max_dist/self.step)
        v = np.linspace(-self.step*u, self.step*u, 2*u + 1)
        if((self.lROW not in v) or (self.rROW not in v)):
            v = np.array(sorted(list(set(v) | {self.lROW, self.rROW})))
        return(v)
    x_sample = property(_get_x_sample, None, None, """Compute an array of horizontal coordinates (ft) for EMF calculations, based on 'max_dist' and 'step'. The left and right ROW edge locations ('lROW' and 'rROW') are always included in the returned array.""")

    def _get_y_sample(self):
        return(self.sample_height*np.ones((len(self.x_sample),), dtype=float))
    y_sample = property(_get_y_sample, None, None, """Return an array the same length as 'x_sample' with all values equal to 'sample_height'""")

    def _get_freq(self):
        return(np.array([c.freq for c in self.conds], dtype=float))
    freq = property(_get_freq, None, None, """Generate an array of Conductor frequencies, in the same order as 'conds'""")

    def _get_x(self): return(np.array([c.x for c in self.conds], dtype=float))
    x = property(_get_x, None, None, """Generate an array of Conductor x coordinates (ft), in the same order as 'conds'""")

    def _get_y(self): return(np.array([c.y for c in self.conds], dtype=float))
    y = property(_get_y, None, None, """Generate an array of Conductor y coordinates (ft), in the same order as 'conds'""")

    def _get_subconds(self):
        return(np.array([c.subconds for c in self.conds], dtype=float))
    subconds = property(_get_subconds, None, None, """Generate an array of Conductor subconductor counts (see Conductor.subconds), in the same order as 'conds'""")

    def _get_d_cond(self):
        return(np.array([c.d_cond for c in self.conds], dtype=float))
    d_cond = property(_get_d_cond, None, None, """Generate an array of Conductor diameters (inches, see Conductor.d_cond), in the same order as 'conds'""")

    def _get_d_bund(self):
        return(np.array([c.d_bund for c in self.conds], dtype=float))
    d_bund = property(_get_d_bund, None, None, """Generate an array of Conductor bundle diameters (inches, see Conductor.d_bund), in the same order as 'conds'""")

    def _get_V(self):
        return(np.array([c.V for c in self.conds], dtype=float))
    V = property(_get_V, None, None, """Generate an array of Conductor voltages (kilovolts, kV), in the same order as 'conds'""")

    def _get_I(self):
        return(np.array([c.I for c in self.conds], dtype=float))
    I = property(_get_I, None, None, """Generate an array of Conductor currents (Amps), in the same order as 'conds'""")

    def _get_phase(self):
        return(np.array([c.phase for c in self.conds], dtype=float))
    phase = property(_get_phase, None, None, """Generate an array of Conductor phase angles (degrees), in the same order as 'conds'""")

    def _get_fields(self):
        if(self._fields is None):
            self._calculate_fields()
        return(self._fields)
    fields = property(_get_fields, None, None, """A pandas DataFrame of magnetic and electric field results at the points defined by 'x_sample' and 'y_sample'. The results are calculated and stored upon the first reference to this property. If Conductor objects in the CrossSection or the CrossSection itself are modified in relevant ways, the results are cleared and will be recalculated when 'fields' is next accessed. All updating/refreshing is done automatically. The DataFrame of results is indexed by the values in 'x_sample' and has the following columns 'Ex', 'Ey', 'Eprod', 'Emax', 'Bx', 'By', 'Bprod', and 'Bmax.' The columns correspond to different components of the electric and magnetic fields (see emf.fields.phasors_to_magnitudes()). The electric fields are reported in units of kV/m and the magnetic fields are reported in units of mG (milliGauss).""")

    def _get_ROW_edge_fields(self):
        return(self.fields.loc[[self.lROW, self.rROW]])
    ROW_edge_fields = property(_get_ROW_edge_fields, None, None, """Slice the 'fields' DataFrame and return another DataFrame with only the results at the left and right right-of-way (ROW) edges, the locations of which are set in 'lROW' and 'rROW'.""")

    def _get_ROW_edge_max(self):
        return(self.ROW_edge_fields[['Bmax','Emax']])
    ROW_edge_max = property(_get_ROW_edge_max, None, None, """Slice the 'fields' DataFrame and return another DataFrame with only the maximum field results at the left and right right-of-way edges, the locations of which are set by the 'lROW' and 'rROW' property""")

    #---------------------------------------------------------------------------

    def __str__(self):
        return(fields_print._str_CrossSection(self))

    def __getitem__(self, key):
        #check in the hot dict
        try:
            idx = self._name2idx[key]
        except(KeyError):
            raise(KeyError('There are no Conductors with the name "%s" in CrossSection "%s"' % (key, self.sheet)))
        else:
            return(self.conds[idx])
        #return None if no xs is found
        return(None)

    def __iter__(self):
        for c in self.conds:
            yield(c)

    def __len__(self):
        'The length of a CrossSection is the number of Conductors in it.'
        return(len(self._conds))

    def add_conductor(self, cond):
        """Add a copy of a Conductor object to the CrossSection, which can be accessed by indexing the CrossSection like a dictionary, using the Conductor's 'name' attribute.
        args:
            cond - Conductor object to add"""
        #check if the name has already been used
        if(cond.name in self._name2idx):
            raise(EMFError("""A Conductor with name "%s" is already in CrossSection "%s". Conductors in a CrossSection must have unique names and another Conductor with the same name cannot be added.""" % (cond.name, self.sheet)))
        #check if the Conductor is complete
        b, s = cond.complete
        if(not b):
            raise(EMFError("""Cannot add Conductor "%s" to the CrossSection because it is not complete. %s""" % (cond.name, s)))
        #see if the conductor is grounded
        if(cond.V == 0):
            #probably shouldn't have current if grounded
            if(cond.I != 0):
                print("""Conductor with name "%s" in CrossSection "%s" is grounded (V = 0) but has nonzero current?""" % (cond.name, self.sheet))
        #check for other conductors in the same location
        for c in self.conds:
            if((cond.x == c.x) and (cond.y == c.y)):
                raise(EMFError("""Can't add conductor "%s" to CrossSection "%s" because conductor "%s" exists in the exact same location (x = %g, y = %g). Conductors must be in unique locations.""" % (cond.name, self.sheet, c.name, c.x, c.y)))
        #add to self.conds and indexing dict
        self._name2idx[cond.name] = len(self.conds)
        self.conds.append(copy.deepcopy(cond))
        #associate xs with the conductor
        self.conds[-1]._xs = self
        #reset the CrossSection's fields
        self._fields = None

    def remove_conductor(self, key):
        """Remove a Conductor object from the CrossSection
        args:
            key - Conductor name or Conductor object to remove"""
        #check input
        if(type(key) is Conductor):
            key = key.name
        #check the key
        try:
            idx = self._name2idx[key]
        except(KeyError):
            raise(KeyError('Conductor name "%s" does not exist in CrossSection "%s".' % (key, self.sheet)))
        else:
            #pop the conductor out of the conductor list
            self.conds.pop(idx)
            #reassemble the naming dictionary
            self._update_name2idx()
            #reset the CrossSection's fields
            self._fields = None

    def _update_name2idx(self):
        self._name2idx = dict(zip(self.names, range(len(self.conds))))

    def _calculate_fields(self):
        """Calculate electric and magnetic fields across the ROW and store the results in the self.fields DataFrame"""
        #check for completeness (all variables that are needed for calculations
        #have been set to allowable values)
        b, s = self.complete
        if(not b):
            raise(EMFError("Cannot calculate fields. %s" % s))
        #pull arrays with conductor and sampling information
        x_sample, y_sample = self.x_sample, self.y_sample
        x, y, I, V, phase = self.x, self.y, self.I, self.V, self.phase
        subconds, d_cond, d_bund = self.subconds, self.d_cond, self.d_bund
        #calculate magnetic field
        Bx, By = fields_calcs.B_field(x, y, I, phase, x_sample, y_sample)
        Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
        #calculate electric field
        Ex, Ey = fields_calcs.E_field(x, y, subconds, d_cond, d_bund, V, phase,
                                        x_sample, y_sample)
        Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
        #store the values
        df = pd.DataFrame({
                'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
                'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax},
                index=x_sample)
        df.index.name = 'Distance (ft)'
        self._fields = df

    def sample(self, *args):
        """Calculate electric and magnetic fields at an arbitrary set of x,y points in the CrossSection
        args:
            x_sample* - iterable of numbers or a single number, the x
                        coordinate(s) of the sample point(s)
            y_sample* - iterable of numbers or a single number, the y
                        coordinate(s) of the sample point(s)

                *if both args are iterables, they must be the same length"""
        #deal with lengths and make sure to work with np arrays
        args = list(args)
        if(len(args) == 1):
            noy = True
            args.append(self.sample_height)
        elif(len(args) != 0):
            raise(EMFError("The sample method of a CrossSection object accepts 1 or 2 arguments, not %d" % len(args)))
        args = [fields_funks._check_to_array(i) for i in args]
        L = max([len(i) for i in args])
        if(L != 1):
            for i in range(len(args)):
                if(len(args[i]) == 1):
                    args[i] = np.array([float(args[i][0])]*L)
                elif(len(args[i]) != L):
                    raise(EMFError("Cannot parse inputs if multiple inputs are iterables (vectors) and they have different lenghts. Inputs must be numbers or iterables with the same length as other iterable inputs."))
        #unpack sample points
        x_sample, y_sample = args
        #pull out other modeling inputs
        x, y, I, V, phase = self.x, self.y, self.I, self.V, self.phase
        subconds, d_cond, d_bund = self.subconds, self.d_cond, self.d_bund
        #calculate magnetic field
        Bx, By = fields_calcs.B_field(x, y, I, phase, x_sample, y_sample)
        Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
        #calculate electric field
        Ex, Ey = fields_calcs.E_field(x, y, subconds, d_cond, d_bund, V, phase,
                                        x_sample, y_sample)
        Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
        #store the values
        if(noy):
            idx = x_sample
        else:
            idx = pd.MultiIndex.from_arrays([x_sample,y_sample], names=['x','y'])
        df = pd.DataFrame(
                {'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
                'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax},
                index=idx)
        if(noy):
            df.index.name = 'Distance (ft)'
        return(df)

    def compare_DAT(self, DAT_path, plot=True, **kw):
        """Load a FIELDS output file (.DAT) to calculate absolute and percentage differences between it and the CrossSection object's results. The results in the DAT file must be sampled at the same x coordinates as those in the CrossSection. A panel of comparison results is returned. If the 'save' or 'path' keywords are used, the comparison results will be saved with plots demonstrating the comparisons.
        args:
            DAT_path - path of FIELDS results file
        optional args:
            plot - bool, if True, two figures are returned in the figs
                   return variable
        kw:
            save - bool, toggle whether output Panel and figures are saved,
                    also saves figures with error and comparison plots
            path - string, destination of saved files, will force save == True
            round - int, round the results in self.fields to a certain number
                    of digits in an attempt to exactly match the FIELDS
                    results, which are printed only to the thousandths digit
            truncate - bool, truncate results after the thousandths digit
        returns:
            pan - pandas Panel with DAT results, results of this code,
                the absolute error between, and the relative error between
            figs - two figures comparing electric and magnetic fields results,
                   only returned if 'plot' is True"""
        #load the .DAT file into a dataframe
        df = FIELDS_io.read_DAT(DAT_path)
        #prepare a dictionary to create a Panel
        if(('round' in kw) and ('truncate' in kw)):
            raise(FLDError("""Cannot both round and truncate for DAT comparison. Choose either rounding or truncation."""))
        elif('round' in kw):
            f = self.fields.round(kw['round'])
        elif('truncate' in kw):
            if(kw['truncate']):
                f = self.fields.copy(deep=True)
                for c in f.columns:
                    for i in f.index:
                        f[c].loc[i] = float('%.3f' % f[c].loc[i])
        else:
            f = self.fields
        #create a panel (the frame names matter in the _plot_comparison function
        #called below)
        frames = ['FIELDS-results', 'emf.fields-results',
                'absolute-difference', 'percent-difference']
        pan = pd.Panel(data={frames[0] : df,
                frames[1]: f,
                frames[2]: f - df,
                frames[3]: 100*(f - df)/f})
        #make plots of the absolute and percent error
        if(plot):
            figs = fields_plots._plot_comparison(self, pan, 'FIELDS', **kw)
        #write data and save figures if called for
        if('path' in kw):
            kw['save'] = True
        if('save' in kw):
            if(kw['save']):
                fn = fields_funks._path_manage(self.sheet + '-DAT-comparison',
                    '.xlsx', **kw)
                xl = pd.ExcelWriter(fn, engine='xlsxwriter')
                for f in frames:
                    pan[f].to_excel(xl, index_label='Distance (ft)', sheet_name=f)
                xl.save()
                print('DAT comparison book saved to: %s' % fn)
        #return the Panel
        if(plot):
            return(pan, figs)
        else:
            return(pan)

    def rand(self, *args):
        """Get a random Conductor  or a list of random Conductors from the CrossSection
        args:
            an input integer determines the number of random xss fetched,
            which might not be unique"""
        c = self.conds
        if(len(args) == 0):
            return(self.conds[np.random.randint(len(c))])
        else:
            r = np.random.randint(len(c), size=args[0])
            return([c[i] for i in r])

    def copy(self):
        """Return a deep copy of the CrossSection object"""
        return(copy.deepcopy(self))

    def export(self, **kw):
        """Write the 'fields' DataFrame to a csv
        kw:
            path - string, destination/filename for saved file"""
        #path management
        fn = fields_funks._path_manage(self.sheet + '-all_results', '.csv', **kw)
        #write results
        self.fields.to_csv(fn, index_label='Distance (ft)')
        print('Cross section fields written to: %s' % fn)

class SectionBook(object):
    """A top level class organizing a group of CrossSection objects. It provides dict-like access to child CrossSection objects by their 'sheet' properties. It also maintains a table of maximum field results at the ROW edges of each CrossSection added (the ROW_edge_max property), provides exporting methods, and can be passed to plotting functions to generate plots comparing groups of CrossSections. CrossSections in a SectionBook are grouped by their 'group' properties. Those with identical groups are automatically grouped."""

    def __init__(self, name, xss=[]):
        """SectionBooks must be initialized with a name string, with the option to pass a list of CrossSection objects in the second arg that will be copied into the new object
        args:
            name - string, name of SectionBook
        optional args:
            xss - list of CrossSection objects to copy into the SectionBook"""

        self._name = name #mandatory identification field
        self._xss = [] #list of cross section objects
        self._sheet2idx = dict() #for CrossSection retrieval
        self._i = _IntegerIndexer(self._xss)
        #add CrossSections if they're passed in
        if(xss):
            for xs in xss:
                self.add_section(xs)

    #---------------------------------------------------------------------------
    #PROPERTIES

    def _get_name(self): return(self._name)
    def _set_name(self, value): self._name = value
    name = property(_get_name, _set_name, None, """Mandatory string identifying the SectionBook object""")

    def _get_xss(self): return(self._xss)
    xss = property(_get_xss, None, None, """List of CrossSection objects added to the SectionBook. Should not be edited directly, but with SectionBook.add_section() and SectionBook.remove_section().""")

    def _get_i(self): return(self._i)
    i = property(_get_i, None, None, """Use SectionBook.i[] to retrieve CrossSection objects by the integer order in which they were added to the SectionBook. For example, SectionBook.i[0] gets the first CrossSection added, SectionBook.i[4] gets the fifth one added, etc.""")

    def _get_sheets(self): return([xs.sheet for xs in self.xss])
    sheets = property(_get_sheets, None, None, """A list of CrossSection sheet strings, in the same order as 'xss'""")

    def _get_unique_group_names(self):
        return(set([xs.group for xs in self.xss]))
    unique_group_names = property(_get_unique_group_names, None, None,
            'A set of unique CrossSection groups')

    def _get_titles(self): return([xs.title for xs in self.xss])
    titles = property(_get_titles, None, None, """A list of CrossSection titles, in the same order as 'xss'""")

    def _get_groups(self):
        u = list(self.unique_group_names) #unique CrossSection groups
        groups = [[] for i in range(len(u))]
        for i in range(len(self.xss)):
            groups[u.index(self.xss[i].group)].append(self.xss[i])
        return(groups)
    groups = property(_get_groups, None, None, 'A list of lists of grouped CrossSection objects')

    def _get_ROW_edge_max(self):
        #gather ROW edge results
        xss = self.xss
        r = range(len(xss))
        #construct DataFrame
        df = pd.DataFrame(data={
            'Bmaxl': [xss[i].fields.at[xss[i].lROW, 'Bmax'] for i in r],
            'Bmaxr': [xss[i].fields.at[xss[i].rROW, 'Bmax'] for i in r],
            'Emaxl': [xss[i].fields.at[xss[i].lROW, 'Emax'] for i in r],
            'Emaxr': [xss[i].fields.at[xss[i].rROW, 'Emax'] for i in r]},
            index=self.sheets)
        return(df)
    ROW_edge_max = property(_get_ROW_edge_max, None, None, """DataFrame with maximum field magnitudes at the right-of-way (ROW) edges of each CrossSection in the SectionBook. The DataFrame is indexed by the CrossSection sheet strings, and has columns 'Bmaxl', 'Bmaxr', 'Emaxl', 'Emaxr'.""")

    def _check_complete(self):
        """Check that all CrossSections in the SectionBook are complete
        returns:
            b - True if all xcs are complete, False if not
            sheet - sheet of first incomplete xs, if any, else None
            cname - name of first incomplete Conductor in first incomplete xs,
                    if any, else None
            v - variable name of first unset variable in first incomplete
                Conductor in first incomplete xs, if any, else None"""
        if(not self.xss):
            return(False, 'There are no CrossSections in SectionBook "%s".' % self.name)
        for xs in self:
            b, s = xs.complete
            if(b is False):
                return(b, s)
        return(True, None)
    complete = property(_check_complete)

    #---------------------------------------------------------------------------

    def __getitem__(self, key):
        """Index the SectionBook by CrossSection names"""
        try:
            idx = self._sheet2idx[key]
        except(KeyError):
            raise(KeyError('There are no CrossSections with the sheet "%s" in SectionBook "%s"' % (key, self.name)))
        else:
            return(self.xss[idx])

    def __iter__(self):
        """Iteration over all CrossSections in the SectionBook"""
        for xs in self.xss:
            yield(xs)

    def __len__(self):
        """Length of SectionBook is the number of CrossSections in it"""
        return(len(self.xss))

    def __str__(self):
        return(fields_print._str_SectionBook(self))

    def rand(self, *args):
        """Get a random CrossSection from the SectionBook or a list of random CrossSection objects
        args:
            an input integer determines the number of random CrossSections returned in a list, which may not be unique"""
        if(len(args) == 0):
            return(self.xss[np.random.randint(len(self))])
        else:
            r = np.random.randint(len(self), size=args[0])
            return([self.xss[i] for i in r])

    def add_section(self, xs, **kw):
        """Copy a CrossSection into the SectionBook. Use this method instead of directly editing the 'xss' list.
        args:
            xs - CrossSection object to add to SectionBook (a copy is added)
        kw:
            sheet - a new sheet name for the added CrossSection"""
        #prevent adding CrossSections with the same sheets
        if(xs.sheet in self._sheet2idx):
            raise(EMFError("""A CrossSection with sheet "%s" already exists in the SectionBook. Duplicate names cannot be added to the lookup dictionary (self._sheet2idx). Use a different name.""" % xs.sheet))
        #check completeness
        b, s = xs.complete
        if(not b):
            raise(EMFError('Cannot add CrossSection "%s" to SectionBook "%s". %s' % (xs.sheet, self.name, s)))
        #copy xs
        xs = copy.deepcopy(xs)
        #get sheet
        if('sheet' in kw):
            xs.sheet = kw['sheet']
        else:
            #add the CrossSection and update indexing dict
            self._sheet2idx[xs.sheet] = len(self.xss)
            self.xss.append(xs)
            #associate self with the copied CrossSection
            self.xss[-1]._sb = self

    def remove_section(self, sheet):
        """Remove a CrossSection from the SectionBook by providing its sheet string
        args:
            sheet - sheet string of CrossSection to remove"""
        #check sheet is in self
        if(sheet not in self._sheet2idx):
            raise(EMFError("""Sheet %s cannot be removed because it does not exist in the SectionBook object""" % sheet))
        #delete from self.xss list
        idx = self._sheet2idx[sheet]
        self.xss.pop(idx)
        #update
        self._update_sheet2idx()

    def _update_sheet2idx(self):
        self._sheet2idx = dict(zip(self.sheets, range(len(self.xss))))

    def export(self, **kw):
        """Write complete sets of model results to an excel workbook with each CrossSection's 'fields' DataFrame in a separate sheet. Also write the ROW_edge_max DataFrame to a csv.
        kw:
            path - string, destination/filename for saved file"""
        self.results_export(**kw)
        self.ROW_edge_export(**kw)

    def results_export(self, **kw):
        """Write all of the cross section results to an excel workbook with each CrossSection's 'fields' DataFrame in a separate sheet.
        kw:
            path - string, destination/filename for saved file"""
        #path management
        fn = fields_funks._path_manage(self.name + '-all_results','.xlsx',**kw)
        #write results
        xlwriter = pd.ExcelWriter(fn, engine='xlsxwriter')
        for xs in self:
            xs.fields.to_excel(xlwriter, sheet_name=xs.sheet,
                    index_label='Distance (ft)')
        print('Full SectionBook results written to: %s' % fn)

    def ROW_edge_export(self, **kw):
        """Write max field results at ROW edges for each cross section to an excel or csv file. Default is csv.
        kw:
            file_type - string, accepts 'csv' or 'excel' (default csv)
            path - string, destination/filename for saved file
            xl - pandas ExcelWriter object, takes precedence over 'path'"""
        #be sure ROW_edge_results are current
        #self.compile_ROW_edge_results()
        #export
        c = ['Bmaxl','Bmaxr','Emaxl','Emaxr']
        h = ['Bmax - Left ROW Edge', 'Bmax - Right ROW Edge',
            'Emax - Left ROW Edge', 'Emax - Right ROW Edge']
        wo = False
        if('xl' in kw):
            wo = kw['xl']
        elif('file_type' in kw):
            file_type = kw['file_type']
            if(file_type[0] == '.'):
                file_type = file_type[1:]
            if(file_type == 'excel'):
                wo = fields_funks._path_manage(self.name + '-ROW_edge_results',
                        '.xlsx', **kw)
        if(wo):
            self.ROW_edge_max.to_excel(wo, index_label='Sheet', columns=c,
                                    header=h, sheet_name='ROW_edge_max')
        else:
            wo = fields_funks._path_manage(self.name + '-ROW_edge_results',
                    '.csv', **kw)
            self.ROW_edge_max.to_csv(wo, index_label='Sheet',
                    columns=c, header=h)
        if(not ('xl' in kw)):
            print('Maximum fields at ROW edges written to: %s' % wo)

class _IntegerIndexer(object):
    """Ancillary class for retrieval of items from a list in a parent object"""

    def __init__(self, L):
        """Accepts a list"""
        self._L = L

    def __getitem__(self, key):
        if(type(key) is not int):
            raise(EMFError("""can only receive integer indices"""))
        return(self._L[key])
