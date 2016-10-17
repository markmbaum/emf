from .. import np, pd, os, copy, interpn

from ..emf_class import EMFError

import subcalc_funks

class Model(object):
	"""Model objects store the results of magnetic field modeling generated
	by the SubCalc program. A Model can be created by passing grid results
	directly or by passing a dictionary of grid results. The function
	subcalc.load_model is the best way to generate a Model object from
	a .REF file containing SubCalc results. The Model object can be saved
	to a more flexible (and smaller) excel file with Model.export(). Then
	the excel file can be read back into a Model object using
	subcalc.load_model again, but targeting the excel file.
		Model objects have a 'Bkey' property that determines which component
	of the magnetic field results is accessed by the Model.B property. For
	example, when Model.Bkey == 'Bmax' all operations involving the grid
	of magnetic field results accessed by Model.B will operate on the
	'max' field component. When Bmax == 'Bz' all operations deal with the
	vertical 'z' component of the field, and so on. This includes plotting.
		Footprints of objets in the Model domain like buildings, the power
	lines, fences, etc. can also be stored in Model objects. Data for these
	objects must be saved in csv template files. The path of the footprint
	csv files can be passed to subcalc.load_model for automatic inclusion in
	a newly generated Model object or it can be passed to an existing Model
	with Model.load_footprints. The footprint data is stored in Footprint
	objects that have very little functionality and are mostly just
	organizational objects.
		Several methods are available for interpolating new values from
	the grid of field results: Model.interp, Model.resample, and
	Model.crosssection.
		Contour plots of the results can be automatically generated with
	subcalc.plot_contour(Model) and colormesh plots can be automatically
	generated with subcalc.plot_pcolormesh(Model)."""

	def __init__(self, *args, **kw):
		"""Grid data must be passed in
		args:
			either:
				X - 2D array of x coordinates
				Y - 2D array of y coordinats
				B - 2D array of magnetic field magnitudes
				info - dict, optional, dictionary of model metadata
			or:
				data - dict, dictionary with X, Y, and B grids, should be
							keyed by 'X','Y','Bmax','Bres','Bx','By','Bz'
				info - dict, optional, dictionary of model metadata
		kw:
			Bkey - str, selects which component of field results to use
					or defines which component is passed as grid arrays
					('Bmax', 'Bres', 'Bx', 'By', 'Bz')."""

		largs = len(args)
		if(largs <= 2):
			#check args[0]
			if(type(args[0]) is not dict):
				raise(EMFError("""
				The first argument to Model() must be a dictionary of
				model result information when passing 1 or 2 arguments to
				initialize the Model, not %s""" % type(args[0])))
			#check keys
			s = set(['X','Y','Bmax','Bres','Bx','By','Bz'])
			if(any([(i not in s) for i in args[0].keys()])):
				print args[0].keys()
				raise(EMFError("""
				If passing a dictionary to initialize a Model object, the
				dict must have the following keys only:
					%s""" % str(s)))
			#store data
			self._grid = args[0]
			#deal with Bkey
			if('Bkey' in kw):
				self.Bkey = kw['Bkey']
			else:
				self.Bkey = 'Bmax'
			#store the info dict if present
			if(largs == 2):
				if(type(args[1]) is not dict):
					raise(EMFError("""
					The fourth argument to Model() must be a dictionary of
					model result information, not %s""" % type(args[1])))
				self.info = args[1]
			else:
				self.info = None
		elif(largs <= 4):
			#check input types
			mgs = """
			If passing three or four arguments to initialize a Model, the first
			three arguments must be 2D numpy arrays representing X, Y, and B
			grids respectively, each with the same shape."""
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
					raise(EMFError("""
					The fourth argument to Model() must be a dictionary of
					model result information, not type %s""" % type(args[3])))
				self.info = args[3]
			else:
				self.info = None

		#other reference objects in the model, like substation boundaries,
		#stored in a list of Footprint objects
		self.footprint_df = None #DataFrame of Footprint information
		self.footprints = []
		#angle of the northern direction with respect to the grid
		#   where 0 degrees is the positive y axis and clockwise is increasing
		self._north_angle = None

	#---------------------------------------------------------------------------
	#properties
	def _get_B(self):
		return(self._grid[self._Bkey])
	B = property(_get_B, None, None, '2D grid of magnetic field results')

	def _get_Bkey(self):
		return(self._Bkey)
	def _set_Bkey(self, value):
		if((not (value in self._grid)) or (value == 'X') or (value == 'Y')):
			k = self._grid.keys()
			k.remove('X')
			k.remove('Y')
			raise(EMFError("""
			Bkey must be set to one of the following elements:
				%s""" % str(k)))
		else:
			self._Bkey = value
	Bkey = property(_get_Bkey, _set_Bkey, None,
			'Component of magnetic field accessed by the B property')

	def _get_X(self):
		return(self._grid['X'])
	X = property(_get_X, None, None, '2D grid of reference grid x coordinates')

	def _get_Y(self):
		return(self._grid['Y'])
	Y = property(_get_Y, None, None, '2D grid of reference grid y coordinates')

	def _get_x(self):
		return(self.X[0,:])
	x = property(_get_x, None, None,
			'Unique x values in model grid (column positions)')

	def _get_y(self):
		return(self.Y[:,0])
	y = property(_get_y, None, None,
			'Unique y values in model grid (row positions)')

	def _get_xmax(self):
		return(np.max(self.x))
	xmax = property(_get_xmax)

	def _get_xmin(self):
		return(np.min(self.x))
	xmin = property(_get_xmin)

	def _get_ymax(self):
		return(np.max(self.y))
	ymax = property(_get_ymax)

	def _get_ymin(self):
		return(np.min(self.y))
	ymin = property(_get_ymin)

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

	def load_footprints(self, footprint_info, **kw):
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

	def cross_section(self, p1, p2, **kw):
		"""Interpolate the field along a line between two points
		args:
			p1 - iterable, an x-y pair
			p2 - iterable, an x-y pair
		kw:
			n - integer, number of points sampled (default 1000)
		returns:
			x - array, x coordinates of interpolated values
			y - array, y coordinates of interpolated values
			B_interp - array, interpolated field values"""
		#check point lengths
		if((len(p1) != 2) or (len(p2) != 2)):
			raise(EMFError('Points must consist of two values (xy pairs).'))
		#check kw
		if('n' in kw):
			n = kw['n']
		else:
			n = 1000
		#create x and y vectors
		x, y = np.linspace(p1[0], p2[0], n), np.linspace(p1[1], p2[1], n)
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
		"""Check if an x,y coordinate pair is inside the Model grid
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
		"""Resample the model grid along a new number of x,y values, a new
		selection of x,y values, or a new number of total values
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
			mod_resample - a new Model object containing the resampled grid"""
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
		#some arrays have to be flipped to conform to conventions of interpn
		B_resample = interpn((self.y[::-1], self.x),
								self.B[::-1,:],
								(Y[::-1,:], X))

		#return with re-flipped results
		mod_resample = Model(X, Y, B_resample[::-1,:],
				copy.deepcopy(self.info),
				Bkey=self.Bkey)
		mod_resample.load_footprints(self.footprint_df)
		return(mod_resample)

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
		"""Export the grid data and accompanying info to an excel file with
		tabs for each Bfield component, another for the info dict, and a
		final one for footprints if they're present
		kw:
			path - string, output destination/filename for workbook"""
		#get appropriate export filename
		fn = os.path.basename(self.info['REF_path'])
		fn = subcalc_funks._path_manage(fn, '.xlsx', **kw)
		#create excel writing object
		xl = pd.ExcelWriter(fn, engine='xlsxwriter')
		#write grid data
		for k in self._grid:
			if((k != 'X') and (k != 'Y')):
				pd.DataFrame(self._grid[k], columns=self.x, index=self.y
						).to_excel(xl, sheet_name=k)
		#write model information if present
		if(self.info is not None):
			pd.DataFrame([self.info[k] for k in self.info], index=self.info.keys(),
					columns=['Parameter Value']).sort_index().to_excel(
							xl, sheet_name='info', index_label='Parameter Name')
		#write footprint DataFrame if present
		if(self.footprint_df is not None):
			self.footprint_df.to_excel(xl, sheet_name='footprints', index=False)
		#save and print
		xl.save()
		print('model saved to: %s' % fn)

class Footprint(object):

	def __init__(self, name, x, y, power_line, of_concern, draw_as_loop, group):
		"""
		args:
			name - string, the name of the Footprint, i.e. "Substation"
			x - iterable, x coordinates of Footprint
			y - iterable, y coordinates of Footprint
			power_line - bool, True indicates the footprint corresponds to
							modeled
			of_concern - bool, True if the Footprint represents an area
						that is potentially concerned about EMF (homes)
			draw_as_loop - bool, True if the footprint should be plotted
							as a closed loop
			group - string, group strings that are identical between
					Footprint objects designate them as part of the same
					group for plotting"""
		#check x and y are the same length
		if(len(x) != len(y)):
			raise(EMFError("""
			Footprints must have the same number of x and y values to form
			spatial coordinates"""))
		#set attributes
		self.name = name #string
		self._x = [float(i) for i in x]
		self._y = [float(i) for i in y]
		self.power_line = power_line #bool
		self.of_concern = of_concern #bool
		self.draw_as_loop = draw_as_loop #bool
		self._group = group #string

	def _get_x(self):
		if(self.draw_as_loop):
			return(self._x + [self._x[0]])
		else:
			return(self._x)
	def _set_x(self, value):
		self._x = value
	x = property(_get_x, _set_x, None, 'x coordinates of Footprint vertices')

	def _get_y(self):
		if(self.draw_as_loop):
			return(self._y + [self._y[0]])
		else:
			return(self._y)
	def _set_y(self, value):
		self._y = value
	y = property(_get_y, _set_y, None, 'y coordinates of Footprint vertices')

	def _get_group(self):
		if(pd.isnull(self._group) or (self._group is None)):
			return('')
		else:
			return(self._group)
	def _set_group(self, value):
		self._group = value
	group = property(_get_group, _set_group, None, 'Group name')

	def __str__(self):
		"""quick and dirty printing"""
		v = vars(self)
		keys = v.keys()
		s = '\n'
		for k in keys:
			s += str(k) + ': ' + str(v[k]) + '\n'
		return(s)
