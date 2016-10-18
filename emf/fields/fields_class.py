from .. import np, pd, copy

from ..emf_class import EMFError

import fields_funks
import fields_plots
import fields_calcs
import FIELDS_io

class Conductor(object):
	"""Class representing a single conductor or power line."""

	def __init__(self, tag):
		"""
		args:
			tag - int or string, essentially just a name"""
		self.tag = tag #conductor labels
		self._freq = 60. #phase frequency
		self._x = None #x coordinate
		self._y = None #y coordinate
		self._subconds = None #number of subconductors per bundle
		self._d_cond = None #conductor diameter
		self._d_bund = None #bundle diameter
		self._V = None #line voltage
		self._I = None #line current
		self._phase = None #phase angle in degrees

	#---------------------------------------------------------------------------
	#PROPERTIES

	def _check_complete(self):
		"""Check if all Conductor variables have been set"""
		d = vars(self)
		keys = d.keys()
		for k in keys:
			if(d[k] is None):
				return(False, k)
		return(True, None)
	complete = property(_check_complete)

	def _check_to_float(self, x, prop):
		if(not fields_funks._is_number(x)):
			raise(EMFError("""
			Conductor property '%s' must be numeric.
			It cannot be set to: %s""" % (prop, repr(x))))
		return(float(x))

	def _check_to_int(self, x, prop):
		if(not fields_funks._is_number(x)):
			raise(EMFError("""
			Conductor property '%s' must be numeric.
			It cannot be set to: %s""" % (prop, repr(x))))
		elif(int(x) != float(x)):
			raise(EMFError("""
			Conductor property '%s' must be an integer.
			It cannot be set to: %s""" % (prop, repr(x))))
		return(int(x))

	def _get_freq(self): return(self._freq)
	def _set_freq(self, value):
		self._freq = self._check_to_float(value, 'freq')
	freq = property(_get_freq, _set_freq)

	def _get_x(self): return(self._x)
	def _set_x(self, value): self._x = self._check_to_float(value, 'x')
	x = property(_get_x, _set_x)

	def _get_y(self): return(self._y)
	def _set_y(self, value): self._y = self._check_to_float(value, 'y')
	y = property(_get_y, _set_y)

	def _get_subconds(self): return(self._subconds)
	def _set_subconds(self, value):
		self._subconds = self._check_to_int(value, 'subconds')
	subconds = property(_get_subconds, _set_subconds)

	def _get_d_cond(self): return(self._d_cond)
	def _set_d_cond(self, value):
		self._d_cond = self._check_to_float(value, 'd_cond')
	d_cond = property(_get_d_cond, _set_d_cond)

	def _get_d_bund(self): return(self._d_bund)
	def _set_d_bund(self, value):
		self._d_bund = self._check_to_float(value, 'd_bund')
	d_bund = property(_get_d_bund, _set_d_bund)

	def _get_V(self): return(self._V)
	def _set_V(self, value): self._V = self._check_to_float(value, 'V')
	V = property(_get_V, _set_V)

	def _get_I(self): return(self._I)
	def _set_I(self, value): self._I = self._check_to_float(value, 'I')
	I = property(_get_I, _set_I)

	def _get_phase(self): return(self._phase)
	def _set_phase(self, value):
		self._phase = self._check_to_float(value, 'phase')
	phase = property(_get_phase, _set_phase)

	#---------------------------------------------------------------------------
	#METHODS

	def __str__(self):
		"""quick and dirty printing"""
		v = vars(self)
		keys = v.keys()
		s = '\n'
		for k in keys:
			s += str(k) + ': ' + repr(v[k]) + '\n'
		return(s)

	def copy(self):
		return(copy.deepcopy(self))

class CrossSection(object):
	"""Class that organizes Conductor objects and stores other input
	information for a power line cross section. Includes plotting methods
	for the fields results and exporting methods for the results."""

	def __init__(self, sheet, *args):
		"""CrossSection must be initialized with a sheet string, and a list
		of Conductor objects can optionally be passed in"""
		self.sheet = sheet #mandatory, short, template worksheet name
		self.tag = None #identifier linking multiple CrossSection objects
		self.title = '' #longer form, used for plotting text
		self.soil_resistivity = 100. #?
		self.max_dist = None #maximum simulated distance from the ROW center
		self.step = 1 #step size for calculations
		self.sample_height = 3. #uniform sample height
		self.lROW = None #exact coordinate of the left ROW edge
		self.rROW = None #exact coordinate of the left ROW edge
		self.conds = [] #list of Conductor objects
		#dictionary mapping Conductor tags to Conductor objects
		self.tag2condidx = {}
		#DataFrame storing results, populated with calculate_fields()
		self.fields = pd.DataFrame(columns=['Bx','By','Bprod','Bmax',
											'Ex','Ey','Eprod','Emax'])
		#add conductors if they're passed in
		if(len(args) == 1):
			for c in args[0]:
				self.add_conductor(c)

	#---------------------------------------------------------------------------
	#PROPERTIES

	def _get_tags(self): return([xs.tag for xs in self.conds])
	tags = property(_get_tags)

	def _get_N_sample(self): return(1 + int(np.ceil(2*self.max_dist/self.step)))
	N_sample = property(_get_N_sample)

	def _get_hot(self): return([c for c in self.conds if c.V != 0])
	hot = property(_get_hot)

	def _get_gnd(self): return([c for c in self.conds if c.V == 0])
	gnd = property(_get_gnd)

	def _get_x_sample(self):
		return(np.linspace(-self.max_dist, self.max_dist, num=self.N_sample))
	x_sample = property(_get_x_sample)

	def _get_y_sample(self):
		return(self.sample_height*np.ones((self.N_sample,), dtype=float))
	y_sample = property(_get_y_sample)

	def _get_x(self): return(np.array([c.x for c in self.conds], dtype=float))
	x = property(_get_x)

	def _get_y(self): return(np.array([c.y for c in self.conds], dtype=float))
	y = property(_get_y)

	def _get_lROWi(self):
		d = np.absolute(self.x_sample - self.lROW)
		return(max(np.where(d == np.min(d))[0]))
	lROWi = property(_get_lROWi)

	def _get_rROWi(self):
		d = np.absolute(self.x_sample - self.rROW)
		return(min(np.where(d == np.min(d))[0]))
	rROWi = property(_get_rROWi)

	def _get_subconds(self):
		return(np.array([c.subconds for c in self.conds], dtype=float))
	subconds = property(_get_subconds)

	def _get_d_cond(self):
		return(np.array([c.d_cond for c in self.conds], dtype=float))
	d_cond = property(_get_d_cond)

	def _get_d_bund(self):
		return(np.array([c.d_bund for c in self.conds], dtype=float))
	d_bund = property(_get_d_bund)

	def _get_V(self):
		return(np.array([c.V for c in self.conds], dtype=float))
	V = property(_get_V)

	def _get_I(self):
		return(np.array([c.I for c in self.conds], dtype=float))
	I = property(_get_I)

	def _get_phase(self):
		return(np.array([c.phase for c in self.conds], dtype=float))
	phase = property(_get_phase)

	def _get_ROW_edge_fields(self):
		return(self.fields.iloc[[self.lROWi, self.rROWi]])
	ROW_edge_fields = property(_get_ROW_edge_fields)

	def _check_complete(self):
		"""Check that all variables in each conductor in the CrossSection
		are set, and not left with the initial None values.
		returns:
			b - True if all Conductors are complete, False if not
			t - tag of first incomplete Conductor, if any, else None
			v - variable name of first unset variable in first incomplete
				Conductor, if any, else None"""
		for c in self.conds:
			b, v = c.complete
			if(b is False):
				return(b, c.tag, v)
		return(True, None, None)
	complete = property(_check_complete)

	#---------------------------------------------------------------------------

	def __str__(self):
		"""quick and dirty printing"""
		v = vars(self)
		keys = v.keys()
		s = '\n'
		for k in keys:
			if(k != 'fields'):
				s += str(k) + ': ' + repr(v[k]) + '\n'
		s += '\ninspect self.fields separately to see field simulation results\n'
		return(s)

	def __getitem__(self, key):
		#check in the hot dict
		try:
			idx = self.tag2condidx[key]
		except(KeyError):
			pass
		else:
			return(self.conds[idx])
		#return None if no xs is found
		return(None)

	def add_conductor(self, cond):
		#check if the Conductor is complete
		b, v = cond.complete
		if(not b):
			raise(EMFError("""
			Cannot add Conductor "%s" to the CrossSection because it
			is not complete. The parameter "%s" is not set."""
			% (cond.tag, v[1:])))
		#check if the tag has already been used
		if(cond.tag in self.tag2condidx):
			raise(EMFError("""
			A Conductor with tag "%s" is already in CrossSection "%s".
			Another Conductor with the same tag cannot be added"""
			% (cond.tag, self.sheet)))
		#see if the conductor is grounded
		if(cond.V == 0):
			#probably shouldn't have current if grounded
			if(cond.I != 0):
				print("""
				Conductor with tag "%s" in CrossSection "%s"
				is grounded (V = 0) but has nonzero current?"""
				% (cond.tag, self.sheet))
		#add to self.conds and indexing dict
		self.tag2condidx[cond.tag] = len(self.conds)
		self.conds.append(copy.deepcopy(cond))

	def remove_conductor(self, key):
		#check in the hot dict
		try:
			idx = self.tag2condidx[key]
		except(KeyError):
			pass
		else:
			self.conds.pop(idx)
			self.tag2condidx.pop(key)

	def calculate_fields(self):
		"""Calculate electric and magnetic fields across the ROW and store the
		results in the self.fields DataFrame"""
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
		self.fields = pd.DataFrame({'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
									'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax},
									index=self.x_sample)
		#return the fields dataframe
		return(self.fields)

	def compare_DAT(self, DAT_path, **kw):
		"""Load a FIELDS output file (.DAT), find absolute and percentage
		differences between it and the CrossSection object's results.
		args:
			DAT_path - path of FIELDS results file
		kw:
			save - boolean, toggle whether output Panel and figures are saved,
				   also saves figures with error and comparison plots
			path - string, destination of saved files, will force save == True
			round - int, round the results in self.fields to a certain
					number of digits in an attempt to exactly match the
					FIELDS results, which are printed only to the
					thousandths digit
			truncate - bool, truncate results after the thousandths digit
		returns:
			pan - pandas Panel with DAT results, results of this code,
				  the absolute error between, and the relative error between"""
		#load the .DAT file into a dataframe
		df = FIELDS_io.read_DAT(DAT_path)
		#check dataframe shape compatibility
		if(df.shape != self.fields.shape):
			raise(EMFError("""
			self.fields in CrossSection named "%s" and the imported .DAT
			DataFrame have different shapes. Be sure to target the correct
			.DAT file and that it has compatible DIST values.""" % self.sheet))
		#prepare a dictionary to create a Panel
		if(('round' in kw) and ('truncate' in kw)):
			raise(FLDError("""
			Cannot both round and truncate for DAT comparison. Choose either
			rounding or truncation."""))
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
		frames = ['FIELDS_DAT_results', 'python_results',
				'Absolute Difference', 'Percent Difference']
		pan = pd.Panel(data={frames[0] : df,
				frames[1] : f,
				frames[2] : f - df,
				frames[3] : 100*(f - df)/f})
		#write data and save figures if called for
		if('path' in kw):
			kw['save'] = True
		if('save' in kw):
			if(kw['save']):
				fn = fields_funks._path_manage(self.sheet + '-DAT_comparison',
					'.xlsx', **kw)
				xl = pd.ExcelWriter(fn, engine='xlsxwriter')
				for f in frames:
					pan[f].to_excel(xl, index_label='x', sheet_name=f)
				xl.save()
				print('DAT comparison book saved to: "%s"' % fn)
				#make plots of the absolute and percent error
				figs = fields_plots._plot_DAT_comparison(self, pan, **kw)
		#return the Panel
		return(pan)

	def sample(self, *args):
		"""Get a random Conductor from the CrossSection
		args:
			an input integer determines the number of random xss fetched"""
		c = self.conds
		if(len(args) == 0):
			return(self.conds[np.random.randint(len(c))])
		else:
			r = np.random.randint(len(c), size=args[0])
			return([c[i] for i in r])

	def copy(self):
		return(copy.deepcopy(self))

class SectionBook(object):
	"""Top level class organizing a group of CrossSection objects. Uses a
	dictionary to track CrossSections in a list and provide a convenient
	__getitem__ method than gets CrossSections by their name. Also tracks
	maximum field results at the ROW edges of each CrossSection added,
	provides a plotting method for CrossSection groups, and provides
	exporting methods."""

	def __init__(self, name, *args):
		self.name = name #mandatory identification field
		self.xss = [] #list of cross section objects
		self.sheet2idx = dict() #mapping dictionary for CrossSection retrieval
		self.i = _IntegerIndexer(self)
		#add CrossSections if they're passed in
		if(len(args) == 1):
			for xs in args[0]:
				self.add_section(xs)

	#---------------------------------------------------------------------------
	#PROPERTIES

	def _get_sheets(self): return([xs.sheet for xs in self.xss])
	sheets = property(_get_sheets, None, None, 'list of CrossSection sheets')

	def _get_tags(self): return([xs.tag for xs in self.xss])
	tags = property(_get_tags, None, None, 'list of CrossSection tags')

	def _get_titles(self): return([xs.title for xs in self.xss])

	def _get_tag_groups(self):
		u = list(set(self.tags)) #get unique CrossSection tags
		tag_groups = []
		tag_groups = [[] for i in range(len(u))]
		for i in range(len(self.xss)):
			tag_groups[u.index(self.xss[i].tag)].append(self.xss[i])
		return(tag_groups)
	tag_groups = property(_get_tag_groups, None, None,
			'a list of lists of grouped CrossSection objects')

	def _get_ROW_edge_max(self):
		#gather ROW edge results
		L = len(self.xss)
		El,Er,Bl,Br = np.zeros((L,)),np.zeros((L,)),np.zeros((L,)),np.zeros((L,))
		for i in range(L):
			xs = self.i[i]
			Bl[i] = xs.fields['Bmax'].iat[xs.lROWi]
			Br[i] = xs.fields['Bmax'].iat[xs.rROWi]
			El[i] = xs.fields['Emax'].iat[xs.lROWi]
			Er[i] = xs.fields['Emax'].iat[xs.rROWi]
		#construct DataFrame
		df = pd.DataFrame(data={
			'sheet': self.sheets, 'title': [xs.title for xs in self],
			'Bmaxl': Bl, 'Emaxl': El, 'Bmaxr': Br, 'Emaxr': Er})
		return(df)
	ROW_edge_max = property(_get_ROW_edge_max, None, None,
			'DataFrame with maximum field magnitudes for each CrossSection')

	def _check_complete(self):
		"""Check that all CrossSections in the SectionBook are complete
		returns:
			b - True if all xcs are complete, False if not
			sheet - sheet of first incomplete xs, if any, else None
			ctag - tag of first incomplete Conductor in first incomplete
					xs, if any, else None
			v - variable name of first unset variable in first incomplete
				Conductor in first incomplete xs, if any, else None"""
		for xs in self:
			b, ctag, v = xs.complete
			if(b is False):
				return(b, xs.sheet, ctag, v)
		return(True, None, None, None)
	complete = property(_check_complete)

	#---------------------------------------------------------------------------

	def __getitem__(self, key):
		"""Index the SectionBook by CrossSection names"""
		try:
			idx = self.sheet2idx[key]
		except(KeyError):
			return(False)
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
		"""quick and dirty printing"""
		v = vars(self)
		keys = v.keys()
		s = '\n'
		for k in keys:
			s += str(k) + ': ' + repr(v[k]) + '\n'
		return(s)

	def sample(self, *args):
		"""Get a random CrossSection from the SectionBook
		args:
			an input integer determines the number of random xss fetched"""
		if(len(args) == 0):
			return(self.xss[np.random.randint(len(self))])
		else:
			r = np.random.randint(len(self), size=args[0])
			return([self.xss[i] for i in r])

	def add_section(self, xs):
		"""Add a CrossSection to the book. Doing so by directly altering
		self.xss will make the CrossSections inaccessible by __getitem__
		and make the group plotting functions impossible, so don't do that.
		Use this method instead."""
		#Prevent adding CrossSections with the same names
		if(xs.sheet in self.sheet2idx):
			raise(EMFError("""
			CrossSection name "%s" already exists in the SectionBook.
			Duplicate names would cause collisions in the lookup dictionary
			(self.sheet2idx). Use a different name.""" % xs.sheet))
		else:
			xs.calculate_fields()
			self.sheet2idx[xs.sheet] = len(self.xss)
			self.xss.append(copy.deepcopy(xs))

	def remove_section(self, sheet):
		"""Remove a CrossSection from the SectionBook by providing its
		sheet string
		args:
			sheet - sheet string of CrossSection to remove"""
		#check sheet is in self
		if(sheet not in self.sheet2idx):
			raise(EMFError("""
			Sheet %s cannot be removed because it does not exist in the
			SectionBook object""" % sheet))
		#delete from self.xss list
		idx = self.sheet2idx[sheet]
		self.xss.pop(idx)
		#update
		self.sheet2idx = dict(zip([xs.sheet for xs in self.xss], xss))

	def export(self, **kw):
		"""Write results to an excel workbook and ROW edge results to a csv
		kw:
			path - string, destination/filename for saved file"""
		self.results_export(**kw)
		self.ROW_edge_export(**kw)

	def results_export(self, **kw):
		"""Write all of the cross section results to an excel workbook"""
		#path management
		fn = fields_funks._path_manage(self.name + '-all_results','.xlsx',**kw)
		#write results
		xlwriter = pd.ExcelWriter(fn, engine='xlsxwriter')
		for xs in self:
			xs.fields.to_excel(xlwriter, sheet_name=xs.sheet, index_label='x')
		print('Full SectionBook results written to: %s' % fn)

	def ROW_edge_export(self, **kw):
		"""Write max field results at ROW edges for each cross section to
		an excel or csv file. Default is csv.
		kw:
			file_type - string, accepts 'csv' or 'excel' (default csv)
			path - string, destination/filename for saved file
			xl - pandas ExcelWriter object, takes precedence over 'path'"""
		#be sure ROW_edge_results are current
		#self.compile_ROW_edge_results()
		#export
		c = ['sheet','title','Bmaxl','Bmaxr','Emaxl','Emaxr']
		h = ['Cross-Section Sheet','Cross-Section Title',
				'Bmax - Left ROW Edge', 'Bmax - Right ROW Edge',
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
			self.ROW_edge_max.to_excel(wo, index=False, columns=c,
									header=h, sheet_name='ROW_edge_max')
		else:
			wo = fields_funks._path_manage(self.name + '-ROW_edge_results',
					'.csv', **kw)
			self.ROW_edge_max.to_csv(wo, index=False, columns=c, header=h)
		if(not ('xl' in kw)):
			print('Maximum fields at ROW edges written to: %s' % wo)

	def update_fields(self):
		"""Executes all of the update functions"""
		for xs in self:
			xs.calculate_fields()

class _IntegerIndexer(object):
	"""Ancillary class for retrieval of CrossSections in SectionBooks by
	their integer position in the SectionBook"""

	def __init__(self, sb):
		self._sb = sb

	def __getitem__(self, key):
		if(type(key) is not int):
			raise(EMFError("""
			SectionBook.i can only receive integer indices."""))
		return(self._sb.xss[key])
