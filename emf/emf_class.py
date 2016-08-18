import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import emf_funks
import emf_plots
import emf_calcs
import fields_io

class EMFError(Exception):
    """Exception class for emf specific errors"""
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return(self.message)

class Conductor:
    """Class representing a single conductor or power line."""

    def __init__(self):
        self.tag = None #conductor label
        self.freq = 60. #phase frequency
        self.x = None #x coordinate
        self.y = None #y coordinate
        self.subconds = None #number of subconductors per bundle
        self.d_cond = None #conductor diameter
        self.d_bund = None #bundle diameter
        self.V = None #line voltage
        self.I = None #line current
        self.phase = None #phase angle

    def __str__(self):
        """quick and dirty printing"""
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            s += str(k) + ': ' + str(v[k]) + '\n'
        return(s)

class CrossSection:
    """Class that organizes Conductor objects and stores other input
    information for a power line cross section. Includes plotting methods
    for the fields results and exporting methods for the results."""

    def __init__(self, name):
        self.name = name #mandatory, short, usually template sheet name
        self.title = '' #must be short, used as FLD file names
        self.tag = None #identifier linking multiple CrossSection objects
        self.subtitle = '' #longer form, used for plotting text
        self.soil_resistivity = 100. #?
        self.max_dist = None #maximum simulated distance from the ROW center
        self.step = None #step size for calculations
        self.sample_height = 3. #uniform sample height
        self.lROW = None #exact coordinate of the left ROW edge
        self.lROWi = None #integer index of x_sample closest to self.lROW
        self.rROW = None #exact coordinate of the left ROW edge
        self.rROWi = None #integer index of x_sample closest to self.rROW
        self.hot = [] #list of Conductor objects with nonzero voltage
        self.gnd = [] #list of Conductor objects with zero voltage
        #arrays with unlabeled conductor data for fast passing to emf_calcs
        #need to be updated with update_arrays() if conductors change
        #these arrays store hot conductor info first, then the grounded lines
        self.x = np.empty((0,)) #x coordinates
        self.y = np.empty((0,)) #y coordinates
        self.subconds = np.empty((0,)) #number of subconductors per bundles
        self.d_cond = np.empty((0,)) #conductor diameters
        self.d_bund = np.empty((0,)) #bundle diameters
        self.V = np.empty((0,)) #ilne voltages
        self.I = np.empty((0,)) #line currents
        self.phase = np.empty((0,)) #phase angles
        self.x_sample = np.empty((0,)) #x coordinates of sample points
        self.y_sample = np.empty((0,)) #y coordinates of sample points
        #DataFrame storing results, populated with calculate_fields()
        self.fields = pd.DataFrame(columns = ['Bx','By','Bprod','Bmax',
                                            'Ex','Ey','Eprod','Emax'])

    def __str__(self):
        """quick and dirty printing"""
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            if(k != 'fields'):
                s += str(k) + ': ' + str(v[k]) + '\n'
        s += '\ninspect self.fields separately to see field simulation results\n'
        return(s)

    def update_arrays(self):
        """Populate the CrossSection objects numpy array attributes"""
        #calculate sample point coordinates
        N = 1 + 2*self.max_dist/self.step
        self.x_sample = np.linspace(-self.max_dist, self.max_dist, num = N)
        self.y_sample = self.sample_height*np.ones((N,), dtype = float)
        #update ROW edge index variables
        #if ROW edge lies between two sample points, use the one closer to zero
        d = np.absolute(self.x_sample - self.lROW)
        self.lROWi = max(np.where(d == np.min(d))[0])
        #if ROW edge lies between two sample points, use the one closer to zero
        d = np.absolute(self.x_sample - self.rROW)
        self.rROWi = min(np.where(d == np.min(d))[0])
        #assemble all the conductor data in arrays for calculations
        conds = self.hot + self.gnd
        self.x = np.array([c.x for c in conds], dtype = float)
        self.y = np.array([c.y for c in conds], dtype = float)
        self.subconds = np.array([c.subconds for c in conds], dtype = float)
        self.d_cond = np.array([c.d_cond for c in conds], dtype = float)
        self.d_bund = np.array([c.d_bund for c in conds], dtype = float)
        self.V = np.array([c.V for c in conds], dtype = float)
        self.I = np.array([c.I for c in conds], dtype = float)
        self.phase = np.array([c.phase for c in conds], dtype = float)

    def calculate_fields(self):
        """Calculate electric and magnetic fields across the ROW and store the
        results in the self.fields DataFrame. This function also initializes
        many CrossSection variables, so it's wise to run it as soon as the
        all CrossSection data is entered or whenever it changes."""
        #populate flat array variables with all Conductor data
        self.update_arrays()
        #calculate magnetic field
        Bx, By = emf_calcs.B_field(self.x, self.y, self.I, self.phase,
            self.x_sample, self.y_sample)
        Bx, By, Bprod, Bmax = emf_calcs.phasors_to_magnitudes(Bx, By)
        #calculate electric field
        Ex, Ey = emf_calcs.E_field(self.x, self.y, self.subconds, self.d_cond,
            self.d_bund, self.V, self.phase, self.x_sample, self.y_sample)
        Ex, Ey, Eprod, Emax = emf_calcs.phasors_to_magnitudes(Ex, Ey)
        #store the values
        self.fields = pd.DataFrame({'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
                                    'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax},
                                    index = self.x_sample)
        #return the fields dataframe
        return(self.fields)

    def ROW_edge_fields(self):
        """Return the values of fields calculations at the left and right
        ROW edges of the cross-section.
        returns:
            left - pandas Series with left ROW edge results
            right - pandas Series with right ROW edge results"""
        return(self.fields.iloc[[self.lROWi, self.rROWi]])

    def compare_DAT(self, DAT_path, **kwargs):
        """Load a FIELDS output file (.DAT), find absolute and percentage
        differences between it and the CrossSection object's results.
        args:
            DAT_path - path of FIELDS results file
        kwargs:
            save - boolean, toggle whether output panel and figures are saved
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
        df = fields_io.read_DAT(DAT_path)
        #check dataframe shape compatibility
        if(df.shape != self.fields.shape):
            raise(EMFError("""
            self.fields in CrossSection named "%s" and the imported .DAT
            DataFrame have different shapes. Be sure to target the correct
            .DAT file and that it has compatible DIST values.""" % self.name))
        #prepare a dictionary to create a Panel
        if(('round' in kwargs) and ('truncate' in kwargs)):
            raise(FLDError("""
            Cannot both round and truncate for DAT comparison. Choose either
            rounding or truncation."""))
        elif('round' in kwargs):
            f = self.fields.round(kwargs['round'])
        elif('truncate' in kwargs):
            if(kwargs['truncate']):
                f = self.fields.copy(deep = True)
                for c in f.columns:
                    for i in f.index:
                        f[c].loc[i] = float('%.3f' % f[c].loc[i])
        else:
            f = self.fields
        comp = {'FIELDS_DAT_results' : df,
                'emf_results' : f,
                'Absolute Difference' : f - df,
                'Percent Difference' : 100*(f - df)/f}
        #wrap comparison DataFrames in a Panel
        pan = pd.Panel(data = comp)
        #write data and save figures if called for
        if('path' in kwargs):
            kwargs['save'] = True
        if('save' in kwargs):
            if(kwargs['save']):
                fn = emf_funks._path_manage(self.name + '-DAT_comparison',
                    '.xlsx', **kwargs)
                pan.to_excel(fn, index_label = 'x')
                print('DAT comparison book saved to: "%s"' % fn)
                #make plots of the absolute and percent error
                figs = emf_plots.plot_DAT_comparison(self, pan, **kwargs)
        #return the Panel
        return(pan)

    def ion(self):
        emf_plots.ion()

    def show(self):
        emf_plots.show()

class SectionBook:
    """Top level class organizing a group of CrossSection objects. Uses a
    dictionary to track CrossSections in a list and provide a convenient
    __getitem__ method than gets CrossSections by their name. Also tracks
    maximum field results at the ROW edges of each CrossSection added,
    provides a plotting method for CrossSection groups, and provides
    exporting methods."""

    def __init__(self, name):
        self.name = name #mandatory identification field
        self.xcs = [] #list of cross section objects
        self.name2idx = dict() #mapping dictionary for CrossSection retrieval
        self.names = [] #list of CrossSection names
        self.tags = [] #collection of CrossSection tags
        self.tag_groups = [[]] #groups of CrossSection indices with identical tags
        #DataFrame of maximum fields at ROW edges
        self.ROW_edge_max = pd.DataFrame(columns = ['name','title',
                                            'Bmaxl','Bmaxr','Emaxl','Emaxr'])

    def __getitem__(self, key):
        """Index the SectionBook by CrossSection names"""
        try:
            idx = self.name2idx[key]
        except(KeyError):
            return(False)
        else:
            return(self.xcs[idx])

    def __iter__(self):
        """Iteration over all CrossSections in the SectionBook"""
        for xc in self.xcs:
            yield(xc)

    def __len__(self):
        """Length of SectionBook is the number of CrossSections in it"""
        return(len(self.xcs))

    def __str__(self):
        """quick and dirty printing"""
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            s += str(k) + ': ' + str(v[k]) + '\n'
        return(s)

    def i(self, idx):
        """Get a CrossSection object by it's numeric index in self.xcs"""
        return(self.xcs[idx])

    def sample(self, *args):
        """Get a random CrossSection from the SectionBook
        args:
            any input integer determines the number of random xcs fetched"""
        if(len(args) == 0):
            return(self.xcs[np.random.randint(len(self))])
        else:
            r = np.random.randint(len(self), size = args[0])
            return([self.xcs[i] for i in r])

    def add_section(self, xc):
        """Add a CrossSection to the book. Doing so by directly altering
        self.xcs will make the CrossSections inaccessible by __getitem__
        and make the group plotting functions impossible, so don't do that
        and use this method instead."""
        #Prevent adding CrossSections with the same names
        if(xc.name in self.names):
            raise(EMFError("""
            CrossSection name "%s" already exists in the SectionBook.
            Duplicate names would cause collisions in the lookup dictionary
            (self.name2idx). Use a different name.""" % xc.name))
        else:
            self.name2idx[xc.name] = len(self.xcs)
            self.xcs.append(xc)
            self.names.append(xc.name)
            self.tags.append(xc.tag)

    def ROW_edge_export(self, **kwargs):
        """Write max field results at ROW edges for each cross section to
        an excel or csv file. Default is csv.
        kwargs:
            file_type - string, accepts 'csv' or 'excel'
            path - string, destination/filename for saved file
            xl - pandas ExcelWriter object, takes precedence over 'path'"""
        #be sure ROW_edge_results are current
        #self.compile_ROW_edge_results()
        #export
        c = ['name','title','Bmaxl','Bmaxr','Emaxl','Emaxr']
        h = ['Cross-Section Name','Cross-Section Title','Bmax - Left ROW Edge',
                'Bmax - Right ROW Edge', 'Emax - Left ROW Edge',
                'Emax - Right ROW Edge']
        wo = False
        if('xl' in kwargs):
            wo = kwargs['xl']
        elif('file_type' in kwargs):
            file_type = kwargs['file_type']
            if(file_type[0] == '.'):
                file_type = file_type[1:]
            if(file_type == 'excel'):
                wo = emf_funks._path_manage(self.name + '-ROW_edge_results',
                        '.xlsx', **kwargs)
        if(wo):
            self.ROW_edge_max.to_excel(wo, index = False, columns = c,
                                    header = h, sheet_name = 'ROW_edge_max')
        else:
            wo = emf_funks._path_manage(self.name + '-ROW_edge_results',
                '.csv', **kwargs)
            self.ROW_edge_max.to_csv(wo, index = False, columns = c, header = h)
        if(not ('xl' in kwargs)):
            print('Maximum fields at ROW edges written to: "%s"' % repr(wo))

    def results_export(self, **kwargs):
        """Write all of the cross section results to an excel workbook"""
        #path management
        fn = emf_funks._path_manage(self.name + '-all_results', '.xlsx',
                **kwargs)
        #write results
        xlwriter = pd.ExcelWriter(fn, engine = 'xlsxwriter')
        for xc in self:
            xc.fields.to_excel(xlwriter, sheet_name = xc.name)
        print('Full SectionBook results written to: "%s"' % fn)

    def update(self):
        """Executes all of the update functions"""
        self._update_fields()
        self._update_ROW_edge_max()
        self._update_tag_groups()

    def ion(self):
        emf_plots.ion()

    def show(self):
        emf_plots.show()

    #----------------------------------
    #functions that update SectionBook variables when CrossSections are done
    #being added or when CrossSection data changes

    def _update_fields(self):
        """run the fields calculations for each CrossSection in the
        SectionBook"""
        for xc in self:
            xc.calculate_fields()

    def _update_ROW_edge_max(self):
        """Execution populates the self.ROW_edge_max DataFrame with
        the most current results of the fields calculation in each
        CrossSection."""
        #gather ROW edge results
        L = len(self.xcs)
        El,Er,Bl,Br = np.zeros((L,)),np.zeros((L,)),np.zeros((L,)),np.zeros((L,))
        titles = []
        for i in range(L):
            xc = self.i(i)
            Bl[i] = xc.fields['Bmax'].iat[xc.lROWi]
            Br[i] = xc.fields['Bmax'].iat[xc.rROWi]
            El[i] = xc.fields['Emax'].iat[xc.lROWi]
            Er[i] = xc.fields['Emax'].iat[xc.rROWi]
            titles.append(xc.title)
        #construct DataFrame
        self.ROW_edge_max = pd.DataFrame(data = {
            'name': self.names, 'title': titles, 'Bmaxl': Bl, 'Emaxl': El,
            'Bmaxr': Br, 'Emaxr': Er}).sort_values('name')

    def _update_tag_groups(self):
        """Generate a list of lists of CrossSection indices with the same tag"""
        u = list(set(self.tags)) #get unique CrossSection tags
        self.tag_groups = [[] for i in range(len(u))]
        for i in range(len(self.xcs)):
            self.tag_groups[u.index(self.xcs[i].tag)].append(i)
