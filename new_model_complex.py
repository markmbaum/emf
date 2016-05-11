import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path

class FLDError(Exception):
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return(self.message)

class Conductor:

    def __init__(self):
        self.name = None #conductor label
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
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            s += str(k) + ': ' + str(v[k]) + '\n'
        return(s)

class CrossSection:
    """Store the input information for a power line cross section and the
    results of electric and magnetic field simulations across the section."""

    def __init__(self, name):
        self.name = name #mandatory
        self.main_title = ''
        self.subtitle = ''
        self.soil_resistivity = 100. #?
        self.max_dist = None #maximum simulated distance from the ROW center
        self.step = None #step size for calculations
        self.sample_height = 3. #uniform sample height
        self.lROW = None #exact coordinate of the left ROW edge
        self.lROWi = None #self.fields index closest to self.lROW
        self.rROW = None #exact coordinate of the left ROW edge
        self.rROWi = None #self.fields index closest to self.rROW
        self.hot = [] #list of Conductor objects with nonzero voltage
        self.gnd = [] #list of Conductor objects with zero voltage
        self.fields = pd.DataFrame(columns = ['Bx','By','Bprod','Bmax',
                                            'Ex','Ey','Eprod','Emax'])
        self.B_color = 'darkgreen' #for plotting
        self.E_color = 'midnightblue' #for plotting

    def __str__(self):
        v = vars(self)
        keys = v.keys()
        s = '\n'
        for k in keys:
            if(k != 'fields'):
                s += str(k) + ': ' + str(v[k]) + '\n'
        s += '\ninspect self.fields separately to see field simulation results\n'
        return(s)

    def calculate_fields(self):
        """Calculate electric and magnetic fields across the ROW"""
        #calculate sample points
        N = 1 + 2*self.max_dist/self.step
        x = np.linspace(-self.max_dist, self.max_dist, num = N)
        y = self.sample_height*np.ones((N,))

        #assemble all the conductor data in arrays for calculations
        conds = self.hot + self.gnd
        x_c = np.array([c.x for c in conds])
        y_c = np.array([c.y for c in conds])
        subc = np.array([c.subconds for c in conds])
        d_c = np.array([c.d_cond for c in conds])
        d_b = np.array([c.d_bund for c in conds])
        V = np.array([c.V for c in conds])
        I = np.array([c.I for c in conds])
        ph = np.array([c.phase for c in conds])

        #calculate electric field
        Ex, Ey = E_field(x_c, y_c, subc, d_c, d_b, V, ph, x, y)
        Ex, Ey, Eprod, Emax = phasors_to_magnitudes(Ex, Ey)

        #calculate magnetic field
        Bx, By = B_field(x_c, y_c, I, ph, x, y)
        Bx, By, Bprod, Bmax = phasors_to_magnitudes(Bx, By)

        #store the values
        self.fields = pd.DataFrame({'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
                                    'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax},
                                    index = x)

        #update ROW edge index variables
        #if ROW edge lies between two sample points, use the one closer to zero
        d = np.absolute((self.fields.index - self.lROW).values)
        self.lROWi = max(self.fields.index[d == np.min(d)])

        #if ROW edge lies between two sample points, use the one closer to zero
        d = np.absolute((self.fields.index - self.rROW).values)
        self.rROWi = min(self.fields.index[d == np.min(d)])

        #return the fields dataframe
        return(self.fields)

    def optimize_phasing(self):
        """Permute the phasing of the non-grounded conductors and find the
        arrangement that results in the lowest fields at the left and right
        edge of the ROW. The number of hot conductors must be a multiple of
        three. The phases of consecutive groups of three conductors are
        swapped around, assuming that those groups represent a single
        three-phase transfer line."""
        #check the number of hot lines
        if(self.N_hot % 3 != 0):
            raise(FLDError("""The number of hot (not grounded) conductors must
                            be a multiple of three."""))
        #number of 3 phase groups
        G = self.N_hot/3
        #number of permutations
        N = 6*G
        #all permutations of a single three phase group
        perm = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]
        #variables to store results of permutations
        #
        #...to be continued?

    def prepare_fig(self, **kwargs):
        """Snippet executed at the beginning of each plotting method"""
        #prepare figure and axis
        plt.rc('font', family = 'calibri')
        k = kwargs
        keys = k.keys()
        if('figure' in keys):
            fig = k['figure']
        else:
            fig = plt.figure()
        ax = plt.gca()
        #get x cutoff, if any
        if('xmax' in keys):
            xmax = abs(k['xmax'])
        else:
            xmax = max(abs(self.fields.index))
        return(fig, ax, xmax)

    def save_fig(self, fig, **kwargs):
        """Snippet executed at the end of each plotting method"""
        #save the fig, or don't
        k = kwargs
        keys = k.keys()
        #force saving if a path is passed in
        if('path' in keys):
            k['save'] = True
        #condition filename and format strings for saving
        if('save' in keys):
            if(k['save']):
                #get filename
                fn = path_manage(self.name, '', **kwargs)
                #get format/extension
                if('format' in keys):
                    fmt = k['format']
                    if('.' in fmt):
                        fmt = fmt[fmt.index('.')+1:]
                else:
                    fmt = 'png'
                #save the plot
                fn += '.' + fmt
                plt.savefig(fn, format = fmt)
                print('plot saved to "%s"' % fn)

    def plot_Emax(self, **kwargs):
        """Plot the maximum electric field along the ROW with conductor
        locations shown in artificial but to-scale locations. Pass in an
        existing figure with keyword argument 'figure' to recycle an object.
        Pass in a plot title with the keyword argument 'title' to specify an
        exact title, otherwise the main_title will be used. Use the kwarg
        'xmax' to cut the plotted fields at a certain distance from the ROW
        center. If the keyword argument 'save' is passed in True, the plot
        will be saved. Use the keyword argument 'path' to specify the path
        or filename of the saved plot. The format can also be specified
        (usually 'png' or 'pdf') with the 'format' keyword."""
        #keyword variables
        k = kwargs
        keys = k.keys()
        #get axes and x cutoff
        (fig, ax, xmax) = self.prepare_fig(**kwargs)
        #plot the field curve
        ax.plot(self.fields['Emax'][-xmax:xmax], '.-', color = self.E_color)
        #plot wires
        y = [c.y for c in self.hot + self.gnd]
        scale = (.25*max(self.fields['Emax'])/min(y))
        (x_cond,y_cond) = [c.x for c in self.hot],[scale*c.y for c in self.hot]
        ax.plot(x_cond, y_cond, 'kd')#hot lines
        (x_cond,y_cond) = [c.x for c in self.gnd],[scale*c.y for c in self.gnd]
        ax.plot(x_cond, y_cond, 'd', color = 'gray')#ground lines
        #plot ROW lines and adjust ylimits to make room for legend
        ax.set_ylim([0, max(self.fields['Emax'])*1.35])
        yl = ax.get_ylim()
        ax.plot([self.lROW]*2, yl, 'k--', [self.rROW]*2, yl, 'k--')
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax.set_ylabel('Maximum Electric Field (kV/m)', fontsize = 14)
        if('title' in keys):
            t = k['title']
        else:
            t = '%s, Maximum Electric Field' % self.main_title
        ax.set_title(t)
        ax.legend(['Electric Field (kV/m)','Conductors','Grounded Conductors',
                    'ROW Edge'], numpoints = 1, fontsize = 12)
        #save the fig, or don't
        self.save_fig(fig, **kwargs)
        #return
        return(fig)

    def plot_Bmax(self, **kwargs):
        """Plot the maximum magnetic field along the ROW with conductor
        locations shown in artificial but to-scale locations. Pass in an
        existing figure with keyword argument 'figure' to recycle an object.
        Pass in a plot title with the keyword argument 'title' to specify an
        exact title, otherwise the main_title will be used. Use the kwarg
        'xmax' to cut the plotted fields at a certain distance from the ROW
        center. If the keyword argument 'save' is passed in, the plot
        will be saved. Use the keyword argument 'path' to specify the path
        or filename of the saved plot. The format can also be specified
        (usually 'png' or 'pdf') with the 'format' keyword."""
        #keyword variables
        k = kwargs
        keys = k.keys()
        #get axes and x cutoff
        (fig, ax, xmax) = self.prepare_fig(**kwargs)
        #plot the field curve
        ax.plot(self.fields['Bmax'][-xmax:xmax], '.-', color = self.B_color)
        #plot wires
        y = [c.y for c in self.hot + self.gnd]
        scale = (.25*max(self.fields['Bmax'])/min(y))
        (x_cond,y_cond) = [c.x for c in self.hot],[scale*c.y for c in self.hot]
        ax.plot(x_cond, y_cond, 'kd')#hot lines
        (x_cond,y_cond) = [c.x for c in self.gnd],[scale*c.y for c in self.gnd]
        ax.plot(x_cond, y_cond, 'd', color = 'gray')#ground lines
        #plot ROW lines and adjust ylimits to make room for legend
        ax.set_ylim([0, max(self.fields['Bmax'])*1.35])
        yl = ax.get_ylim()
        ax.plot([self.lROW]*2, yl, 'k--', [self.rROW]*2, yl, 'k--')
        #set axis text and legend
        ax.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax.set_ylabel('Maximum Magnetic Field (mG)', fontsize = 14)
        if('title' in keys):
            t = k['title']
        else:
            t = '%s, Maximum Magnetic Field' % self.main_title
        ax.set_title(t)
        ax.legend(['Magnetic Field (mG)','Conductors','Grounded Conductors',
                    'ROW Edge'], numpoints = 1, fontsize = 12)
        #save the fig, or don't
        self.save_fig(fig, **kwargs)
        #return
        return(fig)

    def plot_max_fields(self, **kwargs):
        """Plot the maximum fields along the ROW with conductor
        locations shown in artificial but to-scale locations. Pass in an
        existing figure with keyword argument 'figure' to recycle an object.
        Pass in a plot title with the keyword argument 'title' to specify an
        exact title, otherwise the main_title will be used. Use the kwarg
        'xmax' to cut the plotted fields at a certain distance from the ROW
        center. If the keyword argument 'save' is passed in, the plot
        will be saved. Use the keyword argument 'path' to specify the path
        or filename of the saved plot. The format can also be specified
        (usually 'png' or 'pdf') with the 'format' keyword."""
        #keyword variables
        k = kwargs
        keys = k.keys()
        #get axes and x cutoff
        (fig, ax_B, xmax) = self.prepare_fig(**kwargs)
        ax_E = ax_B.twinx()
        #plot the field curves
        hB, = ax_B.plot(self.fields['Bmax'][-xmax:xmax], '.-', color = self.B_color)
        hE, = ax_E.plot(self.fields['Emax'][-xmax:xmax], '.-', color = self.E_color)
        #plot wires
        y = [c.y for c in self.hot + self.gnd]
        scale = (.25*max(self.fields['Bmax'])/min(y))
        (x_cond,y_cond) = [c.x for c in self.hot],[scale*c.y for c in self.hot]
        hhot, = ax_B.plot(x_cond, y_cond, 'kd')#hot lines
        (x_cond,y_cond) = [c.x for c in self.gnd],[scale*c.y for c in self.gnd]
        hgnd, = ax_B.plot(x_cond, y_cond, 'd', color = 'gray')#ground lines
        #plot ROW lines and adjust ylimits to make room for legend
        ax_B.set_ylim([0, max(self.fields['Bmax'])*1.4])
        ax_E.set_ylim([0, max(self.fields['Emax'])*1.4])
        yl = ax_B.get_ylim()
        hROW = ax_B.plot([self.lROW]*2, yl, 'k--', [self.rROW]*2, yl, 'k--')
        #set axis text
        ax_B.set_xlabel('Distance from Center of ROW (ft)', fontsize = 14)
        ax_B.set_ylabel('Maximum Magnetic Field (mG)',
                        fontsize = 14, color = self.B_color)
        ax_E.set_ylabel('Maximum Electric Field (kV/m)',
                        fontsize = 14, color = self.E_color)
        if('title' in keys):
            t = k['title']
        else:
            t = '%s, Maximum Magnetic and Electric Fields' % self.main_title
        ax_B.set_title(t)
        #set color of axis spines and ticklabels
        ax_B.spines['left'].set_color(self.B_color)
        ax_B.spines['right'].set_color(self.B_color)
        ax_E.spines['left'].set_color(self.E_color)
        ax_E.spines['right'].set_color(self.E_color)
        ax_B.tick_params(axis = 'y', colors = self.B_color)
        ax_E.tick_params(axis = 'y', colors = self.E_color)
        #legend
        ax_B.legend([hB, hE, hhot, hgnd, hROW[0]],
                    ['Magnetic Field (mG)','Electric Field (kV/m)','Conductors',
                    'Grounded Conductors','ROW Edge'],
                    numpoints = 1, fontsize = 12)
        #save the fig, or don't
        self.save_fig(fig, **kwargs)
        #return
        return(fig)

    def compare_DAT(self, DAT_path, **kwargs):
        """Load a FIELDS output file, find the absolute and percentage
        differences between it and the CrossSection objects results, and
        write them to an excel file. The default excel file name is the
        CrossSection's main_title with '-DAT_comparison' appended to it. Use
        the keyword argument 'path' to specify a different one."""
        #load the .DAT file into a dataframe
        df = pd.read_table(DAT_path, skiprows = [0,1,2,3,4,5,6],
                            delim_whitespace = True, header = None,
                            names = ['Bx', 'By', 'Bprod', 'Bmax',
                                    'Ex', 'Ey', 'Eprod', 'Emax'],
                            index_col = 0)
        #prepare a dictionary to create a Panel
        comp = {'FIELDS_output' : df,
                'New_model_output' : self.fields,
                'Absolute Difference' : self.fields - df,
                'Percent Difference' : 100*(self.fields - df)/self.fields}
        #write the frames to a spreadsheet
        fn = path_manage(self.name + '-DAT_comparison', '.xlsx', **kwargs)
        pd.Panel(data = comp).to_excel(fn, index_label = 'x')
        print('DAT comparison book saved to: %s' % fn)

class SectionBook:

    def __init__(self, name):
        self.name = name
        self.xcs = [] #list of cross section objects
        self.xcname2idx = dict() #mapping dictionary for CrossSection retrieval
        self.xcnames = [] #list of CrossSection names

    def __iter__(self):
        for xc in self.xcs:
            yield(xc)

    def __getitem__(self, key):
        try:
            idx = self.xcname2idx[key]
        except(KeyError):
            return(False)
        else:
            return(self.xcs[idx])

    def i(self, idx):
        return(self.xcs[idx])

    def add_section(self, xc):
        if(xc.name in self.xcnames):
            raise(FLDError('CrossSection name "%s" already exists.' % xc.name))
        else:
            self.xcname2idx[xc.name] = len(self.xcs)
            self.xcs.append(xc)
            self.xcnames.append(xc.name)

    def export(self, **kwargs):
        """Write all of the cross sections to an excel workbook"""
        #path management
        fn = path_manage(self.name + '-export', '.xlsx', **kwargs)
        #data management
        data = dict(zip(self.xcnames, [xc.fields for xc in self.xcs]))
        pd.Panel(data = data).to_excel(fn, index_label = 'x')
        print('SectionBook written to "%s"' % fn)

def E_field(x_cond, y_cond, subconds, d_cond, d_bund, V_cond, p_cond, x, y):
    """Calculate the approximate electric field generated by a group of
    conductors. Each of the inputs labeled '_cond' should be an numpy
    array of parameters, where each index in those arrays describes a
    unique conductor, i.e. the 0th value in each variable is attributed to
    one power line."""

    #convenient variables/constants
    epsilon = 8.854e-12
    C = 1./(2.*np.pi*epsilon)
    N = len(p_cond)     #number of conductors
    Z = len(x)          #number of sample points or x,y pairs

    #calculate the effective conductor diameters
    d_cond  = d_bund*((subconds*d_cond/d_bund)**(1./subconds))

    #conversions
    x_cond = x_cond*0.3048          #convert to meters
    y_cond = y_cond*0.3048          #convert to meters
    d_cond = d_cond*0.0254          #convert to meters
    V_cond = V_cond/np.sqrt(3)      #convert to ground reference from line-line
                                    #reference, leave in units of kV
    p_cond = p_cond*2*np.pi/360.    #convert to radians
    x      = x*0.3048               #convert to meters
    y      = y*0.3048               #convert to meters

    #compute the matrix of potential coefficients
    P = np.empty((N,N))
    #diagonals
    for i in range(N):
        P[i,i] = C*np.log(4*y_cond[i]/d_cond[i])
    #other elements
    for a in range(N):
        for b in range(N):
            if(a != b):
                n = (x_cond[a] - x_cond[b])**2 + (y_cond[a] + y_cond[b])**2
                d = (x_cond[a] - x_cond[b])**2 + (y_cond[a] - y_cond[b])**2
                P[a,b] = C*np.log(np.sqrt(n/d))

    #initialize complex voltage phasors
    V = V_cond*(np.cos(p_cond) + complex(1j)*np.sin(p_cond))

    #compute real and imaginary charge phasors
    Q = np.linalg.solve(P, V)

    #compute components of the electric field phasors at each point, with each
    #column of E_x and E_y storing components due to each conductor, so that the
    #rows represent a spatial point across the ROW or an x,y pair
    Ex = np.empty((N,Z))
    Ey = np.empty((N,Z))
    #first compute the coefficients without the charges
    for a in range(Z):
        #denominators, squared distance between the point and the conductors
        d1 = (x[a] - x_cond)**2 + (y[a] - y_cond)**2
        d2 = (x[a] - x_cond)**2 + (y[a] + y_cond)**2
        #x component numerator, the same for the conductor and its image
        nx = C*(x[a] - x_cond)
        #y component numerators, different for the conductor and its image
        ny1 = C*(y[a] - y_cond)
        ny2 = C*(y[a] + y_cond)
        #evaluate
        Ex[:,a] = nx/d1 - nx/d2
        Ey[:,a] = ny1/d1 - ny2/d2

    #multiply the charges by the field coefficients calculated above
    Q = np.tile(np.reshape(Q, (N, 1)), (1,Z))
    Ex = Ex*Q
    Ey = Ey*Q

    #sum the real parts and imaginary parts of each phasor for each sample
    #point, yielding the sum of phasors for each conductor, which are the final
    #phasors for each point
    Ex = np.sum(Ex, axis = 0)
    Ey = np.sum(Ey, axis = 0)

    #return phasors, complex numbers, for the x and y components
    return(Ex, Ey)

def B_field(x_cond, y_cond, I_cond, p_cond, x, y):
    """Calculate the approximate magnetic field generated by a group of
    conductors. Each of the variables labeled '_cond' should be an numpy
    array of parameters, where each index in those arrays describes a
    unique conductor, i.e. the 0th value in each variable is attributed to
    one power line."""

    #convenient variables/constants
    mu = 4*np.pi*1e-7       #magnetic permeability constant, in SI units
    C = 1e7*mu/(2*np.pi)    #constant of convenience, converted for milligauss
    N = len(p_cond)         #number of conductors
    Z = len(x)              #number of sample points or x,y pairs

    #conversions
    x_cond = x_cond*0.3048          #convert to meters
    y_cond = y_cond*0.3048          #convert to meters
    p_cond = p_cond*2*np.pi/360.    #convert to radians
    x      = x*0.3048               #convert to meters
    y      = y*0.3048               #convert to meters

    #initialize complex current phasors
    I = I_cond*(np.cos(p_cond) + complex(1j)*np.sin(p_cond))

    #compute magnetic field component phasors for each x,y point
    Bx = np.zeros((Z,), dtype = complex)
    By = np.zeros((Z,), dtype = complex)
    for a in range(Z): #x,y pairs
        for b in range(N): #conductors
            dx = x[a] - x_cond[b]
            dy = y[a] - y_cond[b]
            #magnitude phasor
            B = C*I[b]/np.sqrt(dx**2 + dy**2)
            #break it up into components
            theta = np.arctan(abs(dy/dx))
            #x component calculated with sine and y component with cosine
            #because the field is perpendicular to the line to the conductor,
            Bx[a] -= np.sign(dy)*np.sin(theta)*B
            By[a] += np.sign(dx)*np.cos(theta)*B

    #return phasors, complex numbers, for the x and y components
    return(Bx, By)

def phasors_to_magnitudes(Ph_x, Ph_y):
    """Convert vectors of complex x and y phasors into real quantities,
    namely the amplitude of the field in the x and y directions, the
    product (the hypotenuse of the amplitudes), and the maximum field."""
    #amplitude along each component
    mag_x = np.sqrt(np.real(Ph_x)**2 + np.imag(Ph_x)**2)
    mag_y = np.sqrt(np.real(Ph_y)**2 + np.imag(Ph_y)**2)
    #"product"
    prod = np.sqrt(mag_x**2 + mag_y**2)
    #"max" - attempts at non-brute-force methods have slightly missed the mark,
    #for mysterious reasons
    t = np.linspace(0,2*np.pi,1001)
    maxi = np.zeros((len(Ph_x),))
    for i in range(len(maxi)):
        Rx = mag_x[i]*np.cos(t + np.arctan(np.imag(Ph_x[i])/np.real(Ph_x[i])))
        Ry = mag_y[i]*np.cos(t + np.arctan(np.imag(Ph_y[i])/np.real(Ph_y[i])))
        maxi[i] = max(np.sqrt(Rx**2 + Ry**2))

    #R = np.sqrt(Ph_x**2 + Ph_y**2)
    #R = np.sqrt((np.real(Ph_x) + np.real(Ph_y))**2 + (np.imag(Ph_x) + np.imag(Ph_y))**2)
    #maxi = np.abs(R)

    return(mag_x, mag_y, prod, maxi)

def import_template(file_path):
    """Import conductor data from an excel template, loading each conductor
    into a Conductor object, each Conductor into a CrossSection object, and
    each CrossSection object into a SectionBook object. The SectionBook
    object is returned."""
    #import the cross sections as a dictionary of pandas dataframes
    sheets = pd.read_excel(file_path, sheetname = None, skiprows = [0,1,2,3],
                        parse_cols = 16, header = None)
    #create a SectioBook object to store the CrossSection objects
    xcs = SectionBook(path.basename(file_path[:file_path.index('.')]))
    #convert the dataframes into a list of CrossSection objects
    for k in sheets.keys():
        #load miscellaneous information applicable for the whole CrossSection
        df = sheets[k]
        xc = CrossSection(k)
        misc = df[1]
        xc.main_title = misc[0]
        xc.subtitle = misc[1]
        xc.soil_resistivity = misc[3]
        xc.max_dist = misc[4]
        xc.step = misc[5]
        xc.sample_height = misc[6]
        xc.lROW = misc[7]
        xc.rROW = misc[8]
        #load hot conductors
        for i in range(df[3].dropna().shape[0]):
            cond = Conductor()
            cond.name = df[2].iloc[i]
            cond.freq = misc[2]
            cond.x = df[3].iloc[i]
            cond.y = df[4].iloc[i]
            cond.subconds = df[5].iloc[i]
            cond.d_cond = df[6].iloc[i]
            cond.d_bund = df[7].iloc[i]
            cond.V = df[8].iloc[i]
            cond.I = df[9].iloc[i]
            cond.phase = df[10].iloc[i]
            xc.hot.append(cond)
        #load grounded conductors
        for i in range(df[12].dropna().shape[0]):
            cond = Conductor()
            cond.name = df[11].iloc[i]
            cond.freq = misc[2]
            cond.x = df[12].iloc[i]
            cond.y = df[13].iloc[i]
            cond.subconds = 1.
            cond.d_cond = df[14].iloc[i]
            cond.d_bund = df[14].iloc[i]
            cond.V = 0.
            cond.I = 0.
            cond.phase = 0.
            xc.gnd.append(cond)
        #calculate electric and magnetic fields automatically
        xc.calculate_fields()
        #replace the dataframe with the CrossSection object
        xcs.add_section(xc)
    #return the list of CrossSection objects
    return(xcs)

def path_manage(filename_if_needed, extension, **kwargs):
    """Expects the keyword argument 'path' which leads to a path string,
    otherwise the return will be 'filename_if_needed + '.' + extension'. If
    the path string is a directory, the string should end with a slash."""
    if(extension):
        if(extension[0] != '.'):
            extension = '.' + extension
    if('path' in kwargs.keys()):
        path_string = kwargs['path']
        head,tail = path.split(path_string)
        if(tail):
            if('.' in tail):
                tail = tail[:tail.index('.')]
            fn = head + '/' + tail + extension
        elif(head):
            fn = head + '/' + filename_if_needed + extension
        else:
            fn = filename_if_needed + extension
    else:
        fn = filename_if_needed + extension
    return(fn)

def run(template_path, output_path):
    """Import the templates in an excel file with the path 'template_path'
    then generate a workbook of all fields results with accompanying plot,
    saved to the directory described by 'output_path'. Fine control of the
    output, like x-distance cutoffs for the plots, is given up by the use
    of this function but it's a fast way to generate all the results. Pass
    an empty string to output_path to send results to the active directory."""
    #import templates
    b = import_template(template_path)
    #condition output path
    if(output_path):
        if((output_path[-1] != '/') and (output_path != '\\')):
            output_path += '/'
    #export the results workbook
    b.export(path = output_path)
    #export plots
    for xc in b:
        xc.plot_max_fields(save = True, path = output_path)

run('XC-template.xlsx', '')
