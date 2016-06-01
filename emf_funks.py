import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path

import emf_class
import emf_plots

def phasors_to_magnitudes(Ph_x, Ph_y):
    """Convert vectors of complex x and y phasors into real quantities,
    namely the amplitude of the field in the x and y directions, the
    product (the hypotenuse of the amplitudes), and the maximum field.
    args:
        Ph_x - numpy vector or number, phasor x component
        Ph_y - numpy vector or number, phasor y component"""
    #amplitude along each component, storing squared magnitudes for later
    mag_x_sq = np.real(Ph_x)**2 + np.imag(Ph_x)**2
    mag_x = np.sqrt(mag_x_sq)
    mag_y_sq = np.real(Ph_y)**2 + np.imag(Ph_y)**2
    mag_y = np.sqrt(mag_y_sq)
    #phase angle of each component
    phase_x = np.arctan(np.imag(Ph_x)/np.real(Ph_x))
    phase_y = np.arctan(np.imag(Ph_y)/np.real(Ph_y))
    #"product"
    prod = np.sqrt(mag_x**2 + mag_y**2)
    #maximum resultant value found by setting the time derivative of the
    #squared resultant magnitude to zero (Appendix 8.1 EPRI's "Big Red Book")
    num = mag_x_sq*np.sin(2*phase_x) + mag_y_sq*np.sin(2*phase_y)
    den = mag_x_sq*np.cos(2*phase_x) + mag_y_sq*np.cos(2*phase_y)
    t1 = (0.5)*np.arctan(-num/den)
    t2 = t1 + np.pi/2
    term1 = mag_x_sq*(np.cos(t1 + phase_x))**2
    term2 = mag_y_sq*(np.cos(t1 + phase_y))**2
    ax_mag1 = np.sqrt(term1 + term2)
    term1 = mag_x_sq*(np.cos(t2 + phase_x))**2
    term2 = mag_y_sq*(np.cos(t2 + phase_y))**2
    ax_mag2 = np.sqrt(term1 + term2)
    #pick out the semi-major axis magnitude from the two semi-axis results
    maxx = np.zeros((len(Ph_x),))
    for i in range(len(maxx)):
        if(ax_mag1[i] > ax_mag2[i]):
            maxx[i] = ax_mag1[i]
        else:
            maxx[i] = ax_mag2[i]
    #return the 4 output colums
    return(mag_x, mag_y, prod, maxx)

def load_template(file_path, **kwargs):
    """Import conductor data from an excel template, loading each conductor
    into a Conductor object, each Conductor into a CrossSection object, and
    each CrossSection object into a SectionBook object. The SectionBook
    object is returned.
    args:
        template_path - path to cross section template excel workbook
    kwargs:
        sheets - a list of sheet names to load, default is all sheets"""
    #import the cross sections as a dictionary of pandas dataframes, also
    #getting a list of the ordered sheets
    xl = pd.ExcelFile(file_path)
    sheets = xl.sheet_names
    frames = xl.parse(sheetname = None, skiprows = [0,1,2,3], parse_cols = 16,
                    header = None)
    #remove necessary sheets if the 'sheets' keyword is passed in
    if('sheets' in kwargs.keys()):
        include = kwargs['sheets']
        sheets = [sh for sh in sheets if sh in include]
    #create a SectionBook object to store the CrossSection objects
    xcs = emf_class.SectionBook(path.basename(file_path[:file_path.index('.')]))
    #convert the dataframes into a list of CrossSection objects
    titles = []
    for k in sheets:
        #load miscellaneous information applicable for the whole CrossSection
        df = frames[k]
        xc = emf_class.CrossSection(k)
        misc = df[1]
        xc.title = misc[0]
        xc.tag = misc[1]
        #check for duplicate title inputs
        if(xc.title in titles):
            raise(emf_class.EMFError('Cross-sections should have unique Main Title entries. The Main Title "%s" in the sheet "%s" is used by at least one other sheet.' % (xc.title, k)))
        else:
            titles.append(xc.title)
        xc.subtitle = misc[2]
        xc.soil_resistivity = misc[4]
        xc.max_dist = misc[5]
        xc.step = misc[6]
        xc.sample_height = misc[7]
        xc.lROW = misc[8]
        xc.rROW = misc[9]
        #load hot conductors
        for i in range(df[3].dropna().shape[0]):
            cond = emf_class.Conductor()
            cond.tag = df[2].iat[i]
            cond.freq = misc[2]
            cond.x = df[3].iat[i]
            cond.y = df[4].iat[i]
            cond.subconds = df[5].iat[i]
            cond.d_cond = df[6].iat[i]
            cond.d_bund = df[7].iat[i]
            cond.V = df[8].iat[i]
            cond.I = df[9].iat[i]
            cond.phase = df[10].iat[i]
            xc.hot.append(cond)
        #load grounded conductors
        for i in range(df[12].dropna().shape[0]):
            cond = emf_class.Conductor()
            cond.tag = df[11].iat[i]
            cond.freq = misc[2]
            cond.x = df[12].iat[i]
            cond.y = df[13].iat[i]
            cond.subconds = 1.
            cond.d_cond = df[14].iat[i]
            cond.d_bund = df[14].iat[i]
            cond.V = 0.
            cond.I = 0.
            cond.phase = 0.
            xc.gnd.append(cond)
        #calculate electric and magnetic fields automatically
        xc.calculate_fields()
        #add the CrossSection object to the SectionBook
        xcs.add_section(xc)
    #update the SectionBook's remaining variables
    xcs.update()
    #return the SectionBook object
    return(xcs)

def path_manage(filename_if_needed, extension, **kwargs):
    """This function takes a path string through the kwarg 'path' and
    returns a path string with a file name at it's end, to save a file
    at that location. If the path string is a directory (ends with a slash
    or doesn't but has preceding directory elements) a new path string is
    returned with the 'filename_if_needed' and 'extension' arguments
    appended. If the path string already has a file name at its end, the
    input extension will replace any preexisting one. If no path string is
    passed in via the keyword argument 'path', the returned path is simply
    the input filename_if_needed with the input extension at its end. If the
    last part of the path string is supposed to be the filename, include an
    extension.
    args:
        filename_if_needed - name of file if not in path string
        extension - file extension
    kwargs:
        path - string, destination/filename for saved file(s)"""
    #remove extensions from filename_if_needed
    if('.' in filename_if_needed):
        filename_if_needed = filename_if_needed[:filename_if_needed.index('.')]
    #make sure the extension has a period at its beginning
    if(extension):
        if(extension[0] != '.'):
            extension = '.' + extension
    #construct the path
    if('path' in kwargs.keys()):
        p = path.normcase(kwargs['path'])
        #if there's a filename_if_needed argument and a 'path' keyword, assume
        #the 'path' argument is a directory. Append a slash if it has no
        #leading director(y/ies)
        if(filename_if_needed and p):
            if(all(path.split(p))):
                p += '/'
        #split the path
        head,tail = path.split(p)
        #check that head describes an existing directory if it isn't empty
        if(head and (not path.isdir(head))):
            raise(emf_class.EMFError('"%s" was not recognized as an existing directory, invalid path string' % head))
        #if a file name lies at the end of p, replace its extension
        if(tail):
            if('.' in tail):
                tail = tail[:tail.index('.')]
            if(head):
                return(head + '/' + tail + extension)
            else:
                return(tail + extension)
        #if no file name, but a directory
        elif(head):
            return(head + '/' + filename_if_needed + extension)
    #if 'path' kwarg is empty or missing
    return(filename_if_needed + extension)

def run(template_path, **kwargs):
    """Import the templates in an excel file with the path 'template_path'
    then generate a workbook of all fields results and accompanying plots.
    Use the 'path' keyword argument to specify a destination for the output,
    otherwise it will be save to the active directory, Finer control of the
    output, like x-distance cutoffs for the plots, is given up by the use of
    this function but it's a fast way to generate all the results. Returns
    a SectionBook object of the imported template.
    args:
        template_path - path to cross section template excel workbook
    kwargs:
        sheets - a list of sheet names to load, default is all sheets
        path - string, destination/filename for saved files
        format - string, saved plot format (usually 'png' or 'pdf')
        xmax - cutoff distance from ROW center in plots"""
    #force saving for the plotting functions if there is no 'path' keyword
    kwargs['save'] = True
    #import templates
    b = load_template(template_path)
    #export the full results workbook
    b.full_export(**kwargs)
    #export ROW edge results
    b.ROW_edge_export(**kwargs)
    #export single CrossSection plots
    for xc in b:
        fig = emf_plots.plot_max_fields(xc, **kwargs)
        plt.close(fig)
    #export group comparison plots
    emf_plots.plot_groups(b, **kwargs)
    return(b)
