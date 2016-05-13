import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path

import emf_class
import emf_plots
import emf_calcs

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
    #create a SectionBook object to store the CrossSection objects
    xcs = emf_class.SectionBook(path.basename(file_path[:file_path.index('.')]))
    #convert the dataframes into a list of CrossSection objects
    titles = []
    for k in sheets.keys():
        #load miscellaneous information applicable for the whole CrossSection
        df = sheets[k]
        xc = emf_class.CrossSection(k)
        misc = df[1]
        xc.title = misc[0]
        xc.tag = misc[1]
        #check for duplicate title inputs
        if(xc.title in titles):
            raise(EMFError('Cross-sections should have unique Main Title entries. The Main Title "%s" in the sheet "%s" is used by at least one other sheet.' % (xc.title, k)))
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
            cond.tag = df[2].iloc[i]
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
            cond = emf_class.Conductor()
            cond.tag = df[11].iloc[i]
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
        #add the CrossSection object to the SectionBook
        xcs.add_section(xc)
    #update the SectionBook's remaining variables
    xcs.update()
    #return the SectionBook object
    return(xcs)

def path_manage(filename_if_needed, extension, **kwargs):
    """This function takes a path string through the kwarg 'path' and
    returns a path string with a file name at it's end, to save a file
    at that location. If the path string is a directory (ends with a slash),
    a new path string is returned with the 'filename_if_needed' and
    'extension' arguments appended. If the path string already has a file
    name at its end, the input extension will replace any preexisting one.
    If no path string is passed in via the keyword argument 'path', the
    returned path is simply the input filename_if_needed with the input
    extension at its end. Path strings without slashes at the end are
    assumed to be directories. If the last part of the path string is
    supposed to be the filename, include an extension."""
    #remove extensions from filename_if_needed
    if('.' in filename_if_needed):
        filename_if_needed = filename_if_needed[:filename_if_needed.index('.')]
    #make sure the extension has a period at its beginning
    if(extension):
        if(extension[0] != '.'):
            extension = '.' + extension
    #construct the path
    if('path' in kwargs.keys()):
        p = kwargs['path']
        #if there's a filename_if_needed argument and a 'path' keyword, assume
        #the 'path' argument is a directory and append a slash if it has no
        #leading director(y/ies)
        if(filename_if_needed and p):
            if(all(path.split(p))):
                p += '/'
        #split the path
        head,tail = path.split(p)
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

def run(template_path, output_path):
    """Import the templates in an excel file with the path 'template_path'
    then generate a workbook of all fields results and accompanying plots,
    saved to the directory described by 'output_path'. Finer control of the
    output, like x-distance cutoffs for the plots, is given up by the use
    of this function but it's a fast way to generate all the results. Pass
    an empty string to output_path to send results to the active directory."""
    #import templates
    b = import_template(template_path)
    #export the full results workbook
    b.full_export(path = output_path)
    #export ROW edge results
    b.ROW_edge_export(path = output_path)
    #export plots
    for xc in b:
        fig = emf_plots.plot_max_fields(xc, save = True, path = output_path)
        plt.close(fig)
