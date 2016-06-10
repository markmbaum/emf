import os
import copy
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import emf_class
import emf_calcs
import emf_plots

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
    file_path = check_extention(file_path, 'xlsx', """
        Templates must be excel workbooks. The input target path
            "%s"
        is not recognized as an excel file""" % file_path)
    xl = pd.ExcelFile(file_path)
    sheets = xl.sheet_names
    frames = xl.parse(sheetname = None, skiprows = [0,1,2,3], parse_cols = 16,
                    header = None)
    #remove necessary sheets if the 'sheets' keyword is passed in
    if('sheets' in kwargs.keys()):
        include = kwargs['sheets']
        sheets = [sh for sh in sheets if sh in include]
    #create a SectionBook object to store the CrossSection objects
    xcs = emf_class.SectionBook(os.path.basename(file_path[:file_path.index('.')]))
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
            raise(emf_class.EMFError("""
            Cross-sections should have unique Main Title entries.
            Main Title:
                "%s"
            in the sheet
                "%s"
            is used by at least one other sheet.""" % (xc.title, k)))
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

def optimize_phasing(xc):
    """Permute the phasing of the non-grounded conductors and find the
    arrangement that results in the lowest fields at the left and right
    edge of the ROW. The number of hot conductors must be a multiple of
    three. The phases of consecutive groups of three conductors are
    swapped around, assuming that those groups represent a single
    three-phase transfer circuit."""
    #number of hot wires
    N = len(xc.hot)
    #check the number of hot lines
    if(N % 3 != 0):
        raise(EMFError("""
        The number of hot (not grounded) conductors must be a multiple
        of three for phase optimization."""))
    #number of circuits, groups of 3 hot conductors
    G = N/3
    #all permutations of the phases of each 3 line circuit
    perm = list(itertools.permutations([0,1,2]))
    #all possible arrangements of line phasings, 6 permutations for each circuit
    #so 6^(N/3) total line arrangements
    P = list(itertools.product(perm, repeat = G))
    #flatten the elements of P
    for i in range(len(P)):
        P[i] = [item for sublist in P[i] for item in sublist]
    #turn it into an array and add multiples of 3 to the appropriate columns
    P = np.array(P, dtype = int)
    for i in range(3, P.shape[1], 3):
        P[:,i:i+3] += i
    #variables to find the minima with respect to each field and ROW edge
    B_left_min, B_left_idx = np.inf, -1
    B_right_min, B_right_idx = np.inf, -1
    E_left_min, E_left_idx = np.inf, -1
    E_right_min, E_right_idx = np.inf, -1
    #loop through all possible combinations in P
    x_ROW = np.array([xc.x_sample[xc.lROWi], xc.x_sample[xc.rROWi]])
    y_ROW = np.array([xc.y_sample[xc.lROWi], xc.y_sample[xc.rROWi]])
    phase = np.empty((N,))
    for i in range(len(P)):
        #get a row of indices from P
        idx = P[i,:]
        #calculate fields with index swapped data
        Ex, Ey = emf_calcs.E_field(xc.x[:N], xc.y[:N], xc.subconds[:N],
            xc.d_cond[:N], xc.d_bund[:N], xc.V[:N], xc.phase[idx],
            x_ROW, y_ROW)
        Ex, Ey, Eprod, Emax = emf_calcs.phasors_to_magnitudes(Ex, Ey)
        Bx, By = emf_calcs.B_field(xc.x[:N], xc.y[:N], xc.I[:N],
            xc.phase[idx], x_ROW, y_ROW)
        Bx, By, Bprod, Bmax = emf_calcs.phasors_to_magnitudes(Bx, By)
        #test for minima
        if(Bmax[0] < B_left_min):
            B_left_min = Bmax[0]
            B_left_idx = i
        if(Bmax[1] < B_right_min):
            B_right_min = Bmax[1]
            B_right_idx = i
        if(Emax[0] < E_left_min):
            E_left_min = Emax[0]
            E_left_idx = i
        if(Emax[1] < E_right_min):
            E_right_min = Emax[1]
            E_right_idx = i
    #return results in a DataFrame
    results = pd.DataFrame(data = {
        'Conductor Tag' : [c.tag for c in xc.hot],
        'Optimal Phasing - Bmax Left ROW Edge' : xc.phase[P[B_left_idx,:]],
        'Optimal Phasing - Bmax Right ROW Edge' : xc.phase[P[B_right_idx,:]],
        'Optimal Phasing - Emax Left ROW Edge' : xc.phase[P[E_left_idx,:]],
        'Optimal Phasing - Emax Right ROW Edge' : xc.phase[P[E_right_idx,:]]})

    return(results)

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
    if(not ('path' in kwargs)):
        kwargs['save'] = True
        #also direct output files to the same directory as the template
        kwargs['path'] = os.path.dirname(template_path)
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
        p = os.path.normcase(kwargs['path'])
        #if there's a filename_if_needed argument and a 'path' keyword, assume
        #the 'path' argument is a directory. Append a slash if it has no
        #leading director(y/ies)
        if(filename_if_needed and p):
            if(all(os.path.split(p))):
                p += '/'
        #split the path
        head,tail = os.path.split(p)
        #check that head describes an existing directory if it isn't empty
        if(head and (not os.path.isdir(head))):
            raise(emf_class.EMFError("""
            "%s" was not recognized as an existing directory,
            invalid path string""" % head))
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

def check_extention(file_path, correct_ext, message):
    """Check that a file path has the desired extention, raising an error if
    not and appending the correct extension if no extension is present.
    args:
        file_path - a target file path
        correct_ext - the correct extension for the target path
        message - error message if the extention is wrong
    returns:
        file_path"""
    if(correct_ext[0] == '.'):
        correct_ext = correct_ext[1:]
    if('.' in file_path):
        ext = file_path[file_path.index('.')+1:]
        if(not (correct_ext == ext)):
            raise(emf_class.EMFError(message))
    else:
        file_path += '.' + correct_ext
    return(file_path)
