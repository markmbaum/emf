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
    #import the cross sections as a dictionary of pandas DataFrames, also
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
    if('sheets' in kwargs):
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
        raise(emf_class.EMFError("""
        The number of hot (not grounded) conductors must be a multiple
        of three for phase optimization. Circuits are assumed to be
        three-phase and conductors comprising each circuit are assumed to
        be consecutive groups of three, in the order that they appear in the
        template."""))
    #number of circuits, groups of 3 hot conductors
    G = int(N/3)
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
    print('Optimizing phasing of CrossSection %s, with %d hot conductors'
        % (xc.name, len(xc.hot)))
    print('Number of permutations to test = 6^(%d/3) = %d'
        % (len(xc.hot), P.shape[0]))
    #variables to find the minima with respect to each field and ROW edge
    B_left_min, B_left_idx = np.inf, -1
    B_right_min, B_right_idx = np.inf, -1
    E_left_min, E_left_idx = np.inf, -1
    E_right_min, E_right_idx = np.inf, -1
    #loop through all possible combinations in P
    x_ROW = np.array([xc.x_sample[xc.lROWi], xc.x_sample[xc.rROWi]])
    y_ROW = np.array([xc.y_sample[xc.lROWi], xc.y_sample[xc.rROWi]])
    #array for swapping phases, zeros in the grounded slots
    phase = np.zeros((N + len(xc.gnd),))
    for i in range(len(P)):
        #get a row of indices from P
        idx = P[i,:]
        #swap phases accordingly
        phase[:N] = xc.phase[idx]
        #calculate fields with index swapped data
        Ex, Ey = emf_calcs.E_field(xc.x, xc.y, xc.subconds, xc.d_cond,
            xc.d_bund, xc.V, phase, x_ROW, y_ROW)
        Ex, Ey, Eprod, Emax = emf_calcs.phasors_to_magnitudes(Ex, Ey)
        Bx, By = emf_calcs.B_field(xc.x, xc.y, xc.I, phase, x_ROW, y_ROW)
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

def target_fields(xc, B_l, B_r, E_l, E_r, hot, gnd):
    """Increase conductor y coordinates until fields at ROW edges are below
    thresholds. All selected conductors are adjusted by the same amount.
    If any of the thresholds are empty or false, None is returned for their
    adjustment result.
    args:
        B_l - magnetic field threshold at left ROW edge
        B_r - magnetic field threshold at right ROW edge
        E_l - electric field threshold at left ROW edge
        E_r - electric field threshold at right ROW edge
        hot - indices of hot conductors in self.hot to raise, accepts 'all'
        gnd - indices of ground conductors in self.gnd to raise, accepts 'all'
    returns:
        h_B_l - height adjustment necessary for left magnetic field
        h_B_r - height adjustment necessary for right magnetic field
        h_E_l - height adjustment necessary for left electric field
        h_E_r - height adjustment necessary for right electric field"""
    if(hot == 'all'):
        hot = list(range(len(xc.hot)))
    if(gnd == 'all'):
        gnd = list(range(len(xc.gnd)))
    #maximum number of iterations, lots of breathing room
    max_iter = 1e4
    hlow = 0.0
    hhigh = 1.0e6
    #flattened indices
    conds = np.array(list(hot) + [len(xc.hot) + i for i in gnd])
    #run secant method to find adjustments for each target
    h_B_l, h_B_r, h_E_l, h_E_r = None, None, None, None
    if(B_l):
        h_B_l = bisect(xc, conds, xc.lROWi, B_funk, B_l, hlow, hhigh, max_iter)
    if(B_r):
        h_B_r = bisect(xc, conds, xc.rROWi, B_funk, B_r, hlow, hhigh, max_iter)
    if(E_l):
        h_E_l = bisect(xc, conds, xc.lROWi, E_funk, E_l, hlow, hhigh, max_iter)
    if(E_r):
        h_E_r = bisect(xc, conds, xc.rROWi, E_funk, E_r, hlow, hhigh, max_iter)
    return(h_B_l, h_B_r, h_E_l, h_E_r)

def bisect(xc, conds, sample_idx, funk, target, hlow, hhigh, max_iter):
    #evaluate at the bracketing values
    flow = funk(hlow, target, xc, conds, sample_idx)
    fhigh = funk(hhigh, target, xc, conds, sample_idx)
    #check that the root is bracketed
    if(flow*fhigh > 0.):
        print("""
        The root is not bracketed. Rootfinding with bisection will fail.
            f(h_0 = %g) = %g
            f(h_1 = %g) = %g
        Returning 'out of range'""" % (hlow, flow, hhigh, fhigh))
        return('out of range')
    #evaluate at a midpoint
    hmid = (hhigh + hlow)/2.0
    fmid = funk(hmid, target, xc, conds, sample_idx)
    count = 1
    #iterate
    while((abs((hhigh - hlow)/target) > 1.0e-9) and (count < max_iter)):
        #test and throw half out
        if(fmid*flow > 0.):
            hlow = hmid
        elif(fmid*fhigh > 0.):
            hhigh = hmid
        #evaluate at middle
        hmid = (hhigh + hlow)/2.0
        fmid = funk(hmid, target, xc, conds, sample_idx)
        #increment
        count += 1
    if(count == max_iter):
        raise(emf_class.EMFError("""
        Divergence in bisection method. Iteration limit of %d was exceeded.
        Likely cause is a target field value that is too low. Line adjustments
        would have to be extreme or """
        % max_iter))
    return(hmid)

def B_funk(h, target, xc, conds, sample_idx):
    i = sample_idx
    #apply height adjustment to selected conductors
    y = xc.y.copy()
    y[conds] = y[conds] + h
    #calculate B field at ROW edge
    Bx, By = emf_calcs.B_field(xc.x, y, xc.I, xc.phase,
        xc.x_sample[i:i+1], xc.y_sample[i:i+1])
    Bx, By, Bprod, Bmax = emf_calcs.phasors_to_magnitudes(Bx, By)
    return(Bmax[0] - target)

def E_funk(h, target, xc, conds, sample_idx):
    i = sample_idx
    #apply height adjustment to selected conductors
    y = xc.y.copy()
    y[conds] = y[conds] + h
    #calculate E field at ROW edge
    Ex, Ey = emf_calcs.E_field(xc.x, y, xc.subconds, xc.d_cond,
        xc.d_bund, xc.V, xc.phase, xc.x_sample[i:i+1], xc.y_sample[i:i+1])
    Ex, Ey, Eprod, Emax = emf_calcs.phasors_to_magnitudes(Ex, Ey)
    return(Emax[0] - target)

def run(template_path, **kwargs):
    """Import the templates in an excel file with the path 'template_path'
    then generate a workbook of all fields results and lots of plots.
    Use the 'path' keyword argument to specify a destination for the output,
    otherwise it will be saved to the template's directory. Returns
    a SectionBook object.
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
    sb = load_template(template_path)
    #export the full results workbook
    sb.full_export(**kwargs)
    #export ROW edge results
    sb.ROW_edge_export(**kwargs)
    #export single CrossSection plots
    for xc in b:
        fig = emf_plots.plot_max_fields(xc, **kwargs)
        plt.close(fig)
    #export group comparison plots
    emf_plots.plot_groups(sb, **kwargs)
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
    if('path' in kwargs):
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
