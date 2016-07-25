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
            Main Title: "%s"
            in sheet: "%s"
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
        tags, x, y = [], [], []
        for i in range(df[3].dropna().shape[0]):
            cond = emf_class.Conductor()
            cond.tag = df[2].iat[i]
            #check for conductors with identical tags (names/labels)
            if(cond.tag in tags):
                raise(emf_class.EMFError("""
                Conductors in a Cross Section must have unique tags.
                The conductor tag "%s" in sheet:
                    "%s"
                is used at least twice."""
                % (cond.tag, k)))
            else:
                tags.append(cond.tag)
            cond.freq = misc[2]
            cond.x = df[3].iat[i]
            cond.y = df[4].iat[i]
            #check for conductors with identical x,y coordinates
            if(cond.x in x):
                idx = x.index(cond.x)
                if(cond.y == y[idx]):
                    raise(emf_class.EMFError("""
                Conductors cannot have identical x,y coordinates. Conductor
                "%s" is in the exact same place as conductor "%s"."""
                % (cond.tag, tags[idx])))
            else:
                x.append(cond.x)
                y.append(cond.y)
            cond.subconds = df[5].iat[i]
            cond.d_cond = df[6].iat[i]
            cond.d_bund = df[7].iat[i]
            cond.V = df[8].iat[i]
            cond.I = df[9].iat[i]
            cond.phase = df[10].iat[i]
            xc.hot.append(cond)
        #load grounded conductors
        tags, x, y = [], [], []
        for i in range(df[12].dropna().shape[0]):
            cond = emf_class.Conductor()
            cond.tag = df[11].iat[i]
            #check for conductors with identical tags (names/labels)
            if(cond.tag in tags):
                raise(emf_class.EMFError("""
                Conductors in a Cross Section must have unique tags.
                The conductor tag "%s" in sheet:
                    "%s"
                is used at least twice."""
                % (cond.tag, k)))
            else:
                tags.append(cond.tag)
            cond.freq = misc[2]
            cond.x = df[12].iat[i]
            cond.y = df[13].iat[i]
            #check for conductors with identical x,y coordinates
            if(cond.x in x):
                idx = x.index(cond.x)
                if(cond.y == y[idx]):
                    raise(emf_class.EMFError("""
                Conductors cannot have identical x,y coordinates. Conductor
                "%s" is in the exact same place as conductor "%s"."""
                % (cond.tag, tags[idx])))
            else:
                x.append(cond.x)
                y.append(cond.y)
            cond.subconds = 1.
            cond.d_cond = df[14].iat[i]
            cond.d_bund = df[14].iat[i]
            cond.V = 0.
            cond.I = 0.
            cond.phase = 0.
            xc.gnd.append(cond)
        #add the CrossSection object to the SectionBook
        xcs.add_section(xc)
    #update the SectionBook, which initiates fields calculations and
    #population of lots of other variables in the CrossSection and SectionBook
    #objects
    xcs.update()
    #return the SectionBook object
    return(xcs)

def optimize_phasing(xc, circuits, **kwargs):
    """Permute the phasing of non-grounded conductors and find the
    arrangement that results in the lowest fields at the left and right
    edge of the ROW. The number of hot conductors must be a multiple of
    three. The phases of consecutive groups of three conductors are
    swapped around, assuming that those groups represent a single
    three-phase transfer circuit.
    args:
        xc - target CrossSection object
        circuits - a list of lists, or 'all'. If a list of lists, each
                   sublist contains the integer indices of the conductors
                   that belong to a circuit, indexed from zero. If 'all',
                   circuits are assumed to be consecutive groups of three
                   conductors.
    kwargs:
        save - toggle saving of the results DataFrame to an excel book
        path - location/filename for saved results workbook, forces saving
               even if no 'save' keyword is used.
    returns:
        results - DataFrame listing conductor phasings that optimize
                  electric and magnetic fields at both ROW edges.
        opt - a new SectionBook object containing the permuted phasings that
             optimize the E and B fields at the left and right ROW edges."""

    if(circuits == 'all'):
        #number of hot wires
        N = len(xc.hot)
        #check the number of hot lines
        if(N % 3 != 0):
            raise(emf_class.EMFError("""
            The number of hot (not grounded) conductors must be a multiple
            of three for phase optimization with 'all' circuits. Circuits are
            assumed to be three-phase and conductors comprising each circuit
            are assumed to be consecutive groups of three, in the order that
            they appear in the template. The number of hot conductors is not a
            multiple of three in CrossSection "%s" """ % xc.name))
        #number of circuits, groups of 3 hot conductors
        G = int(N/3)
        #circuits, consecutive groups of three conductors
        circuits = [range(i*3,i*3 + 3) for i in range(G)]

    #all permutations of the phases of each circuit
    perm = []
    for c in circuits:
        perm.append(list(itertools.permutations(c)))
    #all possible arrangements of line phasings, 6 permutations for each circuit
    #so 6^(N/3) total line arrangements
    P = list(itertools.product(*perm))
    #flatten the elements of P
    for i in range(len(P)):
        P[i] = [item for sublist in P[i] for item in sublist]
    #turn P into an array
    P = np.array(P, dtype = int)
    #notifications, remove later
    print("""
    Optimizing phasing of CrossSection "%s":
        Permuting %d total conductors in %d circuits
        Number of permutations to test at ROW edges = %d"""
        % (xc.title, P.shape[1], len(circuits), P.shape[0]))
    #variables to find the minima with respect to each field and ROW edge
    B_left_min, B_left_idx, B_right_min, B_right_idx = np.inf, -1, np.inf, -1
    E_left_min, E_left_idx, E_right_min, E_right_idx = np.inf, -1, np.inf, -1
    #make sure the necessary CrossSection variables are set
    xc.update_arrays()
    #get coordinates of the ROW edges
    x_ROW = np.array([xc.x_sample[xc.lROWi], xc.x_sample[xc.rROWi]])
    y_ROW = np.array([xc.y_sample[xc.lROWi], xc.y_sample[xc.rROWi]])
    #array for swapping phases, zeros in the grounded slots
    phasing = xc.phase.copy()
    #store a flattened version of the conductor indices for swapping
    conds = np.array([item for sublist in circuits for item in sublist],
        dtype = int)
    #loop through all possible combinations in P
    for i in range(len(P)):
        #swap phases accordingly
        phasing[conds] = xc.phase[P[i,:]]
        #calculate fields with index swapped data
        Ex, Ey = emf_calcs.E_field(xc.x, xc.y, xc.subconds, xc.d_cond,
                                    xc.d_bund, xc.V, phasing, x_ROW, y_ROW)
        Ex, Ey, Eprod, Emax = emf_calcs.phasors_to_magnitudes(Ex, Ey)
        Bx, By = emf_calcs.B_field(xc.x, xc.y, xc.I, phasing, x_ROW, y_ROW)
        Bx, By, Bprod, Bmax = emf_calcs.phasors_to_magnitudes(Bx, By)
        #test for minima
        if(Bmax[0] < B_left_min):
            B_left_min, B_left_idx = Bmax[0], i
        if(Bmax[1] < B_right_min):
            B_right_min, B_right_idx = Bmax[1], i
        if(Emax[0] < E_left_min):
            E_left_min, E_left_idx = Emax[0], i
        if(Emax[1] < E_right_min):
            E_right_min, E_right_idx = Emax[1], i
    #return results in a DataFrame
    results = pd.DataFrame(data = {
        'Optimal Phasing - Bmax Left ROW Edge' : xc.phase[P[B_left_idx,:]],
        'Optimal Phasing - Bmax Right ROW Edge' : xc.phase[P[B_right_idx,:]],
        'Optimal Phasing - Emax Left ROW Edge' : xc.phase[P[E_left_idx,:]],
        'Optimal Phasing - Emax Right ROW Edge' : xc.phase[P[E_right_idx,:]]},
        index = [xc.hot[i].tag for i in conds])
    #compile a new sectionbook with the optimal phasings
    opt = emf_class.SectionBook('%s-optimal_phasing' % xc.title)
    names = ['Optimized_for_Bmax_left','Optimized_for_Bmax_right',
            'Optimized_for_Emax_left','Optimized_for_Emax_right']
    titles = ['Bmax_l','Bmax_r','Emax_l','Emax_r']
    subtitles = results.columns
    for n, t, s in zip(names, titles, subtitles):
        #copy the input XC
        new_xc = copy.deepcopy(xc)
        #change the identification fields
        new_xc.name, new_xc.title, new_xc.subtitle = n, t, s
        #swap the conductor phasings
        for c in new_xc.hot:
            t = c.tag
            if(t in results.index):
                c.phase = results.at[t, s]
        #store new_xc in the SectionBook
        opt.add_section(new_xc)
    #update everything in the SectionBook
    opt.update()
    #deal with saving
    if('path' in kwargs):
        kwargs['save'] = True
    if('save' in kwargs):
        if(kwargs['save']):
            fn = path_manage(xc.name + '_phase_optimization', 'xlsx', **kwargs)
            xl = pd.ExcelWriter(fn, engine = 'xlsxwriter')
            results.to_excel(xl, index_label = 'Conductor Tag',
                sheet_name = 'phase_assignments')
            opt.ROW_edge_export(xl = xl)
            for xc in opt:
                xc.fields.to_excel(xl, sheet_name = xc.name)
            xl.save()
            print('Optimal phasing results written to "%s"' % fn)

    return(results, opt)

def target_fields(xc, hot, gnd, B_l, B_r, E_l, E_r, **kwargs):
    """Increase conductor y coordinates until fields at ROW edges are below
    thresholds. All selected conductors are adjusted by the same amount.
    If any of the thresholds are empty or false, None is returned for their
    adjustment result.
    args:
        xc - CrossSection object to perform adjustments on
        hot - indices of hot conductors in self.hot to raise, accepts 'all'
        gnd - indices of ground conductors in self.gnd to raise,
              accepts 'all'
        B_l - magnetic field threshold at left ROW edge*
        B_r - magnetic field threshold at right ROW edge*
        E_l - electric field threshold at left ROW edge*
        E_r - electric field threshold at right ROW edge*

            *an implicitly False input will ignore that field-edge
             combination, return None in the return variable 'h', and
             cause the returned SectionBook to omit that field-edge combo.

    kwargs:
        max_iter - maximum number of bisection iterations allowed
                   default is 1e3
        rel_err - tolerance threshold for relative error (e.g. 0.01 is 1 %)
                  default is 1e-6.
        hhigh - upper limit of the height adjustment, default is 1.0e6
        save - toggle saving of the results DataFrame to an excel book
        path - location/filename for saved results workbook, forces saving
               even if no 'save' keyword is used.
    returns:
        h - height adjustments necessary for E and B fields at left and
            right ROW edges. The ordering is:
                    (B_left, B_right, E_left, E_right)
        adj - a new SectionBook object with the adjusted conductor heights
             for each scenario in a CrossSection"""
    #convert 'all' inputs to numeric indices
    if(hot == 'all'):
        hot = list(range(len(xc.hot)))
    if(gnd == 'all'):
        gnd = list(range(len(xc.gnd)))
    #maximum number of iterations and relative error tolerance
    if('max_iter' in kwargs):
        max_iter = kwargs['max_iter']
    else:
        max_iter = 1e3
    if('rel_err' in kwargs):
        rel_err = kwargs['rel_err']
    else:
        rel_err = 1.0e-6
    if('hhigh' in kwargs):
        hhigh = kwargs['hhigh']
    else:
        hhigh = 1.0e6
    hlow = 0.0
    #flattened indices
    conds = np.array(list(hot) + [len(xc.hot) + i for i in gnd])
    #make sure the necessary CrossSection variables are set
    xc.update_arrays()
    #run secant method to find adjustments for each target
    h_B_l, h_B_r, h_E_l, h_E_r = None, None, None, None
    if(B_l):
        h_B_l = bisect(xc, conds, xc.lROWi, B_funk, B_l, hlow, hhigh,
            max_iter, rel_err)
    if(B_r):
        h_B_r = bisect(xc, conds, xc.rROWi, B_funk, B_r, hlow, hhigh,
            max_iter, rel_err)
    if(E_l):
        h_E_l = bisect(xc, conds, xc.lROWi, E_funk, E_l, hlow, hhigh,
            max_iter, rel_err)
    if(E_r):
        h_E_r = bisect(xc, conds, xc.rROWi, E_funk, E_r, hlow, hhigh,
            max_iter, rel_err)
    #create return variables
    h = (h_B_l, h_B_r, h_E_l, h_E_r)
    adj = emf_class.SectionBook('%s-height_adjusted' % xc.title)
    names = ['Adjusted_for_Bmax_left','Adjusted_for_Bmax_right',
            'Adjusted_for_Emax_left','Adjusted_for_Emax_right']
    titles = ['Bmax_l','Bmax_r','Emax_l','Emax_r']
    subtitles = ['Height Adjusted for %f mG at left ROW edge' % B_l,
                'Height Adjusted for %f mG at left ROW edge' % B_r,
                'Height Adjusted for %f kV/m at left ROW edge' % E_l,
                'Height Adjusted for %f kV/m at left ROW edge' % E_r]
    for n, t, s, a in zip(names, titles, subtitles, h):
        if(a != None):
            #copy the input XC
            new_xc = copy.deepcopy(xc)
            #change the identification fields
            new_xc.name, new_xc.title, new_xc.subtitle = n, t, s
            #adjust conductor heights
            for idx in hot:
                new_xc.hot[idx].y += a
            for idx in gnd:
                new_xc.gnd[idx].y += a
            #store new_xc in the SectionBook
            adj.add_section(new_xc)
    #update everythin in the SectionBook
    adj.update()
    #deal with saving
    if('path' in kwargs):
        kwargs['save'] = True
    if('save' in kwargs):
        if(kwargs['save']):
            fn = path_manage(xc.name + '_height_adjustments', 'xlsx', **kwargs)
            xl = pd.ExcelWriter(fn, engine = 'xlsxwriter')
            pd.DataFrame(data = list(h), index = names,
                columns = ['Height Addition (ft)']).to_excel(xl,
                sheet_name = 'Adjustments', index_label = 'Field - ROW Edge')
            adj.ROW_edge_export(xl = xl)
            for xc in adj:
                xc.fields.to_excel(xl, sheet_name = xc.name)
            xl.save()
            print('Optimal phasing results written to "%s"' % fn)

    return(h, adj)

def bisect(xc, conds, sample_idx, funk, target, hlow, hhigh, max_iter, rel_err):
    #get sample x and y arrays with a single element in each
    x_sample = np.array([xc.x_sample[sample_idx]], dtype = float)
    y_sample = np.array([xc.y_sample[sample_idx]], dtype = float)
    #evaluate at the bracketing values
    flow = funk(hlow, target, xc, conds, x_sample, y_sample)
    fhigh = funk(hhigh, target, xc, conds, x_sample, y_sample)
    #check that the root is bracketed
    if(flow*fhigh > 0.):
        raise(emf_class.EMFError("""
        The root is not bracketed with an upper height adjustment limit
        of %g. Rootfinding with bisection can't be performed.
            f(h_0 = %g) = %g
            f(h_1 = %g) = %g""" % (hhigh, hlow, flow, hhigh, fhigh)))
    #evaluate at a midpoint
    hmid = (hhigh + hlow)/2.0
    fmid = funk(hmid, target, xc, conds, x_sample, y_sample)
    count = 1
    #iterate
    while((abs(fmid/target) > rel_err) and (count < max_iter)):
        #test and throw half out
        if(fmid*flow > 0.):
            hlow = hmid
        elif(fmid*fhigh > 0.):
            hhigh = hmid
        elif(fmid == 0.):
            return(hmid)
        #evaluate at middle
        hmid = (hhigh + hlow)/2.0
        fmid = funk(hmid, target, xc, conds, x_sample, y_sample)
        #increment
        count += 1
    #check if the iteration limit was hit
    if(count == max_iter):
        raise(emf_class.EMFError("""
        Failure in bisection method. The iteration limit of %d was exceeded
        with a relative error threshold of %g. The final estimate was
        %g""" % (max_iter, rel_err, fmid)))
    return(hmid)

def B_funk(h, target, xc, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xc.y.astype(float, copy = True)
    y[conds] += h
    #calculate B field at ROW edge
    Bx, By = emf_calcs.B_field(xc.x, y, xc.I, xc.phase,
        x_sample, y_sample)
    Bx, By, Bprod, Bmax = emf_calcs.phasors_to_magnitudes(Bx, By)
    return(Bmax[0] - target)

def E_funk(h, target, xc, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xc.y.astype(float, copy = True)
    y[conds] += h
    #calculate E field at ROW edge
    Ex, Ey = emf_calcs.E_field(xc.x, y, xc.subconds, xc.d_cond,
        xc.d_bund, xc.V, xc.phase, x_sample, y_sample)
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
    sb.results_export(**kwargs)
    #export ROW edge results
    sb.ROW_edge_export(**kwargs)
    #export single CrossSection plots
    for xc in sb:
        fig = emf_plots.plot_max_fields(xc, **kwargs)
        plt.close(fig)
    #export group comparison plots
    emf_plots.plot_groups(sb, **kwargs)
    return(sb)

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
            "%s"
            was not recognized as an existing directory, invalid path string"""
            % head))
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
        correct_ext - the correct extension for the target path, with or
                      without the period
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

def is_number(s):
    """Check if an element can be converted to a float, returning `True`
    if it can and `False` if it can't"""
    if(s == None):
        return(False)
    else:
        try:
            float(s)
        except ValueError:
            return(False)
        else:
            return(True)
