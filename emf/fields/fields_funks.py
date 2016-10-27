from .. import os, np, pd, shutil, itertools

from ..emf_funks import (_path_manage, _check_extension, _is_number,
                        _check_intable, _flatten, _sig_figs,
                        _path_str_condition)

import fields_class
import fields_calcs
import fields_plots

def drop_template(*args, **kw):
    """Copy the emf.fields template in the current directory or a
    directory specified by an input string
    args:
        drop_path = string, path of copied template file"""
    #check inputs
    if(len(args) > 1):
        raise(fields_class.EMFError("""
        drop_template only accepts zero or one input argument.
        A string can be passed to specify the directory in which the
        template file is copied. With no arguments, the template file
        is copied into the current directory."""))
    elif(len(args) == 1):
        kw = {'path': args[0]}
    #get template file path
    template_path = os.path.dirname(os.path.dirname(__file__))
    template_path = os.path.join(template_path, 'templates')
    template_path = os.path.join(template_path, 'fields-template.xlsx')
    #get drop path
    drop_path = _path_manage('fields-template', 'xlsx', **kw)
    #copy and notify
    shutil.copyfile(template_path, drop_path)
    print('emf.fields template written to: %s' % drop_path)

def run(template_path, **kw):
    """Import the templates in an excel file with the path 'template_path'
    then generate a workbook of all fields results and lots of plots.
    Use the 'path' keyword argument to specify a destination for the output,
    otherwise it will be saved to the template's directory. Returns
    a SectionBook object.
    args:
        template_path - path to cross section template excel workbook
    kw:
        sheets - a list of sheet names to load, default is all sheets
        path - string, destination/filename for saved files
        format - string, saved plot format (usually 'png' or 'pdf')
        xmax - cutoff distance from ROW center in plots"""
    #force saving for the plotting functions if there is no 'path' keyword
    if(not ('path' in kw)):
        kw['save'] = True
        #also direct output files to the same directory as the template
        kw['path'] = os.path.dirname(template_path)
    #import templates
    sb = load_template(template_path)
    #export the full results workbook
    sb.results_export(**kw)
    #export ROW edge results
    sb.ROW_edge_export(**kw)
    #export single CrossSection plots
    for xs in sb:
        fig, ax_E, ax_B = fields_plots.plot_max_fields(xs, **kw)
        fields_plots.close(fig)
    #export group comparison line plots
    fields_plots.plot_groups(sb, **kw)
    #export group ROW comparison bar plots
    fields_plots.plot_groups_at_ROW(sb, **kw)
    return(sb)

def load_template(file_path, **kw):
    """Import conductor data from an excel template, loading each conductor
    into a Conductor object, each Conductor into a CrossSection object, and
    each CrossSection object into a SectionBook object. The SectionBook
    object is returned.
    args:
        template_path - string, path to cross section template excel
                        workbook
    kw:
        sheets - list of strings, a list of sheet names to load, default is
                 all sheets"""
    #import the cross sections as a dictionary of pandas DataFrames, also
    #getting a list of the ordered sheets
    file_path = _check_extension(file_path, 'xlsx', """
        Templates must be excel workbooks. The input target path
            "%s"
        is not recognized as an excel file""" % file_path)
    xl = pd.ExcelFile(file_path)
    sheets = xl.sheet_names
    frames = xl.parse(sheetname = None, skiprows = [0,1,2,3],
            parse_cols = 16, header = None)
    #remove necessary sheets if the 'sheets' keyword is passed in
    if('sheets' in kw):
        include = kw['sheets']
        sheets = [sh for sh in sheets if sh in include]
    #create a SectionBook object to store the CrossSection objects
    basename = os.path.basename(file_path)
    if('.' in basename):
        name = basename[:basename.index('.')]
    else:
        name = basename
    sb = fields_class.SectionBook(name)
    #convert the dataframes into a list of CrossSection objects
    titles = []
    for k in sheets:
        #load miscellaneous information applicable to the whole CrossSection
        df = frames[k]
        xs = fields_class.CrossSection(k)
        misc = df[1].values
        xs.tag = misc[0]
        xs.title = str(misc[1])
        #check for duplicate title inputs
        if(xs.title in titles):
            raise(fields_class.EMFError("""
            Cross-sections should have unique title entries.
            title: "%s"
            in sheet: "%s"
            is used by at least one other sheet.""" % (xs.title, k)))
        else:
            titles.append(xs.title)
        xs.soil_resistivity = float(misc[3])
        xs.max_dist = float(misc[4])
        xs.step = float(misc[5])
        xs.sample_height = float(misc[6])
        xs.lROW = float(misc[7])
        xs.rROW = float(misc[8])
        #load hot conductors
        tags, x, y = [], [], []
        for i in range(df[3].dropna().shape[0]):
            #initialize a Conductor
            cond = fields_class.Conductor(df[2].iat[i])
            #check for conductors with identical tags (names/labels)
            if(cond.tag in tags):
                raise(fields_class.EMFError("""
                Conductors in a Cross Section must have unique tags.
                The conductor tag "%s" in sheet:
                    "%s"
                is used at least twice."""
                % (cond.tag, k)))
            else:
                tags.append(cond.tag)
            #cond.freq = misc[2]
            cond.x = df[3].iat[i]
            cond.y = df[4].iat[i]
            #check for conductors with identical x,y coordinates
            if(cond.x in x):
                idx = x.index(cond.x)
                if(cond.y == y[idx]):
                    raise(fields_class.EMFError("""
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
            xs.add_conductor(cond)
        #load grounded conductors
        tags, x, y = [], [], []
        for i in range(df[12].dropna().shape[0]):
            #initialize a Conductor
            cond = fields_class.Conductor(df[11].iat[i])
            #check for conductors with identical tags (names/labels)
            if(cond.tag in tags):
                raise(fields_class.EMFError("""
                Conductors in a Cross Section must have unique tags.
                The conductor tag "%s" in sheet:
                    "%s"
                is used at least twice."""
                % (cond.tag, k)))
            else:
                tags.append(cond.tag)
            #cond.freq = misc[2]
            cond.x = df[12].iat[i]
            cond.y = df[13].iat[i]
            #check for conductors with identical x,y coordinates
            if(cond.x in x):
                idx = x.index(cond.x)
                if(cond.y == y[idx]):
                    raise(fields_class.EMFError("""
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
            xs.add_conductor(cond)
        #add the CrossSection object to the SectionBook
        #fields automatically updated upon addition to SectionBook
        sb.add_section(xs)
    #return the SectionBook object
    return(sb)

def optimize_phasing(xs, circuits, **kw):
    """Permute the phasing of non-grounded conductors and find the
    arrangement that results in the lowest fields at the left and right
    edge of the ROW. The number of hot conductors must be a multiple of
    three. The phases of consecutive groups of three conductors are
    swapped around, assuming that those groups represent a single
    three-phase transfer circuit.
    args:
        xs - target CrossSection object
        circuits - list of lists of Conductor tags, or 'all'. If a list of
                   lists, each sublist contains the Conductor tags
                   of the Conductors comprising a single circuit.
                   If 'all', circuits are assumed to be consecutive groups
                   of three conductors. (consecutive according to the order
                   in which hot conductors were added to the CrossSection)
    kw:
        save - bool, toggle saving of the results DataFrame to an excel book
        path - string, location/filename for saved results workbook, forces
               saving even if no 'save' keyword is used.
    returns:
        res - pandas DataFrame listing conductor phasings that optimize
              electric and magnetic fields at both ROW edges.
        opt - new SectionBook object containing the permuted phasings that
              optimize the E and B fields at the left and right ROW edges."""

    if(circuits == 'all'):
        #number of hot wires
        N = len(xs.hot)
        #check the number of hot lines
        if(N % 3 != 0):
            raise(fields_class.EMFError("""
            The number of hot (not grounded) conductors must be a multiple
            of three for phase optimization with 'all' circuits. Circuits are
            assumed to be three-phase and conductors comprising each circuit
            are assumed to be consecutive groups of three, in the order that
            they appear in the template. The number of hot conductors is not a
            multiple of three in the CrossSection named: %s""" % xs.sheet))
        #number of circuits, groups of 3 hot conductors
        G = int(N/3)
        #circuits, consecutive groups of three conductors
        circuits = [range(i*3,i*3 + 3) for i in range(G)]
    else:
        #check that all conductor tags are present and refer to hot conds
        gnd = xs.gnd
        for circ in range(len(circuits)):
            for tag in circ:
                if(xs[tag] is None):
                    raise(EMFError("""
                    Unrecognized conductor tag: %s
                    All conductor tags must refer to Conductor objects present
                    in the target CrossSecton object.""" % repr(tag)))
                if(xs[tag] in gnd):
                    raise(EMFError("""
                    Only phasing of non-grounded Conductors can be permuted.
                    Tag "%s" refers to a grounded Conductor""" % repr(tag)))
        #convert the conductor tags to integer indices
        for i in range(len(circuits)):
            for j in range(len(circuits[i])):
                circuits[i][j] = xs.hot.index(xs[circuits[i][j]])
    #all permutations of the phases of each circuit
    perm = []
    for c in circuits:
        perm.append(list(itertools.permutations(c)))
    #all possible arrangements of line phasings, 6 permutations for each circuit
    #so 6^(N/3) total line arrangements. Leave P as a generator to avoid storing
    #a huge, factorial sized array of indices
    P = itertools.product(*perm)
    #variables to find the minima with respect to each field and ROW edge
    B_left_min, B_left_arr, B_right_min, B_right_arr = np.inf, [], np.inf, []
    E_left_min, E_left_arr, E_right_min, E_right_arr = np.inf, [], np.inf, []
    #get coordinates of the ROW edges
    x_ROW = np.array([xs.lROW, xs.rROW], dtype=float)
    y_ROW = xs.sample_height*np.ones((2,), dtype=float)
    #array for swapping phases, zeros in the grounded slots
    phasing = xs.phase.copy()
    #store a flattened version of the conductor indices for swapping
    conds = np.array([i for j in circuits for i in j], dtype=int)
    #loop through all possible arrangements in P
    for arr in P:
        #calculate fields at ROW edges with the new arrangement
        Bmax,Emax,new_arr = _phasing_test(xs, x_ROW, y_ROW, conds, phasing, arr)
        #test for minima
        if(Bmax[0] < B_left_min):
            B_left_min, B_left_arr = Bmax[0], new_arr
        if(Bmax[1] < B_right_min):
            B_right_min, B_right_arr = Bmax[1], new_arr
        if(Emax[0] < E_left_min):
            E_left_min, E_left_arr = Emax[0], new_arr
        if(Emax[1] < E_right_min):
            E_right_min, E_right_arr = Emax[1], new_arr
    #return results in a DataFrame
    results = pd.DataFrame(data={
        'Optimal Phasing - Bmax Left ROW Edge': xs.phase[B_left_arr],
        'Optimal Phasing - Bmax Right ROW Edge': xs.phase[B_right_arr],
        'Optimal Phasing - Emax Left ROW Edge': xs.phase[E_left_arr],
        'Optimal Phasing - Emax Right ROW Edge': xs.phase[E_right_arr]},
        index=[xs.hot[i].tag for i in conds])
    #compile a new sectionbook with the optimal phasings
    fn = _path_str_condition(xs.sheet).replace(' ', '-')
    xs = xs.copy()
    opt = fields_class.SectionBook(xs.sheet + '-optimal_phasing')
    xs.sheet += ' (original)'
    xs.tag = 'Phase Optimized'
    opt.add_section(xs)
    names = ['Optimized for Bmax left','Optimized for Bmax right',
            'Optimized for Emax left','Optimized for Emax right']
    tags = ['Phase Optimized']*4
    titles = results.columns
    for n, ti, ta in zip(names, titles, tags):
        #copy the input xs
        new_xs = xs.copy()
        #change the identification fields
        new_xs.sheet, new_xs.title, new_xs.tag = n, ti, ta
        #swap the conductor phasings
        for c in new_xs.hot:
            t = c.tag
            if(t in results.index):
                c.phase = results.at[t, ti]
        #store new_xs in the SectionBook
        opt.add_section(new_xs)
    #deal with saving
    if('path' in kw):
        kw['save'] = True
    if('save' in kw):
        if(kw['save']):
            fn = _path_manage(fn + '_phase_optimization', 'xlsx', **kw)
            xl = pd.ExcelWriter(fn, engine = 'xlsxwriter')
            results.to_excel(xl, index_label = 'Conductor Tag',
                sheet_name = 'phase_assignments')
            opt.ROW_edge_export(xl = xl)
            df, c, h = _xs_sb_diff(xs, opt)
            df.to_excel(xl, sheet_name = 'ROW_edge_diff', index = False,
                    columns = c, header = h)
            for xs in opt:
                xs.fields.to_excel(xl, sheet_name = xs.sheet)
            xl.save()
            print('Phase optimization results written to: %s' % fn)

    return(results, opt)

def _phasing_test(xs, x_ROW, y_ROW, conds, phasing, arr):
    """Calculate fields at the ROW edges for a given phasing arrangement,
    called by optimize_phasing()
    args:
        xs - CrossSection object with phases to test
        x_ROW - numpy array of the ROW edge x coordinates, to avoid repeated
                creation of this array for fields_calcs to work on
        y_ROW - numpy array of the ROW edge y coordinates, to avoid repeated
                creation of this array for fields_calcs to work on
        conds - array of ints, indices of conductors in xs.hot under
                consideration
        phasing - numpy array, a copy of xs.phase to mess with
        arr - numpy array, permuted phases for the conductors indexed by
              conds, list of lists that will be flattened
    returns:
        Bmax - array of two B field values, the 0th being the field at the
               left ROW edge, 1st is at right ROW edge
        Emax - array of two E field values, the 0th being the field at the
               left ROW edge, 1st is at right ROW edge
        new_arr - flattened version of arr in case it needs to be stored"""
    #flatten the new arrangement
    new_arr = _flatten(arr)
    #swap phases according to the new phasing arrangement
    phasing[conds] = xs.phase[new_arr]
    #calculate fields with index swapped phases
    Ex, Ey = fields_calcs.E_field(xs.x, xs.y, xs.subconds, xs.d_cond,
                                xs.d_bund, xs.V, phasing, x_ROW, y_ROW)
    Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
    Bx, By = fields_calcs.B_field(xs.x, xs.y, xs.I, phasing, x_ROW, y_ROW)
    Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
    #return results
    return(Bmax, Emax, new_arr)

def target_fields(xs, tags, B_l, B_r, E_l, E_r, **kw):
    """Increase conductor y coordinates until fields at ROW edges are below
    thresholds. All selected conductors are adjusted by the same amount.
    If any of the thresholds are empty or false, None is returned for their
    adjustment result.
    args:
        xs - CrossSection object to perform adjustments on
        tags - iterable of Conductor tags, identifying which ones to raise
        B_l - magnetic field threshold at left ROW edge*
        B_r - magnetic field threshold at right ROW edge*
        E_l - electric field threshold at left ROW edge*
        E_r - electric field threshold at right ROW edge*

            *an implicitly False input will ignore that field-edge
             combination, return None in the return variable 'h', and
             cause the returned SectionBook to omit that field-edge combo.

    kw:
        max_iter - maximum number of _bisection iterations allowed
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
    if(tags == 'all'):
        tags = xs.tags
    #maximum number of iterations and relative error tolerance
    if('max_iter' in kw):
        max_iter = kw['max_iter']
    else:
        max_iter = 1e3
    if('rel_err' in kw):
        rel_err = kw['rel_err']
    else:
        rel_err = 1.0e-6
    if('hhigh' in kw):
        hhigh = kw['hhigh']
    else:
        hhigh = 1.0e6
    hlow = 0.0
    #flattened indices
    temp = xs.tags
    conds = np.array([temp.index(i) for i in tags])
    #run secant method to find adjustments for each target
    h_B_l, h_B_r, h_E_l, h_E_r = None, None, None, None
    if(B_l):
        h_B_l = _bisect(xs, conds, xs.lROW, _B_funk, B_l,
                hlow, hhigh, max_iter, rel_err)
    if(B_r):
        h_B_r = _bisect(xs, conds, xs.rROW, _B_funk, B_r,
                hlow, hhigh, max_iter, rel_err)
    if(E_l):
        h_E_l = _bisect(xs, conds, xs.lROW, _E_funk, E_l,
                hlow, hhigh, max_iter, rel_err)
    if(E_r):
        h_E_r = _bisect(xs, conds, xs.rROW, _E_funk, E_r,
                hlow, hhigh, max_iter, rel_err)
    #create return variables
    h = (h_B_l, h_B_r, h_E_l, h_E_r)
    fn = _path_str_condition(xs.sheet).replace(' ', '-')
    xs = xs.copy()
    adj = fields_class.SectionBook('%s-height_adjusted' % xs.sheet)
    xs.sheet += ' (original)'
    xs.tag = 'Height Adjusted'
    adj.add_section(xs)
    sheets = ['Adjusted for Bmax left','Adjusted for Bmax right',
                'Adjusted for Emax left','Adjusted for Emax right']
    titles = ['Height Adjusted for %g mG at left ROW edge' % B_l,
                'Height Adjusted for %g mG at right ROW edge' % B_r,
                'Height Adjusted for %g kV/m at left ROW edge' % E_l,
                'Height Adjusted for %g kV/m at right ROW edge' % E_r]
    xs_tags = ['Height Adjusted']*4
    for s, ti, a, t in zip(sheets, titles, h, xs_tags):
        if(a is not None):
            #copy the input xs
            new_xs = xs.copy()
            #change the identification fields
            new_xs.sheet, new_xs.title, new_xs.tag = s, ti, t
            #adjust conductor heights
            for c_tag in tags:
                new_xs[c_tag].y += a
            #store new_xs in the SectionBook
            adj.add_section(new_xs)
    #deal with saving
    if('path' in kw):
        kw['save'] = True
    if('save' in kw):
        if(kw['save']):
            fn = _path_manage(fn + '_height_adjustments', 'xlsx', **kw)
            xl = pd.ExcelWriter(fn, engine='xlsxwriter')
            pd.DataFrame(data={'Height Addition (ft)': list(h)},
                    index=sheets).to_excel(xl, sheet_name='Adjustments')
            pd.DataFrame(dict(zip(['Height - '+s+' (ft)' for s in adj.sheets],
                    [[adj[s][j].y for j in tags] for s in adj.sheets])),
                    index=tags).to_excel(xl,
                            sheet_name='Adjusted Conductor Heights',
                            index_label='Conductor Tag')
            adj.ROW_edge_export(xl = xl)
            df, c, h = _xs_sb_diff(xs, adj)
            df.to_excel(xl, sheet_name='ROW_edge_diff', index=False,
                    columns=c, header=h)
            for xs in adj:
                xs.fields.to_excel(xl, sheet_name=xs.sheet,
                        index_label='Distance (ft)')
            xl.save()
            print('Optimal phasing results written to: %s' % fn)

    return(h, adj)

def _bisect(xs, conds, x_sample, funk, target, hlow, hhigh, max_iter, rel_err):
    #get sample x and y arrays with a single element in each
    x_sample = np.array([x_sample], dtype=float)
    y_sample = xs.sample_height*np.ones((2,), dtype=float)
    #evaluate at the bracketing values
    flow = funk(hlow, target, xs, conds, x_sample, y_sample)
    fhigh = funk(hhigh, target, xs, conds, x_sample, y_sample)
    #check that the root is bracketed
    if(flow*fhigh > 0.):
        raise(fields_class.EMFError("""
        The root is not bracketed with an upper height adjustment limit
        of %g. Rootfinding with bisection can't be performed.
            f(h_0 = %g) = %g
            f(h_1 = %g) = %g""" % (hhigh, hlow, flow, hhigh, fhigh)))
    #evaluate at a midpoint
    hmid = (hhigh + hlow)/2.0
    fmid = funk(hmid, target, xs, conds, x_sample, y_sample)
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
        fmid = funk(hmid, target, xs, conds, x_sample, y_sample)
        #increment
        count += 1
    #check if the iteration limit was hit
    if(count == max_iter):
        raise(fields_class.EMFError("""
        Failure in _bisection method. The iteration limit of %d was exseeded
        with a relative error threshold of %g. The final estimate was
        %g""" % (max_iter, rel_err, fmid)))
    return(hmid)

def _B_funk(h, target, xs, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xs.y.astype(float, copy = True)
    y[conds] += h
    #calculate B field at ROW edge
    Bx, By = fields_calcs.B_field(xs.x, y, xs.I, xs.phase,
        x_sample, y_sample)
    Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
    return(Bmax[0] - target)

def _E_funk(h, target, xs, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xs.y.astype(float, copy = True)
    y[conds] += h
    #calculate E field at ROW edge
    Ex, Ey = fields_calcs.E_field(xs.x, y, xs.subconds, xs.d_cond,
        xs.d_bund, xs.V, xs.phase, x_sample, y_sample)
    Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
    return(Emax[0] - target)

def _xs_sb_diff(xs, sb):
    """Compute the difference in ROW edge fields of all the CrossSections
    in a SectionBook object to those of a single CrossSection.
    args:
        xs - CrossSection, single xs to compare to all xss in sb
        sb - SectionBook, contains xss to compare to xs
    returns:
        df - DataFrame with columns for the names of CrossSections in sb
             and with the difference between ROW edge values, computed by
             subtracting xs values from sb values (sb.i[idx] - xs)
        c - column names
        h - refined column names (header names)"""
    #gather ROW edge differences
    L = len(sb)
    El,Er,Bl,Br = np.zeros((L,)),np.zeros((L,)),np.zeros((L,)),np.zeros((L,))
    for i in range(L):
        Bl[i] = (sb.i[i].fields.at[sb.i[i].lROW, 'Bmax']
                    - xs.fields.at[xs.lROW, 'Bmax'])
        Br[i] = (sb.i[i].fields.at[sb.i[i].rROW, 'Bmax']
                    - xs.fields.at[xs.rROW, 'Bmax'])
        El[i] = (sb.i[i].fields.at[sb.i[i].lROW, 'Emax']
                    - xs.fields.at[xs.lROW, 'Emax'])
        Er[i] = (sb.i[i].fields.at[sb.i[i].rROW, 'Emax']
                    - xs.fields.at[xs.rROW, 'Emax'])
    #create and return DataFrame
    df = pd.DataFrame(data = {
        'sheet': sb.sheets, 'Bmaxl': Bl, 'Emaxl': El, 'Bmaxr': Br, 'Emaxr': Er}
                ).sort_values('sheet')
    c = ['sheet', 'Bmaxl', 'Bmaxr', 'Emaxl', 'Emaxr']
    h = ['Cross-Section Sheet',
            'Bmax Diff - Left ROW Edge', 'Bmax Diff - Right ROW Edge',
            'Emax Diff - Left ROW Edge', 'Emax Diff - Right ROW Edge']
    return(df, c, h)
