from .. import os, np, pd, shutil, itertools, _resource_filename

from emf.emf_funks import (_path_manage, _check_extension, _is_number, _is_int,
                        _check_intable, _flatten, _sig_figs,
                        _path_str_condition, _check_to_array)

from . import fields_class
from . import fields_calcs
from . import fields_plots

def drop_template(*args, **kw):
    """Copy the emf.fields template into the current directory or a directory specified by an input string
    args:
        drop_path - string, path of copied template file"""
    #check inputs
    if(len(args) > 1):
        raise(fields_class.EMFError("""drop_template only accepts zero or one input argument. A string can be passed to specify the directory in which the template file is copied. With no arguments, the template file is copied into the current directory."""))
    elif(len(args) == 1):
        kw = {'path': args[0]}
    #get template file path
    fn_template = _resource_filename(__name__, 'fields-template.xlsx')
    #get drop path
    fn_drop = _path_manage('fields-template', 'xlsx', **kw)
    #check for existing files
    if(os.path.isfile(fn_drop)):
        raise(fields_class.EMFError('A file with the path "%s" already exists. Move/delete it or pass a new path string to drop_template().' % fn_drop))
    #copy and notify
    shutil.copyfile(fn_template, fn_drop)
    print('emf.fields template written to: %s' % fn_drop)

def run(template_path, **kw):
    """Import the templates in an excel file with the path 'template_path' then generate a workbook of all fields results and lots of plots. Use the 'path' keyword argument to specify a destination for the output, otherwise it will be saved to the template's directory. Returns a SectionBook object.
    args:
        template_path - path to cross section template excel workbook
    kw:
        sheets - a list of sheet names to load, default is all sheets
        path - string, destination/filename for saved files
        format - string, saved plot format (usually 'png' or 'pdf')
        xmax - cutoff distance from ROW center in plots
    returns:
        sb - SectionBook object created from the template file"""
    #force saving for the plotting functions if there is no 'path' keyword
    if(not ('path' in kw)):
        kw['save'] = True
        #also direct output files to the same directory as the template
        template_dir = os.path.dirname(template_path)
        if(template_dir):
            kw['path'] = template_dir
    #import templates
    if('sheets' in kw):
        sheets = kw['sheets']
    else:
        sheets = 'all'
    sb = load_template(template_path, sheets)
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

def load_template(file_path, sheets='all'):
    """Import conductor data from an excel template, loading each conductor in the template into a Conductor object, each Conductor object into a CrossSection object, and each CrossSection object into a SectionBook object. The SectionBook object is returned.
    args:
        template_path - string, path to cross section template excel
                        workbook
    optional args:
        sheets - list of strings, a list of sheet names to load, default is
                'all' sheets
    returns:
        sb - SectionBook object containing each sheet of the excel template
             as a CrossSection object"""

    #import the cross sections as a dictionary of pandas DataFrames, also
    #getting a list of the ordered sheets
    file_path = _check_extension(file_path, 'xlsx', """Templates must be excel workbooks. The input target path "%s" is not recognized as an excel file""" % file_path)
    xl = pd.ExcelFile(file_path)
    if(sheets == 'all'):
        sheets = xl.sheet_names
    else:
        for sh in sheets:
            if(sh not in xl.sheet_names):
                raise(fields_class.EMFError('Input sheet "%s" was not found in the target workbook' % sh))
    frames = xl.parse(sheetname=None, skiprows=4, parse_cols=16, header=None)
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
        xs.group = misc[0]
        xs.title = str(misc[1])
        #check for duplicate title inputs
        if(xs.title in titles):
            raise(fields_class.EMFError("""Cross-sections should have unique title entries. Title "%s" in sheet "%s" is used by at least one other sheet.""" % (xs.title, k)))
        else:
            titles.append(xs.title)
        xs.soil_resistivity = misc[3]
        xs.max_dist = misc[4]
        xs.step = misc[5]
        xs.sample_height = misc[6]
        xs.lROW = misc[7]
        xs.rROW = misc[8]
        #load hot conductors
        for i in range(df[3].dropna().shape[0]):
            #initialize a Conductor
            cond = fields_class.Conductor(df[2].iat[i])
            #cond.freq = misc[2]
            cond.x = df[3].iat[i]
            cond.y = df[4].iat[i]
            cond.subconds = df.iat[i,5]
            cond.d_cond = df.iat[i,6]
            cond.d_bund = df.iat[i,7]
            cond.V = df.iat[i,8]
            cond.I = df.iat[i,9]
            cond.phase = df.iat[i,10]
            xs.add_conductor(cond)
        #load grounded conductors
        for i in range(df[11].dropna().shape[0]):
            #initialize a Conductor
            cond = fields_class.Conductor(df.iat[i,11])
            #check for conductors with identical names (names/labels)
            #cond.freq = misc[2]
            cond.x = df.iat[i,12]
            cond.y = df.iat[i,13]
            cond.subconds = 1.
            cond.d_cond = df.iat[i,14]
            cond.d_bund = df.iat[i,14]
            cond.V = 0.
            cond.I = df.iat[i,15]
            cond.phase = df.iat[i,16]
            xs.add_conductor(cond)
        #add the CrossSection object to the SectionBook
        #   Fields are automatically computed or updated once the 'fields'
        #   property of a CrossSection is accessed (directly or indirectly).
        sb.add_section(xs)
    #return the SectionBook object
    return(sb)

def optimize_phasing(xs, circuits='all', **kw):
    """Permute the phasing of non-grounded conductors and find the arrangement that results in the lowest fields at the left and right edge of the ROW. The number of hot conductors must be a multiple of three. The phases of consecutive groups of three conductors are swapped around, assuming that those groups represent a single three-phase transfer circuit.
    args:
        xs - target CrossSection object
    optional args:
        circuits - list of lists of Conductor names, or 'all'. If a list of
                    lists, each sublist contains the Conductor names of the
                    Conductors comprising a single circuit. If 'all',
                    circuits are assumed to be consecutive groups of three
                    conductors. (consecutive according to the order in which
                    hot, or non-grounded, conductors were added to the
                    CrossSection)
    kw:
        save - bool, toggle saving of the results DataFrame to an excel book
        path - string, location/filename for saved results workbook,
                forces saving even if no 'save' keyword is used.
    returns:
        res - pandas DataFrame listing conductor phasings that optimize
                electric and magnetic fields at both ROW edges.
        opt - new SectionBook object containing the permuted phasings that
                optimize the E and B fields at the left and right ROW edges."""

    if(circuits == 'all'):
        #number of hot wires
        hot = xs.hot
        N = len(xs.hot)
        #check the number of hot lines
        if(N % 3 != 0):
            raise(fields_class.EMFError("""The number of hot (not grounded) conductors must be a multiple of three for phase optimization with 'all' circuits. Circuits are assumed to be three-phase and conductors comprising each circuit are assumed to be consecutive groups of three, in the order that they appear in the template. The number of hot conductors is not a multiple of three in the CrossSection named: %s""" % xs.sheet))
        #number of circuits, groups of 3 hot conductors
        G = int(N/3)
        #circuits, consecutive groups of three conductors
        circuits = [[] for i in range(G)]
        for i in range(G):
            for j in range(3):
                circuits[i].append(hot[i*3 + j].name)
    else:
        #check that all conductor names are present and refer to hot conds
        gnd = xs.gnd
        for circ in circuits:
            for name in circ:
                if(xs[name] is None):
                    raise(fields_class.EMFError("""Unrecognized conductor name: %s All conductor names must refer to Conductor objects in the target CrossSecton object.""" % repr(name)))
                if(xs[name] in gnd):
                    raise(fields_class.EMFError("""Only phasing of non-grounded Conductors can be permuted. Name "%s" refers to a grounded Conductor""" % repr(name)))
    #convert the conductor names to integer indices in xs.conds
    for i in range(len(circuits)):
        for j in range(len(circuits[i])):
            circuits[i][j] = xs._name2idx[circuits[i][j]]
    #all permutations of the phases of each circuit
    perm = []
    for c in circuits:
        perm.append(list(itertools.permutations(c)))
    #all possible arrangements of line phasings, 6 permutations for each circuit
    #so 6^(N/3) total line arrangements. Leave P as a generator to avoid storing
    #a huge, factorial sized array of indices
    P = itertools.product(*perm)
    #variables to find the minima with respect to each field and ROW edge
    B_left_min, B_right_min = np.inf, np.inf
    E_left_min, E_right_min = np.inf, np.inf
    #get coordinates of the ROW edges
    x_ROW = np.array([xs.lROW, xs.rROW], dtype=float)
    y_ROW = xs.sample_height*np.ones((2,), dtype=float)
    #pull conductor data
    x, y, I, V, phase = xs.x, xs.y, xs.I, xs.V, xs.phase
    subconds, d_cond, d_bund = xs.subconds, xs.d_cond, xs.d_bund
    #phasing array to swap the elements around in
    phase_swap = phase.copy()
    #store a flattened version of the conductor indices for swapping
    conds = np.array([i for j in circuits for i in j], dtype=int)
    #loop through all possible arrangements in P
    for arr in P:
        #flatten the new arrangement
        new_arr = _flatten(arr)
        #swap phases according to the new phasing arrangement
        phase_swap[conds] = phase[new_arr]
        #calculate fields with index swapped phases
        Ex, Ey = fields_calcs.E_field(x, y, subconds, d_cond, d_bund, V,
                phase_swap, x_ROW, y_ROW)
        Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
        Bx, By = fields_calcs.B_field(x, y, I, phase_swap, x_ROW, y_ROW)
        Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
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
        'Optimal Phasing - Bmax Left ROW Edge': phase[B_left_arr],
        'Optimal Phasing - Bmax Right ROW Edge': phase[B_right_arr],
        'Optimal Phasing - Emax Left ROW Edge': phase[E_left_arr],
        'Optimal Phasing - Emax Right ROW Edge': phase[E_right_arr]},
        index=[xs.conds[i].name for i in conds])
    #compile a new sectionbook with the optimal phasings
    fn = _path_str_condition(xs.sheet).replace(' ', '-')
    xs = xs.copy()
    opt = fields_class.SectionBook(xs.sheet + '-optimal_phasing')
    xs.sheet += ' (original)'
    xs.group = 'Phase Optimized'
    opt.add_section(xs)
    names = ['Optimized for Bmax left','Optimized for Bmax right',
            'Optimized for Emax left','Optimized for Emax right']
    for n, ti in zip(names, results.columns):
        #copy the input xs
        new_xs = xs.copy()
        #change the identification fields
        new_xs.sheet, new_xs.title, new_xs.group = n, ti, 'Phase Optimized'
        #swap the conductor phasings
        for c in new_xs.hot:
            t = c.name
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
            xl = pd.ExcelWriter(fn, engine='xlsxwriter')
            results.to_excel(xl, index_label='Conductor Name',
                sheet_name='phase_assignments')
            opt.ROW_edge_export(xl=xl)
            abs_dif, rel_dif = sb_xs_compare(opt, opt.sheets[0])
            abs_dif.to_excel(xl, sheet_name='Abs Difference at ROW Edges')
            rel_dif.to_excel(xl, sheet_name='Rel Difference at ROW Edges')
            for xs in opt:
                xs.fields.to_excel(xl, sheet_name=xs.sheet,
                        index_label='Distance (ft)')
            xl.save()
            print('Phase optimization results written to: %s' % fn)

    return(results, opt)

def target_fields(xs, names='all', B_l=False, B_r=False, E_l=False, E_r=False,
        max_iter=1e3, rel_err=1.0e-6, hhigh=1.0e6, **kw):
    """Increase conductor y coordinates until max fields at ROW edges are below thresholds. All selected conductors are adjusted by the same amount. If any of the thresholds are empty or false, None is returned for their adjustment result. If
    args:
        xs - CrossSection object to perform adjustments on
    optional args:
        names - iterable of Conductor names, identifying which ones to raise,
               or 'all' to include all Conductors in the CrossSection
        B_l - magnetic field threshold at left ROW edge*
        B_r - magnetic field threshold at right ROW edge*
        E_l - electric field threshold at left ROW edge*
        E_r - electric field threshold at right ROW edge*

            *an implicitly False input will ignore that field-edge
            combination, return None in the return variable 'h', and cause
            the returned SectionBook to omit that field-edge combo.

        max_iter - maximum number of _bisection iterations allowed,
                    default is 1e3
        rel_err - tolerance threshold for relative error (e.g. 0.01 is 1 %),
                    default is 1e-6.
        hhigh - upper limit of the height adjustment, default is 1.0e6
    kw:
        save - toggle saving of the results DataFrame to an excel book
        path - location/filename for saved results workbook, forces saving
                even if no 'save' keyword is used.
    returns:
        h - height adjustments necessary for E and B fields at left and
            right ROW edges. The ordering is: (B_left, B_right, E_left, E_right)
        adj - a new SectionBook object with the adjusted conductor heights
                for each scenario in a CrossSection"""
    #convert 'all' inputs to numeric indices
    if(names == 'all'):
        names = xs.names
    #flattened indices
    temp = xs.names
    conds = np.array([temp.index(i) for i in names])
    #maximum number of iterations and relative error tolerance
    hlow = 0.0
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
    xs.group = 'Height Adjusted'
    adj.add_section(xs)
    sheets = ['Adjusted for Bmax left','Adjusted for Bmax right',
                'Adjusted for Emax left','Adjusted for Emax right']
    titles = ['Height Adjusted for %g mG at left ROW edge' % B_l,
                'Height Adjusted for %g mG at right ROW edge' % B_r,
                'Height Adjusted for %g kV/m at left ROW edge' % E_l,
                'Height Adjusted for %g kV/m at right ROW edge' % E_r]
    xs_groups = ['Height Adjusted']*4
    for s, ti, a, t in zip(sheets, titles, h, xs_groups):
        if(a is not None):
            #copy the input xs
            new_xs = xs.copy()
            #change the identification fields
            new_xs.sheet, new_xs.title, new_xs.group = s, ti, t
            #adjust conductor heights
            for c_name in names:
                new_xs[c_name].y += a
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
                    [[adj[s][j].y for j in names] for s in adj.sheets])),
                    index=names).to_excel(xl,
                            sheet_name='Adjusted Conductor Heights',
                            index_label='Conductor Names')
            adj.ROW_edge_export(xl=xl)
            abs_dif, rel_dif = sb_xs_compare(adj, adj.sheets[0])
            abs_dif.to_excel(xl, sheet_name='Abs Difference at ROW Edges')
            rel_dif.to_excel(xl, sheet_name='Rel Difference at ROW Edges')
            for xs in adj:
                xs.fields.to_excel(xl, sheet_name=xs.sheet,
                        index_label='Distance (ft)')
            xl.save()
            print('Height adjustment results written to: %s' % fn)

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
    Failure in bisection method (_bisect). The iteration limit of %d
    was exceeded with a relative error threshold of %g. The absolute error
    of the final estimate was %g. Try reducing the relative error
    threshold (rel_err) or increasing the maximum iteration limit
    (max_iter).""" % (max_iter, rel_err, fmid)))
    return(hmid)

def _B_funk(h, target, xs, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xs.y.copy()
    y[conds] += h
    #calculate B field at ROW edge
    Bx, By = fields_calcs.B_field(xs.x, y, xs.I, xs.phase,
        x_sample, y_sample)
    Bx, By, Bprod, Bmax = fields_calcs.phasors_to_magnitudes(Bx, By)
    return(Bmax[0] - target)

def _E_funk(h, target, xs, conds, x_sample, y_sample):
    #adjust conductor heights
    y = xs.y.copy()
    y[conds] += h
    #calculate E field at ROW edge
    Ex, Ey = fields_calcs.E_field(xs.x, y, xs.subconds, xs.d_cond,
        xs.d_bund, xs.V, xs.phase, x_sample, y_sample)
    Ex, Ey, Eprod, Emax = fields_calcs.phasors_to_magnitudes(Ex, Ey)
    return(Emax[0] - target)

def sb_xs_compare(sb, sheet):
    """Compute the absolute and percentage difference between the maximum fields of one CrossSection in a SectionBook and all the other CrossSections
    args:
        sb - SectionBook object
        sheet - sheet string of the CrossSection in sb that should be
                compared with all the other CrossSections
    returns:
        abs_dif - DataFrame with the absolute difference of the maximum
                  fields at ROW edges
        rel_dif - DataFrame with the relative difference of the maximum
                  fields at ROW edges (multiply by 100 to get percentage dif)"""

    rem = sb.ROW_edge_max
    abs_dif = rem - rem.loc[sheet]
    rel_dif = (rem - rem.loc[sheet])/rem
    return(abs_dif, rel_dif)
