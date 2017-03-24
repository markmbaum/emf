"""This module contains functions for reading the files generated/used by the ENVIRO program, which is part of the same software suite as the SUBCALC program. ENVIRO is very similary to FIELDS. It models the electric and magnetic fields in the vicinity of groups of parallel power lines. There are, however, some differences between ENVIRO and FIELDS that might impact modeling results:
 * ENVIRO asks for more physical parameters of the power lines (resistance and reactance)
 * ENVIRO can take rain, fog, and snow into account somehow
 * ENVIRO doesn't appear to be able to model underground lines
Aside from those differences, ENVIRO appears to perform the same calculations as FIELDS and emf.fields

ENVIRO makes use of several input/output files with different extensions:
 * i01, unlabeled input data
 * o01, formatted model results
 * o02, formatted input data
 * o03, table of audible noise results without column headers
 * o04, table of electric field results without column headers
 * o05, table of magnetic field results without column headers

The functions in this module read some of the above file types.
"""

from .. import os, np, pd, plt

import fields_funks
import fields_class
import fields_plots

def _advance_lines(ifile, n):
    for i in range(n):
        line = ifile.readline()
    return(line)

def read_o01(fn):
    """Read a formatted output file of the ENVIRO program and extract the electric and magnetic fields across the model's transect, returning them in a DataFrame
    args:
        fn - str, path to the target file, must have a ".o01" extension
    returns:
        df - DataFrame of fields across the transect"""

    #check extension
    if((fn[-4:] != '.o01') and (fn[-4:] != '.O01')):
        raise(EMFError('Target file must have an .o01 or .O01 extension'))
    #open file for reading
    ifile = open(fn, 'r')
    #read to the electric field profile
    line = ifile.readline()
    while('ELECTRIC FIELD PROFILE' not in line):
        line = ifile.readline()
    #pick out the sample height
    line = ifile.readline()
    sample_height = float(line.split()[2])
    line = _advance_lines(ifile, 8)
    #capture the electric field profile
    dist, Emax, Ey, Ex, Eprod = [], [], [], [], []
    while(line):
        s = [fields_funks._check_intable(i) for i in line.split()]
        dist.append(s[0])
        Emax.append(s[2])
        Ey.append(s[4])
        Ex.append(s[5])
        Eprod.append(np.nan)
        line = ifile.readline().strip()
    #read to the magnetic field profile
    line = ifile.readline()
    while('MAGNETIC FIELD PROFILE' not in line):
        line = ifile.readline()
    line = _advance_lines(ifile, 10)
    #capture the electric field profile
    Bmax, By, Bx, Bprod = [], [], [], []
    while(line):
        s = [fields_funks._check_intable(i) for i in line.split()]
        Bmax.append(s[2])
        By.append(s[4])
        Bx.append(s[5])
        Bprod.append(s[6])
        line = ifile.readline().strip()

    ifile.close()

    df = pd.DataFrame(dict(
            Emax=Emax, Ex=Ex, Ey=Ey, Eprod=Eprod,
            Bmax=Bmax, Bx=Bx, By=By, Bprod=Bprod),
            index=dist)
    df.index.name = 'Distance (ft)'

    return(df)

def read_o01s(dir_name):
    """Extract fields from all "o01" files encountered in a directory, returning them in a dictionary keyed by the filenames (without extensions)
    args:
        dir_name - str, path to directory to search for o01 files
    returns:
        dfs - dict of DataFrames created with read_o01, a single DataFrame
              for each o01 file encountered"""

    dfs = dict()
    for fn in os.listdir(dir_name):
        if(len(fn) > 4):
            if(fn[-4:].lower() == '.o01'):
                k = fn.replace('.o01', '').replace('.O01', '')
                dfs[k] = read_o01(os.path.join(dir_name, fn))

    return(dfs)

def read_o02(fn):
    """Read a formatted input file of the ENVIRO program, which occur with ".o02" extensions. The relevant inputs are loaded into a CrossSection object. Because the inputs may specify a group of samples with non-uniform spacing, a DataFrame containing the fields at the points specified in the model (generated from the returned CrossSection) is also returned.
    args:
        fn - str, path to the formatted input file, must have an o02 extension
    returns:
        xs - CrossSection object constructed from parameters in the input file
        df - DataFrame containing fields at the sample points given in the
             input file"""

    #check extension
    #check extension
    if((fn[-4:] != '.o02') and (fn[-4:] != '.O02')):
        raise(EMFError('Target file must have an .o02 or .O02 extension'))
    #open file for reading
    ifile = open(fn, 'r')
    line = _advance_lines(ifile, 2)

    #read general parameters
    gen = dict()
    while('====' not in line):
        try:
            idx = line.index(':')
        except(ValueError):
            pass
        else:
            gen[line[:idx].replace('.','').lower()] = line[idx+1:].strip()
        line = ifile.readline()

    #read bundle specifications
    line = _advance_lines(ifile, 3)
    C = dict()
    params = ['V', 'I', 'phase', 'subconds', 'x', 'y']
    while('====' not in line):
        s = line.split()
        if(s):
            #for each bundle creat a subdictionary containing modeling params
            #   and the bundle name
            k = s[0]
            C[k] = dict()
            C[k]['name'] = s[0] + '-' + s[1] + s[2]
            p = [fields_funks._check_intable(i) for i in (s[3:7] + s[8:])]
            C[k]['params'] = dict(zip(params, p))
        line = ifile.readline()

    #read regular bundle specifications
    line = _advance_lines(ifile, 4)
    while('====' not in line):
        s = line.split()
        if(s):
            k = s[0]
            C[k]['params']['d_cond'] = float(s[3])
            C[k]['params']['d_bund'] = float(s[2]) + float(s[3])
        line = ifile.readline()

    #read irregular bundle specifications
    line = ifile.readline()
    while('====' not in line):
        #read to the next bundle
        while(('====' not in line) and ('No.' not in line)):
            line = ifile.readline()
        #check if a bundle is There
        if('No.' in line):
            #advance one line to the table of irregular bundle specs
            line = ifile.readline()
            #read bundle specs and number into a dict
            irr = {'x': [], 'y': [], 'd_cond': [], 'count': 0}
            while(line.strip() != ''):
                irr['count'] += 1
                s = line.split()
                k = s[0]
                irr['x'].append(float(s[2]))
                irr['y'].append(float(s[3]))
                irr['d_cond'].append(float(s[4]))
                line = ifile.readline()
            #replace the irregular bundle with the new conductors
            p = C[k]['params']
            n = C[k]['name']
            if(p['subconds'] != irr['count']):
                raise(fields_class.EMFError('Number of subconductors indicated in the "BUNDLE SPECIFICATIONS" section for bundle "%s" does not match the number of subconductors printed in the "IRREGULAR BUNDLES" section. Cannot parse.' % k))
            del(C[k])
            for i in range(irr['count']):
                newk = k + '-' + str(i+1)
                C[newk] = dict()
                C[newk]['name'] = n + 'sub' + str(i+1)
                C[newk]['params'] = dict(p)
                C[newk]['params']['x'] += irr['x'][i]/12.0
                C[newk]['params']['y'] += irr['y'][i]/12.0
                C[newk]['params']['d_cond'] = irr['d_cond'][i]
                C[newk]['params']['subconds'] = 1
                C[newk]['params']['d_bund'] = irr['d_cond'][i]
                C[newk]['params']['I'] /= irr['count']

    #read lateral profile info
    line = _advance_lines(ifile, 3)
    lat = []
    while('====' not in line):
        s = line.split()
        if(s):
            start = float(s[1])
            stop = float(s[2])
            step = float(s[3])
            lat.append(list(np.linspace(start, stop, 1 + (stop - start)/step)))
        line = ifile.readline()
    x_sample = np.unique(np.array(fields_funks._flatten(lat)))

    #ignore the weather section, close the file
    ifile.close()

    #create a CrossSection object
    xs = fields_class.CrossSection(os.path.basename(fn).replace('.o02',''),
        [fields_class.Conductor(C[k]['name'], C[k]['params']) for k in C])
    xs.title = gen['title']
    xs.sample_height = gen['field sensor height (feet)']
    xs.max_dist = np.max(np.abs(x_sample))
    if(len(lat) == 3):
        xs.lROW = lat[0][-1]
        xs.rROW = lat[-1][0]
    else:
        xs.lROW = -xs.max_dist
        xs.rROW = xs.max_dist

    return(xs, xs.sample(x_sample))

def compare_o01_o02(fn_o01, fn_o02, **kw):
    """Compare the results of an ENVIRO model to those of emf.fields. Inputs to the ENVIRO model are read from an o02 file into an emf.fields CrossSection object to compute fields. The computed fields are compared to the results of the same ENVIRO model, contained in a o01 file.
    args:
        fn_o01 - str, the path to a o01 file created by ENVIRO, which
                 contains formatted model results
        fn_o02 - str, the path to a o02 file created by ENVIRO, which
                 contains formatted input parameters of the model
    returns:
        pan - pandas Panel object with the results of both programs and
              DataFrames with the absolute and percentage error between them
        figs - list of matplotlib Figure objects, one comparing magnetic
               field results and another comparing electric field results"""

    #create output filename for later
    fn_out = os.path.basename(fn_o01).replace('.o01','').replace('.O01','')

    #read in the enviro output
    df1 = read_o01(fn_o01)
    #read enviro input parameters and calculate fields with emf.fields
    #  for comparison
    xs, df2 = read_o02(fn_o02)

    #get the horizontal distance samples
    x = df1.index.values
    #check that they match the samples in the other df
    if(not np.array_equal(x, df2.index.values)):
        raise(fields_class.EMFError('The x samples found in the target o01 file do not match those found in the target o02 file. Make sure that the target files are were generated by the same ENVIRO model.'))

    #create a panel (the frame names matter in the _plot_comparison function
    #called below)
    frames = ['ENVIRO-results', 'emf.fields-results',
            'absolute-difference', 'percent-difference']
    pan = pd.Panel(data={
            frames[0] : df1,
            frames[1]: df2,
            frames[2]: df2 - df1,
            frames[3]: 100*(df2 - df1)/df1})

    #make plots of the absolute and percent error
    figs = fields_plots._plot_comparison(xs, pan, 'ENVIRO', **kw)

    if('path' in kw):
        kw['save'] = True
    if('save' in kw):
        if(kw['save']):
            fn = fields_funks._path_manage(fn_out, '.xlsx', **kw)
            xl = pd.ExcelWriter(fn, engine='xlsxwriter')
            for f in frames:
                pan[f].to_excel(xl, sheet_name=f, index_label='Distance (ft)')
            xl.save()
            print('ENVIRO comparison book saved to: %s' % fn)

    return(pan, figs)
