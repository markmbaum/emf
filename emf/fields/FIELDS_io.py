"""The FIELDS_io module provides a handful of functions for streamlining the file management process required by the old modeling program called FIELDS. There are functions for reading the DAT files generated by FIELDS and converting them into standard csv files (targeting a single DAT file (convert_DAT) or all DAT files in a given directory and all its subdirectories (convert_DAT_crawl, which has an additional option to bundle all DAT files found in a given directory inside single xlsx files for easy access.)). There are also functions that write the data in existing CrossSection objects to formatted text files with FLD extensions, to be used as input files with FIELDS instead of using the program's troublesome and error-prone menus to input data. These functions are to_FLD, to_FLDs, and to_FLDs_crawl. CrossSection objects are generated by reading from a template excel file, using the load_template function. Converting all the sheets in a template excel file to FLD files is a one line process: fields.to_FLDs(template_file.xlsx)"""

from .. import os, pd, glob

import fields_funks
import fields_class

#-------------------------------------------------------------------------------
#FUNCTIONS FOR GENERATING INPUT .FLD FILES

def _is_int(x):
    """check if a number is an integer, will break on non-numeric entries"""
    if(int(x) == x):
        return(True)
    else:
        return(False)

#function formats entries and writes them into a .FLD file targeted by ofile
def _write_FLD_entries(ofile, *entries):
    """format entries and writes them into a .FLD file targeted by ofile"""
    for entry in entries:
        if(fields_funks._is_number(entry)):
            if(_is_int(entry)):
                w = '{:< d} \n'.format(int(entry))
            else:
                w = '{:< .2f} \n'.format(float(entry))
            if('.' in w):
                idx = w.index('.')
                if('0' in w[:idx]):
                    idx = w.index('0')
                    if(not fields_funks._is_number(w[idx-1])):
                        w = w[:idx] + w[idx+1:]
        else:
            w = str(entry) + '\n'
        ofile.write(w)

def to_FLD(xs, **kw):
    """Create an FLD input file for FIELDS from a CrossSection object
    args:
        xs - CrossSection object
    kw:
        path - output file destination"""
    #check input
    if(not isinstance(xs, fields_class.CrossSection)):
        raise(fields_class.EMFError("""Input argument to to_FLD() must be a CrossSection object, not an input of type: %s Use to_FLDs() for a SectionBook object.""" % str(type(xs))))
    #get a filename
    fn = fields_funks._path_manage(xs.sheet, 'FLD', **kw)
    #write the .FLD file
    ofile = open(fn, 'w')
    #miscellaneous stuff first
    _write_FLD_entries(ofile, xs.sheet, xs.title, 60, xs.soil_resistivity,
            xs.max_dist, xs.step, xs.sample_height, xs.lROW, xs.rROW)
    #number of conductors and ground wires
    Lconds = len(xs.hot)
    _write_FLD_entries(ofile, Lconds)
    Lgrounds = len(xs.gnd)
    _write_FLD_entries(ofile, Lgrounds)
    #write the hot and gnd conductor data in the same format
    for c in xs.hot + xs.gnd:
        _write_FLD_entries(ofile, c.name, c.x, c.y, c.subconds, c.d_cond,
                c.d_bund, 'ED!(I)', c.I, c.V, c.phase)
    #write the ground wire data a second time, in a different format
    for c in xs.gnd:
        _write_FLD_entries(ofile, c.name, c.x, c.y, c.d_cond, 0, 0)
    #close/save
    ofile.close()
    print('FLD file generated: %s' % fn)

def to_FLDs(*args, **kw):
    """Load or recieve a template workbook of CrossSections and convert them all to FLD files
    args:
        can either be a path string to a target template workbook or an
        existing SectionBook object
    kw:
        path - output destination for FLD files"""
    if(type(args[0]) == str):
        #load the template
        sb = fields_funks.load_template(args[0])
        if('path' not in kw):
            kw['path'] = os.path.dirname(args[0])
    elif(isinstance(args[0], fields_class.SectionBook)):
        sb = args[0]
    else:
        raise(fields_class.EMFError("""Input argument to to_FLDs() must be a filepath or a SectionBook."""))
    #check for duplicate sheets
    sheets = []
    for xs in sb:
        if(xs.sheet in sheets):
            raise(fields_class.EMFError("""Can't create FLD files because of duplicate CrossSection names. Name "%s" is used at least twice.""" % xs.sheet))
        else:
            sheets.append(xs.sheet)
    #generate FLD files
    for xs in sb:
        to_FLD(xs, **kw)

def to_FLDs_crawl(dir_name, **kw):
    """crawl a directory and all of its subdirectories for excel workbooks that can be passed to create_FLDs(). The FLD files are generated in the same directory as the template book they come from."""

    #get input directory's file and subdir names
    dir_contents = glob.glob(dir_name)

    #loop over the dir_contents
    for dir_element in dir_contents:
        if(dir_element[-5:] == '.xlsx'):
            #operate on the file
            try:
                to_FLDs(dir_element, path = os.path.dirname(dir_element))
            except(KeyError, ValueError, IOError, fields_class.EMFError) as e:
                print('failure to write FLD files from:\n\t%s' % dir_element)
                if(type(e) is fields_class.EMFError):
                    print('\n\tBecause of EMFError:' + str(e))
        else:
            #if the element is a subdirecory, crawl it
            if(os.path.isdir(dir_element)):
                to_FLDs_crawl(os.path.join(dir_element, '*'))

#------------------------------------------------------------------------------
#FUNCTIONS FOR READING/CONVERTING .FLD FILES

def read_FLD(file_path):
    """Read a FLD file in to a CrossSection object
    args:
        file_path - string, path to FLD file
    returns:
        xs - CrossSection object representing the information in the FLD file"""
    #check extension
    fields_funks._check_extension(file_path, 'FLD', """Input file must have a '.FLD' extention""")
    #read the FLD file into a list of lines without whitespace on the right
    with open(file_path, 'r') as ifile:
        fld = [i.rstrip() for i in ifile.readlines()]
    #initialize a CrossSection object
    xs = fields_class.CrossSection(fld[0])
    xs.title = fld[1]
    xs.max_dist = float(fld[4])
    xs.step = float(fld[5])
    xs.sample_height = float(fld[6])
    xs.lROW = float(fld[7])
    xs.rROW = float(fld[8])
    Lconds = int(fld[9])
    Lgrounds = int(fld[10])
    #create hot Conductors for the CrossSection
    n = 11
    i = 0
    while(i < Lconds):
        xs.add_conductor(fields_class.Conductor(fld[n],
                dict(x=fld[n+1], y=fld[n+2], subconds=fld[n+3],
                d_cond=fld[n+4], d_bund=fld[n+5], I=fld[n+7], V=fld[n+8],
                phase=fld[n+9])
            )
        )
        n += 10
        i += 1
    #create ground Conductors
    i = 0
    while(i < Lgrounds):
        #sometimes there's a blank conductor?
        if(fld[n]):
            xs.add_conductor(fields_class.Conductor(fld[n],
                    dict(x=fld[n+1], y=fld[n+2], subconds=1,
                    d_cond=fld[n+3], d_bund=fld[n+3], I=0, V=0, phase=0)
                )
            )
            n += 6
            i += 1
        else:
            n += 1
            go = True
            while(go):
                if(fld[n]):
                    if(fld[n][0] != ' ') and ('ED!(I)' not in fld[n]):
                        go = False
                        n -= 1
                n += 1

    return(xs)

#------------------------------------------------------------------------------
#FUNCTIONS FOR CONVERTING OUTPUT .DAT FILES TO CSV/excel FILES

def read_DAT(file_path):
    """Read a DAT file, which can have some funky extra characters if the numbers are too large (percent signs)
    args:
        file_path - string, path to DAT file
    returns:
        pandas DataFrame with 'Bx','By','Bprod','Bmax','Ex','Ey','Eprod',
        and 'Emax' columns, and the distance ('x') in the index"""

    #check that the target file is a DAT
    fields_funks._check_extension(file_path, 'DAT', """Input file must have a '.DAT' extension.""")
    #load data
    und_message = 'Electric Field cannot be computed for underground circuit'
    und_only = False
    with open(file_path,'r') as ifile:
        #read through the header
        line = ifile.readline()
        while(not ('----' in line)):
            line = ifile.readline()
            if(und_message in line):
                und_only = True
        #get the data
        data = [[] for i in range(9)]
        line = ifile.readline()
        while(line):
            if(line):
                #check strange % signs printed with large numbers
                if(line[0][0] == '%'):
                    line += ifile.readline()
            l = [i.replace('%','') for i in line.split()]
            if(l and all([fields_funks._is_number(i) for i in l])):
                for i in range(5):
                    data[i].append(float(l[i]))
                if(not (und_only)):
                    for i in range(5,9):
                        data[i].append(float(l[i]))
                else:
                    for i in range(5,9):
                        data[i].append(0.)
            line = ifile.readline()
    return(pd.DataFrame(data=
        dict(zip(['Bx','By','Bprod','Bmax','Ex','Ey','Eprod','Emax'],data[1:])),
        index = data[0]))

def convert_DAT(file_path, **kw):
    """read a DAT file and write it to a csv
    args:
        file_path - target DAT file
    kw:
        path - output location/name for csv file"""
    #get the DAT data
    df = read_DAT(file_path)
    #write it to a csv
    fn = fields_funks._path_manage(os.path.basename(file_path)[:file_path.index('.')],
        'csv', **kw)
    with open(fn, 'w') as ofile:
        df.to_csv(ofile, index_label = 'Distance (ft)')
        print('DAT converted to csv: "%s"' % fn)

def convert_DAT_crawl(dir_name, **kw):
    """crawl a directory and all of its subdirectories for .DAT files that can be passed to DAT_to_csv() for output re-formatting and optional plotting.
    args:
        dir_name - Directory to initiate the crawl in.
                    To designate the current directory, use '*'
    kw:
        bundle - bool, if True, all DAT files found in the same directory
                are written to a common excel workbook. If False or absent,
                the DAT files are simply written to individual csv files."""

    #get input directory's file and subdir names
    dir_contents = glob.glob(dir_name)

    #get the bundling option
    if('bundle' in kw):
        bundle = kw['bundle']
    else:
        bundle = False
    #if bundle is True, set a flag for the existence of an ExcelWriter object
    if(bundle):
        xl_flag = False

    #loop over the dir_contents
    for dir_element in dir_contents:
        #if the file is a .DAT, convert it
        if(dir_element[-4:] == '.DAT'):
            #if bundling is on, read the DAT and write it to the excel book
            if(bundle):
                #if an ExcelWriter doesn't exist, initialize one and set flag
                if(not xl_flag):
                    xl_flag = True
                    fn = os.path.join(os.path.dirname(dir_element),
                            'converted_DATs.xlsx')
                    xl = pd.ExcelWriter(fn, engine = 'xlsxwriter')
                #use the DAT's filename as the sheet name, without extension
                sn = os.path.basename(dir_element).replace('.DAT','')
                #read the DAT into a DataFrame and send it to the ExcelWriter
                read_DAT(dir_element).to_excel(xl, sheet_name = sn,
                    index_label = 'Distance (ft)')
            else:
                #without bundling, write an individual csv file
                convert_DAT(dir_element, path = os.path.dirname(dir_name))
        else:
            #if there's a period in the dir_element, it's not a directory
            if(os.path.isdir(dir_element)):
                convert_DAT_crawl(os.path.join(dir_element, '*'), **kw)
    #close/save the excel bundles DAT results
    if(bundle and xl_flag):
        xl.save()
        print('Converted DAT files written to "%s"' % fn)
