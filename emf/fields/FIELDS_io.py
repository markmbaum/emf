from .. import os
from .. import glob
from .. import pd

import fields_funks
import fields_class

#-------------------------------------------------------------------------------
#FUNCTIONS FOR GENERATING INPUT .FLD FILES

def _is_int(x):
    """check if a number is an integer, will break on non-numeric entries"""
    if(int(x) == x):
        return True
    else:
        return False

#function formats entries and writes them into a .FLD file targeted by ofile
def _write_entries(ofile, *entries):
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

def to_FLD(xc, **kwargs):
    """Create an FLD input file for FIELDS out of a CrossSection object
    args:
        xc - CrossSection object
    kwargs:
        path - output file destination"""
    #check input
    if(not isinstance(xc, fields_class.CrossSection)):
        raise(fields_class.EMFError("""
        Input argument to to_FLD() must be a CrossSection object, not an
        input of type: %s
        Use to_FLDs() for a SectionBook object.""" % str(type(xc))))
    #get a filename
    fn = fields_funks._path_manage(xc.sheet, 'FLD', **kwargs)
    #write the .FLD file
    ofile = open(fn, 'w')
    #miscellaneous stuff first
    _write_entries(ofile, xc.sheet, xc.title, 60, xc.soil_resistivity,
            xc.max_dist, xc.step, xc.sample_height, xc.lROW, xc.rROW)
    #number of conductors and ground wires
    Lconds = len(xc.hot)
    _write_entries(ofile, Lconds)
    Lgrounds = len(xc.gnd)
    _write_entries(ofile, Lgrounds)
    #write the hot and gnd conductor data in the same format
    for c in xc.hot + xc.gnd:
        _write_entries(ofile, c.tag, c.x, c.y, c.subconds, c.d_cond, c.d_bund,
                'ED!(I)', c.I, c.V, c.phase)
    #write the ground wire data a second time, in a different format
    for c in xc.gnd:
        _write_entries(ofile, c.tag, c.x, c.y, c.d_cond, 0, 0)
    #close/save
    ofile.close()
    print('FLD file generated: %s' % fn)

def to_FLDs(*args, **kwargs):
    """Load or recieve a template workbook of CrossSections and convert them
    all to FLD files
    args:
        can either be a path string to a target template workbook or an
        existing SectionBook object
    kwargs:
        path - output destination for FLD files, default is the same directory
               as the template book if a path string is passed or the current
               directory if a SectionBook is passed."""
    if(type(args[0]) == str):
        #load the template
        sb = fields_funks.load_template(args[0])
        kwargs['path'] = os.path.dirname(args[0]) + '/'
    elif(isinstance(args[0], fields_class.SectionBook)):
        sb = args[0]
    else:
        raise(fields_class.EMFError("""
        Input argument to to_FLDs() must be a filepath or a SectionBook."""))
    #check for duplicate sheets
    sheets = []
    for xc in sb:
        if(xc.sheet in sheets):
            raise(fields_class.EMFError("""
            Can't create FLD files because of duplicate CrossSection names.
            Name "%s" is used at least twice.""" % xc.sheet))
        else:
            sheets.append(xc.sheet)
    #generate FLD files
    for xc in sb:
        to_FLD(xc, **kwargs)

def to_FLDs_crawl(dir_name, **kwargs):
    """crawl a directory and all of its subdirectories for excel workbooks
    that can be passed to create_FLDs(). The FLD files are generated in the
    same directory as the template book they come from."""

    #get input directory's file and subdir names
    dir_contents = glob.glob(dir_name)
    dir_name = dir_name.rstrip('\\/*').lstrip('\\/*')

    #loop over the dir_contents
    for dir_element in dir_contents:
        if(dir_element[-5:] == '.xlsx'):
            #operate on the file
            try:
                to_FLDs(dir_element, path = os.path.dirname(dir_element) + '/')
            except(KeyError, ValueError, IOError, fields_class.EMFError) as e:
                print('failure to write FLD files from:\n\t%s' % dir_element)
                if(type(e) is fields_class.EMFError):
                    print('\n\tBecause of EMFError:' + str(e))
        else:
            #if there's a period in the dir_element, it's not a directory
            if(not ('.' in dir_element)):
                to_FLDs_crawl(dir_element + '\\*')

#------------------------------------------------------------------------------
#FUNCTIONS FOR CONVERTING OUTPUT .DAT FILES TO CSV/EXCEL FILES

def read_DAT(file_path):
    """Read a DAT file, which can have some funky extra characters if the
    numbers are too large (percent signs)
    args:
        file_path - string, path to DAT file
    returns:
        pandas DataFrame with 'Bx','By','Bprod','Bmax','Ex','Ey','Eprod',
        and 'Emax' columns, and the distance ('x') in the index"""

    #check that the target file is a DAT
    fields_funks._check_extension(file_path, 'DAT', """
        Input file must have a '.DAT' extension.""")
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
    return(pd.DataFrame(data =
        dict(zip(['Bx','By','Bprod','Bmax','Ex','Ey','Eprod','Emax'],data[1:])),
        index = data[0]))

def convert_DAT(file_path, **kwargs):
    """read a DAT file and write it to a csv
    args:
        file_path - target DAT file
    kwargs:
        path - output location/name for csv file"""
    #get the DAT data
    df = read_DAT(file_path)
    #write it to a csv
    fn = fields_funks._path_manage(os.path.basename(file_path)[:file_path.index('.')],
        'csv', **kwargs)
    with open(fn, 'w') as ofile:
        df.to_csv(ofile, index_label = 'dist (ft)')
        print('DAT converted to csv: "%s"' % fn)

def convert_DAT_crawl(dir_name, **kwargs):
    """crawl a directory and all of its subdirectories for .DAT files that
    can be passed to DAT_to_csv() for output re-formatting and optional
    plotting.
    args:
        dir_name - the directory to initiate the crawl in. To designate the
                   current directory, use '*'
    kwargs:
        bundle - bool, if True, all DAT files found in the same directory are
                 written to a common excel workbook. If False or absent,
                 the DAT files are simply written to csv files."""

    #get input directory's file and subdir names
    dir_contents = glob.glob(dir_name)

    #get the bundling option
    if('bundle' in kwargs):
        bundle = kwargs['bundle']
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
                    fn = os.path.dirname(dir_element) + '/converted_DATs.xlsx'
                    xl = pd.ExcelWriter(fn, engine = 'xlsxwriter')
                #use the DAT's filename as the sheet name, without extension
                sn = os.path.basename(dir_element).replace('.DAT','')
                #read the DAT into a DataFrame and send it to the ExcelWriter
                read_DAT(dir_element).to_excel(xl, sheet_name = sn,
                    index_label = 'dist (ft)')
            else:
                #without bundling, write an individual csv file
                convert_DAT(dir_element, path = os.path.dirname(dir_name) + '/')
        else:
            #if there's a period in the dir_element, it's not a directory
            if(os.path.isdir(dir_element)):
                convert_DAT_crawl(dir_element + '/*', **kwargs)
    #close/save the excel bundles DAT results
    if(bundle and xl_flag):
        xl.save()
        print('Converted DAT files written to "%s"' % fn)
