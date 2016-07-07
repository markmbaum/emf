import os
import glob
import numpy as np
import pandas as pd

import emf_funks
import emf_class

#-------------------------------------------------------------------------------
#FUNCTIONS FOR GENERATING INPUT .FLD FILES

def is_int(x):
    """check if a number is an integer, will break on non-numeric entries"""
    if(int(x) == x):
        return True
    else:
        return False

#function formats entries and writes them into a .FLD file targeted by ofile
def write_entry(ofile, entry):
    """format entries and writes them into a .FLD file targeted by ofile"""
    if(emf_funks.is_number(entry)):
        if(is_int(entry)):
            w = '{:< d} \n'.format(int(entry))
        else:
            w = '{:< .2f} \n'.format(float(entry))
        if('.' in w):
            idx = w.index('.')
            if('0' in w[:idx]):
                idx = w.index('0')
                if(not emf_funks.is_number(w[idx-1])):
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
    #get a filename
    fn = emf_funks.path_manage(xc.title, 'FLD', **kwargs)
    #write the .FLD file
    ofile = open(fn, 'w')
    #miscellaneous stuff first
    write_entry(ofile, xc.title)
    write_entry(ofile, xc.subtitle)
    write_entry(ofile, 60)
    write_entry(ofile, xc.soil_resistivity)
    write_entry(ofile, xc.max_dist)
    write_entry(ofile, xc.step)
    write_entry(ofile, xc.sample_height)
    write_entry(ofile, xc.lROW)
    write_entry(ofile, xc.rROW)
    #number of conductors and ground wires
    Lconds = len(xc.hot)
    write_entry(ofile, Lconds)
    Lgrounds = len(xc.gnd)
    write_entry(ofile, Lgrounds)
    #write the hot and gnd conductor data in the same format
    for c in xc.hot + xc.gnd:
        write_entry(ofile, c.tag)
        write_entry(ofile, c.x)
        write_entry(ofile, c.y)
        write_entry(ofile, c.subconds)
        write_entry(ofile, c.d_cond)
        write_entry(ofile, c.d_bund)
        write_entry(ofile, 'ED!(I)')
        write_entry(ofile, c.I)
        write_entry(ofile, c.V)
        write_entry(ofile, c.phase)
    #write the ground wire data a second time, in a different format
    for c in xc.gnd:
        write_entry(ofile, c.tag)
        write_entry(ofile, c.x)
        write_entry(ofile, c.y)
        write_entry(ofile, c.d_cond)
    #close/save
    ofile.close()
    print('FLD file generated: "%s"' % fn)

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
        sb = emf_funks.load_template(args[0])
        kwargs['path'] = os.path.dirname(args[0]) + '/'
    else:
        sb = args[0]
    #check for duplicate titles and subtitles
    titles, subtitles = [], []
    for xc in sb:
        if(xc.title in titles):
            raise(emf_class.EMFError("""
            Can't create FLD files because of duplicate CrossSection titles.
            Title "%s" is used at least twice.""" % xc.title))
        else:
            titles.append(xc.title)
        if(xc.subtitle in subtitles):
            raise(emf_class.EMFError("""
            Can't create FLD files because of duplicate CrossSection subtitles.
            Subtitle "%s" is used at least twice.""" % xc.subtitle))
        else:
            subtitles.append(xc.subtitle)
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

    #loop over the dir_contents, extracting if .LST is at the end and attempting
    #find subdirectories if not
    for dir_element in dir_contents:
        if(dir_element[-5:] == '.xlsx'):
            #operate on the file
            try:
                to_FLDs(dir_element, path = os.path.dirname(dir_element) + '/')
            except(KeyError, ValueError, IOError):
                print('failure to write FLD files from:\n\t%s' % dir_element)
        else:
            #if there's a period in the dir_element, it's not a directory
            if(not ('.' in dir_element)):
                to_FLDs_crawl(dir_element + '\\*')

#------------------------------------------------------------------------------
#FUNCTIONS FOR CONVERTING OUTPUT .DAT FILES TO CSV FILE AND PLOTTING

def read_DAT(file_path):
    """Read a DAT file, which can have some funky extra characters if the
    numbers are too large (percent signs)"""

    #check that the target file is a DAT
    emf_funks.check_extention(file_path, 'DAT', """
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
        x,Bx,By,Bprod,Bmax,Ex,Ey,Eprod,Emax = [],[],[],[],[],[],[],[],[]
        line = ifile.readline()
        while(line):
            if(line):
                if(line[0][0] == '%'):
                    line += ifile.readline()
            temp = [i.replace('%','') for i in line.split()]
            if(temp and all([emf_funks.is_number(i) for i in temp])):
                x.append(float(temp[0]))
                Bx.append(float(temp[1]))
                By.append(float(temp[2]))
                Bprod.append(float(temp[3]))
                Bmax.append(float(temp[4]))
                if(not (und_only)):
                    Ex.append(float(temp[5]))
                    Ey.append(float(temp[6]))
                    Eprod.append(float(temp[7]))
                    Emax.append(float(temp[8]))
                else:
                    Ex.append(0.)
                    Ey.append(0.)
                    Eprod.append(0.)
                    Emax.append(0.)
            line = ifile.readline()
    return(pd.DataFrame(data = {
        'Ex':Ex,'Ey':Ey,'Eprod':Eprod,'Emax':Emax,
        'Bx':Bx,'By':By,'Bprod':Bprod,'Bmax':Bmax}, index = x))

def convert_DAT(file_path, **kwargs):
    """read a DAT file and write it to a csv
    args:
        file_path - target DAT file
    kwargs:
        path - output location/name for csv file"""
    #get the DAT data
    df = read_DAT(file_path)
    #write it to a csv
    fn = emf_funks.path_manage(os.path.basename(file_path)[:file_path.index('.')],
        'csv', **kwargs)
    with open(fn, 'w') as ofile:
        df.to_csv(ofile)
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
                sn = os.path.basename(dir_element)[:-4]
                #read the DAT into a DataFrame and send it to the ExcelWriter
                read_DAT(dir_element).to_excel(xl, sheet_name = sn)
            else:
                #without bundling, write an individual csv file
                convert_DAT(dir_element, path = dir_name + '/')
        else:
            #if there's a period in the dir_element, it's not a directory
            if(os.path.isdir(dir_element)):
                convert_DAT_crawl(dir_element + '\\*', **kwargs)
    if(bundle and xl_flag):
        xl.save()
        print('Converted DAT files written to "%s"' % fn)
