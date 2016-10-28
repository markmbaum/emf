import os
import glob

#Flatten a list of lists. Only works on a list and the first level of
#sublists, not recursively on all levels
flatten = lambda L: [item for sublist in L for item in sublist]

def get_docstring(lines, line_number):
    i = line_number
    while((i < len(lines)) and ('"""' not in lines[i])):
        i += 1
    if(i < len(lines)):
        docstring = lines[i][lines[i].index('"""') + 3:]
        while('"""' not in docstring):
            i += 1
            docstring += lines[i]
        return(docstring.replace('"""', ''))
    else:
        return('')

def is_unindented(line):
    """Check if a line of text is unindented. Empty strings return False."""
    if(len(line.lstrip()) != len(line)):
        return(True)
    else:
        return(False)

def is_funk_def(line):
    if('def' in line):
         if(':' in line):
            if('(' in line):
                if(')' in line):
                    return(True)
    return(False)

def doc_funk(fn, lines, line_number):

    L = lines
    i = line_number

    if(not is_funk_def(L[i])):
        return(None)

    F = {'funkname': L[i][:L[i].index('(')].replace('def', '').strip(),
            'docstring': get_docstring(L, i),
            'filename': fn,
            'call': L[i].replace('def', '').replace(':', '').strip()}

    F['filename'] = os.path.abspath(os.path.join(__file__, F['filename']))

    return(F)

def funk_to_html(F):

    tab = '&nbsp; '*2
    ds = F['docstring'].replace('\n','<br>%s' % (tab*2))
    ds = ds.replace('    ','').replace('\t','')

    return(
    """<p class="funkname">%s</p>
        <p class="funk-call">%s %s</p>
        <p class="source-file">%s %s</p>

        <p class="funk-docstring">
            %s %s
        </p>""" % (
        F['funkname'], tab, F['call'], tab, F['filename'], tab*2, ds))



def is_class_def(line):
    if('class ' in line):
        if(':' in line):
            return(True)
    return(False)

def doc_class(fn, lines, line_number):

    L = lines
    i = line_number

    if(not is_class_def(L[i])):
        return(None)

    C = {'classname': L[i][:L[i].index(':')].replace('class', '').strip(),
            'inheritance': '',
            'docstring': get_docstring(L, i),
            'filename': fn,
            'constructor': '',
            'methods': []}

    if('(' in L[i]):
        C['classname'] = L[i][:L[i].index('(')].replace('class', '').strip()
        C['inheritance'] = L[i][L[i].index('(')+1:L[i].index(')')]
    else:
        C['classname'] = L[i][:L[i].index(':')].replace('class', '').strip()

    i += 1
    while((i < len(L)) and is_unindented(L[i])):
        if(is_funk_def(L[i])):
            F = doc_funk(fn, L, i)
            F['call'] = F['call'].replace('self, ', '').replace('self,', '')
            if(F['funkname'] == '__init__'):
                F['funkname'] = C['classname'] + ' constructor'
                F['call'] = F['call'].replace('__init__', C['classname'])
                C['constructor'] = F
            else:
                C['methods'].append(F)
        i += 1

    return(C)

def class_to_html(C):
    pass

fns = (glob.glob(os.path.join('..', 'emf', '*'))
        + glob.glob(os.path.join('..', 'emf', 'fields', '*')))

classes = []
funks = []

for fn in fns:

    if(fn[-3:] == '.py'):
        ifile = open(fn, 'r')
        in_class_def = False
        lines = ifile.readlines()
        for i in range(len(lines)):

            line = lines[i]

            if(is_class_def(line)):
                in_class_def = True
                classes.append(doc_class(fn, lines, i))
            elif(is_unindented(line)):
                in_class_def = False

            if(not in_class_def):
                if(is_funk_def(line)):
                    funks.append(doc_funk(fn, lines, i))

        ifile.close()

with open('test.txt', 'w') as ofile:
    ofile.write(funk_to_html(funks[-5]))
