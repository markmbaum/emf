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

    ds = ''
    for line in F['docstring'].split('\n'):
        line = line.replace('    ', '\t')

        if(' - ' in line):
            line = line.split('-')
            left, right = line[0], ''.join(line[1:])
            ds += """
                    <div class="docstr-split-div">
                        <div class="docstr-left" style="width=%s">
                            %s -&nbsp;
                        </div>
                        <div class="docstr-right">
                            %s
                        </div>
                    </div>""" % (
                        '%drem' % (len(left) + 3),
                        left,
                        right.replace('\t', ''))
        else:
            line = line.replace('args:', '<span class="args">arguments</span>:')
            line = line.replace('kw:', '<span class="kw">keyword arguments</span>:')
            line = line.replace('returns:', '<span class="returns variables">returns</span>:')
            indents = 0
            if(':' in line):
                while(line[indents] == '\t'):
                    indents += 1
            ds += ('<p class="funk-docstr-line" style="margin-left: %dcm">%s</p>' %
                    (indents, line))

    return(
    ("""<div class="funk" id="%s">
            <p class="funk-name">%s<span class="funk-call">%s</span></p>
            <div class="funk-docstr">%s</div>
        </div>""" % (
        F['funkname'], F['funkname'], F['call'].replace(F['funkname'],''), ds)
            ).replace('\n', '').replace('    ', '').replace('\t',''))

def is_property_def(line):
    if('property' in line):
        if('=' in line):
            if('(' in line):
                if(')' in line):
                    return(True)
    return(False)

def doc_property(lines, line_number):

    L = lines
    i = line_number

    propstr = L[i]
    while(')' not in L[i]):
        propstr += L[i]
        i += 1
    propstr += L[i]

    docstring = ''
    if(propstr.count("'") > 1):
        docstring = propstr[:propstr.rfind("'")]
        docstring = docstring[docstring.rfind("'")+1:]

    P = {'propname': propstr.split('=')[0].strip(),
        'docstring': docstring}

    return(P)

def property_to_html(P):

    return("""
    <div class="docstr-split-div">
        <div class="docstr-left" style="width=%s">
            %s -&nbsp;
        </div>
        <div class="docstr-right">
            %s
        </div>
    </div>""" % (
        '%drem' % len(P['propname']),
        P['propname'],
        P['docstring']))

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
            'methods': [],
            'properties': []}

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
                F['funkname'] = C['classname']
                F['call'] = F['call'].replace('__init__', C['classname'])
                C['constructor'] = F
            else:
                C['methods'].append(F)

        if(is_property_def(L[i])):
            C['properties'].append(doc_property(L, i))

        i += 1

    return(C)

def class_to_html(C):

    props = sorted(C['properties'], key=lambda x: x['propname'])
    methods = sorted(C['methods'], key=lambda x: x['funkname'])
    methods = [i for i in methods if i['funkname'][0] != '_']

    return(("""
    <div class="class" id="%s">
        <p class="class-name">%s</p>
        <p class="class-inheritance">%s</p>
        %s
        <div class="class-properties">
            %s
        </div>
        <div class="class-methods">
            %s
        </div>
    </div>""" %
        (C['classname'],
        C['classname'],
        C['inheritance'],
        funk_to_html(C['constructor']),
        ''.join([property_to_html(i) for i in props]),
        ''.join([funk_to_html(i) for i in methods])
        )).replace('\n', '').replace('    ', '').replace('\t',''))


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

with open('text.txt', 'w') as ofile:
    ofile.write(class_to_html(classes[2]))
