import inspect

def clean_html(s):
    return(s.replace('\n', '').replace('\t', ' ').replace('    ', ' '))

def get_public_vars(obj):
    v = vars(obj)
    v = dict(zip([k for k in v if (k[0] != '_')],
                [v[k] for k in v if (k[0] != '_')]))
    return(v)

def nav_html(mod_name, submodules, functions, classes, other):
    """generate an unordered list of links for the navigation menu or index"""
    #wrap in a single dict for iteration
    keys = ['submodules', 'functions', 'classes', 'other']
    d = dict(zip(keys, [submodules, functions, classes, other]))
    html = '<div id="nav"><h2>Contents</h2><ul id="nav-ul">'
    for k in keys:
        if(d[k]):
            html += """
            <li class="nav-li"><a href="#%s">%s</a>
                <ul class="subnav-ul">""" % (k, k)
            for m in sorted(d[k].keys()):
                if(k == 'submodules'):
                    href = '%s-%s.html' % (mod_name, m)
                else:
                    href = '#%s' % m
                html += '<li class="subnav-li"><a href="%s">%s</a></li>' % (href, m)
            html += '</ul>'
    html += '</ul></div>'
    return(clean_html(html))

def get_docstr(obj):
    """Pull and format the docstring of an object for html"""
    tab = '&ensp;'*2
    docstr = obj.__doc__
    if(docstr):
        docstr = docstr.strip()
        docstr = docstr.replace('\n', '<br>'
                ).replace('\t', tab).replace('    ', tab)
        return('<p class="docstr">%s</p>' % docstr)
    else:
        return('')

def get_argstr(function):
    #get the argspec named tuple
    argspec = inspect.getargspec(function)
    #remove 'self'
    if('self' in argspec[0]):
        argspec[0].pop(argspec[0].index('self'))
    #apply spans for syntax coloring
    for i in range(len(argspec[0])):
        argspec[0] = '<span class="function-arg">' + argspec[0] + '</span>'
    #join default arg names with their values
    if(argspec[3] is not None):
        for i in range(1, len(argspec[3]) + 1):
            argspec[0][-i] += '=%s' % str(argspec[3][-i])
    #concatenate the complete string
    argstr = '('
    if(argspec[0]):
        argstr += ', '.join(argspec[0])
    if(argspec[1]):
        if(argspec[0]):
            argstr += ', '
        argstr += '*' + argspec[1]
    if(argspec[2]):
        if(argspec[1] or argspec[0]):
            argstr += ', '
        argstr += '**' + argspec[2]
    argstr += ')'
    return(argstr)

def function_to_html(F, name):
    """Create an html div for a function"""
    #open a div
    html = '<div class="function" id="%s">' % name
    #add the function name and argstr
    html += '<p class="function-str"><span class="function-name">%s</span>%s</p>' % (name, get_argstr(F))
    #add the docstr
    html += get_docstr(F)
    #add function name
    html += '</div>'
    return(clean_html(html))

def functions_to_html(functions):
    """Creat an html div for a group of functions, alphabetically"""
    #start the html
    html = ''
    if(functions):
        html += '<div class="section" id="functions"><h2>Functions</h2>'
        for k in sorted(functions.keys()):
            html += function_to_html(functions[k], k)
        html += '</div>'
    return(clean_html(html))

def property_to_html(P, name, classname):
    html = '<div class="property" id="%s-%s">' % (name, classname)
    html += '<p class="property-name">%s</p>' % name
    html += get_docstr(P) + '</div>'
    return(clean_html(html))

def class_to_html(C, name):
    """Create an html string for a class/type object"""
    #get it's variables, properties, and methods
    v = get_public_vars(C)
    #open a div for this individual class
    html = '<div class="class" id="%s">' % name
    #add the class name
    html += '<h3>%s</h3>' % name
    #add the class docstr
    if(C.__doc__):
        html += get_docstr(C)
    #add the class-constructor
    if('__init__' in vars(C)):
        init = vars(C)['__init__']
        html += function_to_html(init, name).replace('id="%s"' % name, '')
    #allocate separate strings for properties, methods
    properties, methods = '', ''
    for k in sorted(v.keys()):
        #check if property
        if(inspect.isdatadescriptor(v[k])):
            if(properties == ''):
                properties += '<div class="properties" id="%s-properties"><h4>Properties</h4>' % name
            properties += property_to_html(v[k], k, name)
        #check if method
        if(inspect.ismethod(v[k]) or (inspect.isfunction(v[k]))):
            if(methods == ''):
                methods += '<div class="methods"><h4>Methods</h4>'
            methods += function_to_html(v[k], k).replace(
                    'id="%s"' % k, 'id=%s-%s' % (name, k))
    if(properties != ''):
        properties += '</div>'
    if(methods != ''):
        methods += '</div>'
    html += properties + methods
    #close the individual class's div
    html += '</div>'
    return(clean_html(html))

def classes_to_html(classes):
    """Convert a dictionary of class/type objects into html documentation strings
    args:
        classes - dict mapping class/type names to class/type objects"""
    #start html string
    html = ''
    #proceed alphabetically
    if(classes):
        #section div for all classes
        html = '<div class="section" id="classes"><h2>Classes</h2>'
        for k in sorted(classes.keys()):
            html += class_to_html(classes[k], k)
        html += '</div>'
    return(clean_html(html))

def doc_module(mod, **kw):
    '''document a module by providing its name in the form of a string or a module object
    args:
        mod - str, module name or module object
    kw:
        all_submods - bool, determines whether submodules that are not in the module's __all__ list are recursively documented'''

    if('all_submods' in kw):
        all_submods = kw['all_submods']
    else:
        all_submods = False

    #import the module if a string is passed in
    if(type(mod) is str):
        exec('import %s as mod' % mod)
    #store the module name
    mod_name = mod.__name__

    #get the module's dictionary/namespace/whatever
    v = vars(mod)
    #check for submods in __all__
    if('__all__' in v):
        submods_in_all = v['__all__']
        for submod_name in submods_in_all:
            doc_module(v['%s' % submod_name])
    else:
        submods_in_all = []
    #remove private keys
    v = get_public_vars(mod)
    #separate the objects in v based on their types
    submodules, functions, classes, other = dict(), dict(), dict(), dict()
    for k in v:
        t = type(v[k])
        if(inspect.ismodule(v[k])):
            if(all_submods or (k in submods_in_all)):
                submodules[k] = v[k]
        elif(inspect.isfunction(v[k])):
            functions[k] = v[k]
        elif(inspect.isclass(v[k])):
            classes[k] = v[k]
        else:
            other[k] = v[k]
    #generate html
    html = ("""
    <html lang="en">
    	<head>
    		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    		<title>%s</title>
            <link href="https://fonts.googleapis.com/css?family=Anonymous+Pro|Open+Sans" rel="stylesheet">
            <link rel="stylesheet" type="text/css" href="style.css">
    	</head>
        <body onresize="arrange()" onload="arrange()">

            <!--Navigation Bar-->
            %s

            <div id="main">
                <h1>%s</h1>

                <!--classes section-->
                %s

                <!--functions section-->
                %s
            </div>

        </body>

        <!--Control the width and height of the #main element-->
        <script type="text/javascript">

            function arrange(){

                var nav = document.getElementById('nav');
                var main = document.getElementById('main');

                var w_nav = nav.offsetWidth;
                var w_body = document.body.offsetWidth;
                var w = (w_body - w_nav - 3).toString() + 'px';
                main.style.width = w;

                var h = main.offsetHeight;
                nav.style.height = h;
                //console.log('arrange');
                };

            arrange();

        </script>

    </html>
    """ %
    (mod_name,
    nav_html(mod_name, submodules, functions, classes, other),
    mod_name,
    classes_to_html(classes),
    functions_to_html(functions)))

    with open(mod_name.replace('.', '-') + '.html', 'w') as ofile:
        ofile.write(html)

    return(html)

doc_module('emf')
