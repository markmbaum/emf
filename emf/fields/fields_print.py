from .. import pd, textwrap

wrapper = textwrap.TextWrapper(width=70, expand_tabs=True)

def _table_iterable_fill(field, a):
    wrapper.subsequent_indent = ' '*len(field)
    v = []
    for i in a:
        if((type(i) is unicode) or (type(i) is str)):
            v.append(repr(i))
        else:
            v.append(str(i))
    return(wrapper.fill(field + ', '.join(v)))

def _str_Conductor(z):

    if(z._xs is not None):
        xs = repr(z._xs.sheet)
    else:
        xs = repr(None)

    return("""  Conductor object
    name:                    %s
    parent CrossSection:     %s
    frequency (Hz):          %s
    x coordinate (ft):       %s
    y coordinate (ft):       %s
    subconductors:           %s
    conductor diameter (in): %s
    bundle diameter (in):    %s
    voltage (V):             %s
    current (I):             %s
    phase angle (deg):       %s""" %
        (repr(z.name), xs, repr(z.freq), repr(z.x), repr(z.y),
        repr(z.subconds), repr(z.d_cond), repr(z.d_bund),
        repr(z.V), repr(z.I), repr(z.phase)))

def _str_CrossSection(z):

    b, t, v = z.complete
    if(b):
        idx = z.fields.index
        locs = sorted(list(set([idx.min(), z.lROW, z.rROW, idx.max()])))
        f = z.fields.loc[locs]
        f = pd.concat((f['Bmax'], f['Emax']), axis=1)
        f = ' '*6 + str(f).replace('\n', '\n' + ' '*6)
    else:
        f = """      Cannot compute fields.
      Attribute "%s"
        in Conductor "%s" is unset.""" % (t, v[1:])

    if(z._sb is not None):
        sb = repr(z._sb.name)
    else:
        sb = repr(None)

    return("""  CrossSection object

    sheet:                         %s
    parent SectionBook:            %s
    group:                         %s
    title:                         %s
    soil resistivity (?):          %s
    max distance from center (ft): %s
    step size (ft):                %s
    sample height (ft):            %s
    left ROW edge (ft):            %s
    right ROW edge (ft):           %s

    conductor information (%d Conductors)
%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s

    fields sample (see CrossSection.fields for all EMF results)\n%s""" % (
        z.sheet, sb, z.group, z.title, z.soil_resistivity,
        z.max_dist, z.step, z.sample_height, z.lROW, z.rROW, len(z.conds),
        _table_iterable_fill('      names:                 ', z.names),
        _table_iterable_fill('      frequencies (Hz):      ', z.freq),
        _table_iterable_fill('      x coordinates (ft):    ', z.x),
        _table_iterable_fill('      y coordinates (ft):    ', z.y),
        _table_iterable_fill('      subconductors:         ', z.subconds),
        _table_iterable_fill('      diameters (in):        ', z.d_cond),
        _table_iterable_fill('      bundle diameters (in): ', z.d_bund),
        _table_iterable_fill('      voltages (V):          ', z.V),
        _table_iterable_fill('      currents (A):          ', z.I),
        _table_iterable_fill('      phase angles (deg):    ', z.phase), f))

def _str_SectionBook(z):

    b, sheet, cname, v = z.complete
    if(b):
        f = ' '*6 + str(z.ROW_edge_max).replace('\n', '\n' + ' '*6)
    else:
        f = """      ROW edge max fields unavailable.
      Attribute "%s"
        in Conductor "%s"
          in CrossSection "%s" is unset""" % (v[1:], cname, sheet)

    return("""  SectionBook object
      name:         %s\n%s\n%s

    maximum fields at CrossSection ROW edges:\n%s""" %
    (repr(z.name),
    _table_iterable_fill('      sheets:        ', z.sheets),
    _table_iterable_fill('      unique groups: ', set([xs.group for xs in z.xss])), f))
