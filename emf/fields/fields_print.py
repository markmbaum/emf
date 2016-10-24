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
    tag:                 %s
    parent CrossSection: %s
    frequency:           %s Hz
    x coordinate:        %s ft
    y coordinate:        %s ft
    subconductors:       %s
    conductor diameter:  %s in
    bundle diameter:     %s in
    voltage:             %s V
    current:             %s A
    phase angle:         %s deg""" %
        (repr(z.tag), xs, repr(z.freq), repr(z.x), repr(z.y),
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

    return("""  CrossSection object

    sheet:                    %s
    tag:                      %s
    title:                    %s
    soil resistivity:         %s ?
    max distance from center: %s ft
    step size:                %s ft
    sample height:            %s ft
    left ROW edge:            %s ft
    right ROW edge:           %s ft

    conductor information
%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s

    fields sample (see CrossSection.fields for all EMF results)\n%s""" % (
        z.sheet, z.tag, z.title, z.soil_resistivity,
        z.max_dist, z.step, z.sample_height, z.lROW, z.rROW,
        _table_iterable_fill('      tags:                  ', z.tags),
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

    b, sheet, ctag, v = z.complete
    if(b):
        f = ' '*6 + str(z.ROW_edge_max).replace('\n', '\n' + ' '*6)
    else:
        f = """      ROW edge max fields unavailable.
      Attribute "%s"
        in Conductor "%s"
          in CrossSection "%s" is unset""" % (v[1:], ctag, sheet)

    return("""  SectionBook object
      name:        %s\n%s\n%s

    fields at CrossSection ROW edges:\n%s""" %
	(repr(z.name),
	_table_iterable_fill('      sheets:      ', z.sheets),
	_table_iterable_fill('      unique tags: ', z.tags), f))
