def _str_Conductor(c):

    if(z._xs is not None):
        xs = repr(c._xs.sheet)
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
        f = ' '*6 + str(z.fields.iloc[[0, z.lROWi, z.rROWi, len(z.fields) - 1]])
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
      tags:             %s
      frequencies:      %s Hz
      x coordinates:    %s ft
      y coordinates:    %s ft
      subconductors:    %s
      diameters:        %s in
      bundle diameters: %s in
      voltages:         %s V
      currents:         %s A
      phase angles:     %s deg

    fields\n%s""" % (
        z.sheet, z.tag, z.title, z.soil_resistivity,
        z.max_dist, z.step,z.sample_height, z.lROW, z.rROW,
        ', '.join([str(i) for i in z.tags]),
        ', '.join([str(i) for i in z.freq]),
        ', '.join([str(i) for i in z.x]),
        ', '.join([str(i) for i in z.y]),
        ', '.join([str(i) for i in z.subconds]),
        ', '.join([str(i) for i in z.d_cond]),
        ', '.join([str(i) for i in z.d_bund]),
        ', '.join([str(i) for i in z.V]),
        ', '.join([str(i) for i in z.I]),
        ', '.join([str(i) for i in z.phase]), f))

def _str_SectionBook(z):

    b, sheet, ctag, v = z.complete
    if(b):
        f = ' '*6 + str(z.ROW_edge_max)
    else:
        f = """      ROW edge max fields unavailable.
      Attribute "%s"
        in Conductor "%s"
          in CrossSection "%s" is unset""" % (v[1:], ctag, sheet)

    return("""  CrossSection object
    name: %s
    sheets: %s
    unique tags: %s

    fields at CrossSection ROW edges:\n%s""" %
    (z.name, z.sheets, z.tags, f))
