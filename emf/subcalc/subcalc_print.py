from ..emf_funks import _sig_figs

perc_g = lambda x: ('%g' % x)

def _element_fill_join(elements, width):
    """Create a multiline string with a maximum width from a list of strings,
    without breaking lines within the elements"""

    s = ''
    if(elements):
        L = 0
        for i in range(len(elements) - 1):
            s += str(elements[i]) + ', '
            L += len(elements[i]) + 2
            if(L + len(elements[i+1]) >= width - 2):
                s += '\n'
                L = 0
        s += elements[-1]

    return(s)

def _element_fill_join_indent(elements, indent, conv_funk):

    elements = [conv_funk(i) for i in elements]
    s = _element_fill_join(elements, 80 - indent).replace('\n', '\n' +  ' '*indent)
    return(s)

def _str_Model(mod):

    return(
        '\n    '.join(
            ['Model object',
            'name: %s' % repr(mod.name),
            'x limits: %g to %g ft' % (mod.xmin, mod.xmax),
            'y limits: %g to %g ft' % (mod.ymin, mod.ymax),
            'total samples: %d' % mod.N,
            'sample spacing: %g ft' % mod.spacing,
            'number of Tower objects: %d' % len(mod.towers),
            'number of Tower groups: %d' % len(mod.tower_group_names),
            'tower group names: %s' % _element_fill_join_indent(
                    mod.tower_group_names, 23, repr),
            'number of Condutor objects: %s' % len(mod.conductors),
            'total number of wire segments: %s' % len(mod.segments)]
        )
    )

def _str_Results(res):

    spacing = tuple([str(_sig_figs(i, 6)) for i in res.spacing])
    spacing = '%s ft along x axis, %s ft along y axis' % res.spacing

    return(
        '\n    '.join(
            ['Results object',
            'name: %s' % repr(res.name),
            'components/Bkeys: %s' % (', '.join([repr(i) for i in res.Bkeys])),
            '%s range: %g to %g mG' % (res.Bkey, res.Bmin, res.Bmax),
            'x limits: %g to %g ft' % (res.xmin, res.xmax),
            'y limits: %g to %g ft' % (res.ymin, res.ymax),
            'total samples: %d' % res.N,
            'sample spacing: %s' % spacing,
            'number of Footprints: %d' % len(res.footprints),
            'Footprint groups: %s' % _element_fill_join_indent(
                    res.footprint_group_names, 22, repr)]
        )
    )

def _str_Tower(t):

    return(
        '\n    '.join(
            ['Tower object',
            'group: %s' % repr(t.group),
            'sequence number: %d' % t.seq,
            'x coordinate: %g ft' % t.tower_x,
            'y coordinate: %g ft' % t.tower_y,
            'rotation: %g degrees' % t.rot,
            'wire distances from tower (ft): %s' % _element_fill_join_indent(
                    t.h, 36, perc_g),
            'wire distances from ground (ft): %s' % _element_fill_join_indent(
                    t.v, 37, perc_g),
            'wire x coordinates (ft): %s' % _element_fill_join_indent(
                    t.conductor_x, 25, perc_g),
            'wire y coordinates (ft): %s' % _element_fill_join_indent(
                    t.conductor_y, 25, perc_g),
            'currents (Amps): %s' % _element_fill_join_indent(
                    t.I, 17, perc_g),
            'phase angles (degrees): %s' % _element_fill_join_indent(
                    t.phase, 24, perc_g)]
        )
    )

def _str_Conductor(c):

    return(
        '\n    '.join(
            ['Conductor object',
            'name: %s' % repr(c.name),
            'x (ft): %s' % _element_fill_join_indent(c.x, 12, perc_g),
            'y (ft): %s' % _element_fill_join_indent(c.y, 12, perc_g),
            'current: %g Amps' % c.I,
            'phase angle: %g degrees' % c.phase]
        )
    )

def _str_Footprint(f):

    return(
        '\n    '.join(
            ['Footprint object',
            'name: %s' % repr(f.name),
            'group: %s' % repr(f.group),
            'x (ft): %s' % _element_fill_join_indent(f.x, 12, perc_g),
            'y (ft): %s' % _element_fill_join_indent(f.y, 12, perc_g),
            'power line?: %s' % str(f.power_line),
            'of concern?: %s' % str(f.of_concern),
            'draw as loop?: %s' % str(f.draw_as_loop)]
        )
    )
