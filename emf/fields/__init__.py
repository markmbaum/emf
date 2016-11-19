"""
emf.fields was originally a small stand-alone package that streamlined the use of an old electromagnetic field (EMF) modeling program called FIELDS, which predicts electric and magnetic fields near parallel sets of power lines by assuming the conductors are infinitely long and computing the fields along a transect perpendicular to the power lines (a cross section model). The old FIELDS program is very difficult to use (details on why below), so it made sense to transplant as much of the modeling process as possible into other programs. Initially, to compute a cross section model, emf.fields would read power line information from excel templates, create input files that could be run through FIELDS to get EMF results, then read the FIELDS output files and provide plots, formatted results, etc. The package did everything except perform the actual EMF calculations. emf.fields still contains functions streamlining the use of FIELDS, but at this point the package does all the things the old FIELDS program does and more, including the EMF calculations, and it does them better.
"""

from fields_class import (Conductor,
                        CrossSection,
                        SectionBook)

from fields_funks import (drop_template,
						load_template,
                        optimize_phasing,
                        target_fields,
                        run)

from fields_calcs import (E_field,
                        B_field,
                        phasors_to_magnitudes,
                        EPSILON,
                        electric_prefactor,
                        MU,
                        magnetic_prefactor)

from fields_plots import (plot_Bmax,
                        plot_Emax,
                        plot_max_fields,
                        plot_xs,
                        plot_groups,
                        plot_groups_at_ROW,
                        ion,
                        show,
                        close)

from FIELDS_io import (to_FLD,
                        to_FLDs,
                        to_FLDs_crawl,
                        read_DAT,
                        convert_DAT,
                        convert_DAT_crawl)

del(fields_class,
    fields_funks,
    fields_calcs,
    fields_plots,
    FIELDS_io,
    fields_print)
