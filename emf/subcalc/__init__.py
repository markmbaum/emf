"""emf.subcalc is currently similar to the original versions of emf.fields because it supplements another modeling program without fully replacing it. In this case, the other modeling program is SUBCALC (developed by Enertech, sponsored by EPRI). SUBCALC predicts EMF over a fixed-height 2 dimensional grid and can model many non-parallel segments of power lines. emf.subcalc's primary functions are to:

    - read the text file output of SUBCALC models and associate the results
      with Footprint objects (outlines of nearby objects in the model domain, like houses)

    - generate contour and colormesh plots of the results, with maximum fields
      annotated along Footprint objects

    - convert results to excel files (much smaller file size)

    - interpolate the grid (at points, along lines, or complete resampling)

It would be nice to build the calculations into emf.subcalc and totally replace SUBCALC, but that hasn't been done yet.
"""

from subcalc_class import (Model,
                        Footprint)

from subcalc_funks import (drop_footprint_template,
						read_REF,
                        load_model,
                        convert_REF)

from subcalc_plots import (plot_contour,
                        plot_pcolormesh,
                        ion,
                        show,
                        close)

del(subcalc_class, subcalc_funks, subcalc_plots)
