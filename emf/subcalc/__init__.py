"""
emf.subcalc is currently similar to the original versions of emf.fields because it supplements another modeling program without fully replacing it yet. In this case, the other modeling program is SUBCALC (developed by Enertech, sponsored by EPRI). SUBCALC predicts EMF over a fixed-height 2 dimensional grid and can model many non-parallel segments of power lines. It is generally used to model EMF near electrical substations, where circuit arrangements can be complex. The functions needed to fully replace SUBCALC have been written and tested, but there is not currently a framework for building models and using them efficiently. Alas, the primary functions of current versions of emf.subcalc are to:

  * read the text file output of SUBCALC models and associate the results with Footprint objects (outlines of nearby objects in the model domain, like houses)

  * generate contour and colormesh plots of the results with some options for annotations like maximum field markers, footprint outlines, and cross sections

  * convert results to excel files (much smaller file size and easier to work with than formatted text)

  * interpolate the results grid (at points, along line segments, along paths, or complete resampling)

  * manipulate the results grid (zoom, rereference)

As mentioned above, the code for fully replacing SUBCALC and giving emf.subcalc the tools to create standalone models is well underway.
"""

from subcalc_class import (Model,
                        Tower,
                        Conductor,
                        Results,
                        Footprint)

from subcalc_funks import (drop_tower_template,
                        drop_footprint_template,
                        load_towers,
                        read_REF,
                        mesh_dict_grids,
                        load_results,
                        convert_REF,
                        cumulative_distance)

from subcalc_calcs import (B_field_segment,
                        grid_segment_results,
                        phasors_to_magnitudes)

from subcalc_plots import (plot_cross_sections,
                        plot_contour,
                        plot_pcolormesh,
                        plot_path,
                        ion,
                        show,
                        close,
                        be_concerned)
