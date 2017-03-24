"""
Like emf.fields, emf.subcalc was originally intended to supplement/extend the use of another modeling program. In this case, the other modeling program is SUBCALC (developed by Enertech, sponsored by EPRI). SUBCALC predicts magnetic fields over a fixed-height 2 dimensional grid and can model many non-parallel segments of power lines. It is generally used to model magnetic fields near electrical substations, where circuit arrangements can be complex. Just like emf.fields, emf.subcalc is now a complete replacement of the program it was originally intended to supplement. Although it provides no fancy UI, emf.subcalc is capable of calculating the magnetic fields produced by complex, three dimensional arrangements of current carrying wire segments.

While the emf.subcalc package is capable of calculating the magnetic fields generated by complex, three dimensional arrangements of current carrying wire segments, it is not capable of calculating fields generated by other electrical equiment like buses, circuit breakers, power transformers, air core reactors, and capacitor banks. The SUBCALC program, which emf.subcalc emulates and interfaces with, is capable of modeling that equipment. The SUBCALC program also provides a UI enabling a user to draw towers onto a map view of a model domain and edit the wire configurations of each tower. The UI is generally helpful for creating small models, but it can be frustrating when creating more complicated models and when creating slightly different versions of the same model. emf.subcalc provides no UI and would certainly be harder to learn, but it can perform the calculations that most modeling scenarios require and it provides better analysis and plotting methods than SUBCALC.

Importantly though, the old SUBCALC program is not free. It's not terribly expensive either, but one must navigate an annoying web of contacts and phone numbers to get an annually renewed license. In many tests, emf.subcalc has produced the exact same results as SUBCALC for the same inputs.
"""

from subcalc_class import (Model,
                        Tower,
                        Conductor,
                        Results,
                        Footprint)

from subcalc_funks import (drop_tower_template,
                        drop_footprint_template,
                        mesh_dict_grids,
                        load_results,
                        load_towers,
                        cumulative_distance)

from subcalc_calcs import (B_field_segment,
                        grid_segment_results,
                        phasors_to_magnitudes)

from subcalc_plots import (plot_cross_sections,
                        plot_contour,
                        plot_pcolormesh,
                        plot_path,
                        plot_wires_3D,
                        ion,
                        show,
                        close,
                        be_concerned)

from SUBCALC_io import (read_INP,
                        convert_INP,
                        read_REF,
                        convert_REF)
