# FLD

This repo contains code for the approximation of electric and magnetic fields
near parallel groups of power lines. It's meant to replace the old and
difficult program called FIELDS, originally developed by Southern California
Edison Co., which has been the standard tool for these simple modeling efforts.
The FIELDS software runs through a DOSBOX application and is long out of date.

This code is intended to exactly replicate the results of FIELDS by following
the conceptual guidelines laid out in the Electric Power Research Institute's
"Big Red Book" (exact title/citation unknown at the moment). The FIELDS source
code is not released, so a line-by-line replication of the calculations isn't
possible. However, the physics are pretty manageable and this version of
the code has been able to reproduce FIELDS results to a very high degree of
accuracy. In most cases, this code is able to exactly match FIELDS results.
The largest difference between the results of this code and those of
FIELDS has been about 2%. In most cases the error is vanishingly small,
on the order of round-off error.

The FIELDS method of calculating EMF near transmission lines is not improved by
this code. It seems like the FIELDS approach has a lot of inertia and people
are hesitant to diverge from it. However, this code does significantly improve
upon the usability and analytical capabilities of FIELDS, mostly by making the
functions that calculate electric and magnetic fields accessible. It uses data
structures and I/O methods from the pandas library to interface with excel
templates that store the conductor and cross section data. It implements
three classes to organize the imported data and the EMF results in a flexible,
hierarchical system. For the most routine modeling scenarios, this code enables
a one line effort (after filling in template excel sheets) to generate full sets
of electric and magnetic field results, double-axis plots of both electric and
magnetic fields, plots comparing the electric and magnetic fields of grouped
cross sections, and a table of maximum field magnitudes at the right-of-way
edges of each cross-section. The run() function does all this and only requires
the path of the excel workbook of templates.
