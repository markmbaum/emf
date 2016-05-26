# FLD

### Computing Electric and Magnetic Fields Near Parallel Sets of High Voltage Transmission Lines

This repo contains code for the approximation of electric and magnetic fields
near parallel groups of power lines. It's meant to replace the old and
difficult program called FIELDS, originally developed by Southern California
Edison Co., which has been the standard tool for these simple modeling efforts.
The FIELDS software runs through a DOSBOX application and is long out of date.

This code is intended to exactly replicate the results of FIELDS by following
the conceptual guidelines laid out in the Electric Power Research Institute's
"Red Book" (some more information on this source below). The FIELDS source
code isn't released, so a line-by-line replication of the calculations isn't
possible. However, this version of the code has been able to reproduce FIELDS results to a very high degree of accuracy.

Some testing has shown error between this code and FIELDS results up to 3 %. I've been told that the FIELDS program runs 16-bit BASIC, lower precision than modern languages, and that BASIC has known accuracy issues with the sine and cosine functions. Additionally, the results of FIELDS simulations are saved to output files with the values rounded or truncated to the thousandths digit. These factors might explain the error. However, most testing has shown extremely small error, on the order of actual computational roundoff.

The FIELDS method of calculating EMF near transmission lines is not improved by
this code. It seems like the FIELDS approach has a lot of inertia and people
are hesitant to diverge from it. However, this code does significantly improve
upon the usability and analytical capabilities of FIELDS, mostly by making the
functions that calculate electric and magnetic fields accessible. This code relies on:
* [data structures](http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe)
and [I/O methods](http://pandas.pydata.org/pandas-docs/stable/io.html) from the [pandas library](http://pandas.pydata.org/pandas-docs/stable/index.html) to interface with excel templates, store the results of fields simulations, and write results to output files
* fast, explicit, low-level [numpy](http://www.numpy.org/) arrays and functions to perform the actual fields calculations
* the [matplotlib](http://matplotlib.org/) plotting package to automatically generate useful plots of the simulation results
* three custom classes (Conductor, CrossSection, and SectionBook) to organize the imported data and the EMF results in a flexible, hierarchical, system that has room for expansion.

For the most routine modeling scenarios, this code enables a one line effort (after filling in template excel sheets) to generate full sets of electric and magnetic field results, double-axis plots of both electric and
magnetic fields, plots comparing the electric and magnetic fields of grouped
cross sections, and a table of maximum field magnitudes at the right-of-way
edges of each cross-section. The run() function does all this and only requires
the path of the excel workbook of templates.

######EPRI's "Red Book"

The editon of EPRI's "Red Book" that I worked from is titled
"Transmission Line Reference Book: 345 kV and Above/Second Edition." Section
8.3 outlines the calculation of electric fields. Section 8.4 outlines the
calculation of magnetic fields. Appendix 8.1 details how to calculate the
maximum field magnitude from horizontal and vertical component phasors, which
might sound almost trivial but is more involved than I thought it would be.
