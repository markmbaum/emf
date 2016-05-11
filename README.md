# FLD

This repo contains code for the approximation of electric and magnetic fields
near parallel groups of power lines. It's meant to replace the old and
difficult to use program called FIELDS, originally developed by California
Edison, which is been the standard tool for these simple modeling efforts.
The FIELDS software runs through a DOSBOX application and is long out of date.

This code is intended to exactly replicate the results of FIELDS by following
the conceptual guidelines laid out in the Electric Power Research Institute's
"Big Red Book" (exact title/citation unknown at the moment). The FIELDS source
code is not released, so a line-by-line replication of the calculation isn't
possible. However, the physics are definitely manageable and this version of
the code has been able to reproduce FIELDS results to a very high degree of
accuracy. In most cases, this code is able to exactly match FIELDS results.
However, because FIELDS only prints results to three digits past the decimal
point, a question remains about the internal precision of the calculations.
The largest difference between the results of this code and those of 
FIELDS has been about 2% and difference are vanishingly small in most cases,
on the order of round-off error.

The FIELDS method of calculating EMF near transmission lines is not improved by
this code. It seems like the FIELDS approach has a lot of inertia and people
are hesitant to diverge from it. However, this code does significantly improve
upon the usability and analytical capabilities of FIELDS, mostly by making the
functions that calculate electric and magnetic fields accessible. It interfaces
with excel templates that store the conductor and cross section data. In the
most routine modeling scenarios, this code enables a one line effort (after
filling in template worksheets) to generate full sets of electric and magnetic
field results and accompanying double-axis plots of both fields.
