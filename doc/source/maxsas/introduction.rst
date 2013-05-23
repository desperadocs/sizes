.. $Id$

INTRODUCTION
============

This program provides an accurate and reliable analysis of small-angle
scattering data.  Using an iterative approach it calculates the size
distribution of scatterers of a specified form, with the maximum
entropy, whose scattering pattern is consistent with the data.  The
consistency test is the ChiSquared statistic.  A full description of
the program's methodology and of its rigorous validation can be found
in [Culverwell, 1986] and [Potton, 1988a].

The MaxSas code is an interactive program, and users must be prepared
to answer the questions that the program will ask.  The questions are
grouped into two major sections: input data selection and scattering
model selection.  The input data selection is by file name and
Q-vector range. The scattering model selection is by form factor,
diameter range, and number of histogram bins.  It may be necessary to
run this program several times, with slight changes in user-adjustable
parameters, in order to obtain a satisfactory analysis.

The running of the program will be demonstrated by example, using
synthetic scattering data calculated from a hypothetical bimodal
particle size distribution (ref 3).   This distribution consists of a
narrow Gaussian peak centered at 80 Angstroms with a standard
deviation of 20 A, together with a broader secondary Gaussian peak,
half the size of the primary one, centered at 200 A with a width
of 60 A.   In order that they can become familiar with the questions
that the program will ask them, prospective users are advised to run
the program using this data set, reproducing verbatim the example
responses quoted in this guide.
