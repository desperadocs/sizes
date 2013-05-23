.. $Id$

SOURCE CODE
==================================================

This Maximum Entropy program was originally written in BASIC by G.J.
Daniell (Department of Physics, Southampton University, UK) and later
translated into FORTRAN and adapted for SAS analysis by J.A. Potton.
Further modifications have been made by I.D. Culverwell, G.P. Clarke
and A.J. Allen (UKAEA Harwell Laboratory,UK) and P.R. Jemian
(Northwestern University, USA).

There is only one source code module, MaxSas.For.  Compile and link it
with the fastest floating point math that you can get your hands on.

Unfortunately, some data storage had to be placed in COMMON because of
the limitation of the Language Systems MPW version 1.2.1 FORTRAN
compiler for the Apple Macintosh.  Because of this compiler's
eccentricacies, there is one compiler-dependent line of code very near
the first executable statement.  If you use this compiler, un-comment
this line so that you get a chance to see the output.  (Compiler
dependence, ugh!)

As it stands on 7 February 1990, the code will now compile on:

*    Digital Equipment Corporation VAX 11/785,VMS version 5.2
*    Apple Macintosh, Language Systems FORTRAN v. 1.2.1
*    Apple Macintosh, Microsoft (Absoft) FORTRAN v. 2.2 compiler
*    MS-DOS (e.g. IBM-PC), Microsoft FORTRAN v. 5.0

Of course the program RUNS on these computers as well.  Quite well!

Most of the comments in the source code have been added by P.R.
Jemian. Where they exist, they are usually quite explanatory.  Where
they do not exist, consult the references of [Skilling, 1984] for the
operation of MaxEnt.


Complete listing of **MaxSas.for**
--------------------------------------

.. literalinclude:: ../../../fortran/maxsas/maxsas.for
    :tab-width: 4
    :linenos:
    :language: guess

