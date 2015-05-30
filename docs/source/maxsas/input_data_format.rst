INPUT DATA FORMAT
==================

The input data file should be in the same format as in the supplied
synthetic scattering data (file: BIMODAL.SAS).  Formally the data
should exist in the file as ordered triples of scattering vector (in
1/A), intensity and estimated error of intensity.  Typically, the data
will be in three columns, separated by "white space" of spaces and/or
tabs.  The intensity may be in any units as the program will ask for a
conversion factor between these units and 1/cm units.  The numbers are
read as floating point numbers, hence 0.0015, 1.5E-3, 0.15E-2, and
15E-4 will be identical.  However, do not be tempted to force double
precision input (such as 1.5D-3) because this may provoke a nasty
error condition.

The name of the data set should any standard filename for the
operating system on which the program is run.  A typical name would be
BIMODAL.SAS (the synthetic scattering described earlier) which
describes bimodal SAS data.  If the data file resides in a directory
other than the default, then you will have to specify that as part of
the file name.  Example file names, including a "full path
description" follow for the more popular operating systems:

Digital Equipment Corporation VAX running VMS::

        DISK$MPD_USERS:[CULVERWEL.MAXE]BIMODAL.SAS

Apple Macintosh::

        Hard Drive:MaxSas Folder:test case:BIMODAL.SAS

MS-DOS computer (such as the IBM-PC)::

        C:\MAXSAS\TEST\BIMODAL.SAS
