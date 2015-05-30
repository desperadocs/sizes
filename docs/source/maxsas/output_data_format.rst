OUTPUT FORMAT
==================

One output file is created by the program - a data file whose name is
user supplied (e.g. BIMODAL.MAX).

The first few lines of output of a typical output file are as follows::

	 Results of maximum entropy analysis of SAS
	    version 3.1 (PRJ)                , edited:7 February 1990          
	
	 input file:    Bimodal.Sas                             
	 output file:   Bimodal.Out                             
	 --------------------------------------------
	
	 N(D) dD is number of particles/cm**3
	    whose size is between D and D + dD
	
	        D, A      V(D)*N(D), 1/A          N(D), 1/A/cm^3          dD, A
	        ----      --------------          --------------          -----
	       10.00         8.81630E-15             1.68379E+07              10.0000
	       20.00         2.16369E-14             5.16543E+06              10.0000
	       30.00         6.36171E-14             4.49999E+06              10.0000

The four columns are separated by TABs and also spaces.

=============   ===========================================================
term            description
=============   ===========================================================
D               dimension of scatterer in Angstroms
V(D)*N(D)       volume-weighted differential number distribution
N(D)            differential number distribution
dD              width of the bin in Angstroms in terms of "D"
=============   ===========================================================

V(D)*N(D) is the distribution solved by the program.  N(D) is the
distribution most often reported by other analysis techniques (TEM,
mercury porosimetry, etc.).  N(D) is calculated from the second one by
dividing by the particle volume.  This table continues as::


	      380.00         8.11590E-07             2.82480E+10              10.0000
	      390.00         5.74524E-07             1.84976E+10              10.0000
	      400.00         3.71052E-07             1.10728E+10              10.0000
	
	
	
	
	      Q 1/A           I 1/cm         ^I 1/cm         dI 1/cm               z
	      -----           ------         -------         -------            ----
	 7.4936E-03       2.1200E+00      2.2330E+00      1.6700E-01       -0.676923
	 9.9975E-03       1.9000E+00      1.9596E+00      1.0200E-01       -0.584780
	 1.2501E-02       1.7840E+00      1.6605E+00      6.6900E-02        1.845438



These columns are also separated by a combination of spaces and tabs

=====   ===========================================================
term    description
=====   ===========================================================
Q       input Q-vector in 1/Angstrom units
I       input intensity in 1/cm units
^I      intensity calculated from the distribution above
dI      input error in 1/cm, scaled by the error scaling factor
z       standardized residual, z = (I - ^I)/dI
=====   ===========================================================

::

	 9.2452E-02       4.3770E-03      4.2886E-03      1.6700E-04        0.529397
	 9.5008E-02       3.4510E-03      3.7178E-03      1.5000E-04       -1.778752
	 9.7564E-02       3.5330E-03      3.2077E-03      1.3500E-04        2.409922
	 9.9937E-02       2.9560E-03      2.7793E-03      1.2300E-04        1.436465
	
	
	Input data: Bimodal.Sas                             
	Contrast =       1.0000000 x 10^28 m^-4.
	 spheroid: D x D x D*   1.00000000000000     
	Data conversion factor to 1/cm =  1.00000E+00

After all the intensities have been printed, the same summary appears
in the output file as appeared on the screen shown above.  Only the
first few lines are shown here.

The program creates two extra pseudo-particle sizes, one that is
"smaller" than the smallest bin in the distribution and one that is
"larger" than the largest bin in the distribution.  The scattering
from these pseudo-particles is approximated at low angles by the
Guinier relation and at high angles by the Porod relation.  The object
of the user is to define the dimension range large enough that neither
of these pseudo-particles develops any significant volume fraction.
The percentage of the distribution that was assigned to each of these
sizes is reported.  Provided that both are small compared to the total
volume fraction of scatterers (also given in the summary) then the
program has detected essentially all of the particles contributing to
the observed scattering.  If they are not, it may be worth re-running
the program with a different diameter range in order to detect this
extraneous volume fraction.
