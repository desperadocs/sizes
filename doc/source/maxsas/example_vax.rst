.. $Id$

EXAMPLE OF THE PROGRAM EXECUTION ON A DEC Vax
======================================================

The following is an excerpt from the execution of MaxSas as it is
analyzing the supplied test distribution BIMODAL.SAS, a bimodal
distribution with two Gaussian peaks: one at 80 Angstroms with a sigma
of 20 A and the other at 200 A with a sigma of 60 A.  The Gaussian at
200 A is half the height of the one at 80 A.  For a further discussion
of this distribution, see [Culverwell, 1986].  The exact volume
fraction of the original distribution was not specified.


In the section to follow, all user responses will be all in upper case
where appropriate and will be followed by a {CR}, signifying that the
user has pressed the return key.  {CR} by itself signifies that the
user has accepted the default answer to the question, as shown in
<default>.  Almost all questions have a default response, the only
exception being the output file names which may take any value EXCEPT
blank or the same as the input file name.  Initially, all the defaults
are preset to the proper answers for the supplied test distribution
BIMODAL.SAS.  If you type in a new value, that will value will become
the default the next time that the question is asked.

The example to follow will be interrupted from time to time for
explanations.  These will be isolated between rows of "====" signs.
For further explanations of each question which is asked, see the
appropriate section appearing elsewhere in this document.

Now here's the excerpt from the beginning of the run.  Remember that
all user responses are on a single line and are terminated with a
{CR}.

::

	$ RUN MaxSas{CR}
	
	Size distributions from SAS data using the maximum entropy criterion
	   version: 3.1 (PRJ)                ,   7 February 1990          
	 Input file? <Quit>
	Bimodal.Sas{CR}
	 Output file?
	Bimodal.Out{CR}
	 Minimum q-vector? [1/A] <  1.000000000000000E-008>
	{CR}
	 Maximum q-vector? [1/A] <   100.000000000000     >
	{CR}
	 Scattering contrast? [10^28 m^-4] <   1.00000000000000     >
	{CR}
	 Factor to convert data to 1/cm? <   1.00000000000000     >
	{CR}
	 Error scaling factor? <   1.00000000000000     >
	{CR}
	 Background? <  0.000000000000000E+000>
	{CR}
	 Spheroids: D x D x vD,  Aspect ratio (v)?  <   1.00000000000000     >
	{CR}
	 Bin step scale? (1=Linear, 2=Log) <           1>
	{CR}
	 Number of histogram bins? <          40>
	{CR}
	 Maximum value of D? [A] <   400.000000000000     >
	{CR}
	 Minimum value of D? [A] <   10.0000000000000     >
	{CR}
	 Maximum number of iterations? <          20>
	{CR}
	 Reading from file: Bimodal.Sas                             
	         38 points were read from the file
	         38 points were selected from the data
	 Preparation of the GRID function...
	 Setting BASE constant at   1.000000000000000E-012
	 MaxEnt routine beginning ...


To recap what has happened so far, the program was started and the
input data file BIMODAL.SAS was specified for analysis.  All the
default answers appeared to be acceptable so the program read 38
points from the file, of which all 38 points were retained for the
analysis.  The program informs you that it has set an internal array
called BASE to 1.0E-12.  BASE is the initial guess for the
distribution.  If there is a value in the output distribution that is
comparable with this number, then that particular histogram bin has no
significant information.  Consider then that BASE is the "featureless"
distribution and rest assured that it is quite flat.

At this point, the program has gotten all the adjustable parameters
that it needs and is now proceeding to attempt to solve the problem.

The fitting routine is an iterative one. If there are N histogram bins
and P data points then the computation time for each iteration will be
on the order of (very approximate) Order(2N+P). (This is why both N
and P are bounded.)

All sorts of information will begin to scroll on the screen at an
alarming rate as the Maximimum Entropy routine, or MaxEnt, (see
[Skilling, 1984] for details on what's happening) begins to extract
statistically significant information from the SAS data.  With each
iteration, several different types of screen plots are drawn, each
describing the data extracted so far. The different plots are titled:

* LOG (ChiSq) vs. iteration number
* Entropy vs. iteration number
* Residuals
* Distribution

The first two plots will not appear until the third iteration.  They
are intended to keep you informed about the progress of the MaxEnt
routine.  The residuals (difference between MaxEnt fit and input data
normalized to the input errors) are plotted as a function of data
point subscript number, not Q-vector. The distribution is weighted by
the bin width, dD(i) = D(i+1) - D(i), and is also plotted as a
function of diametral bin subscript number. All plots will be scaled
to fit within the screen boundaries.

The last information to appear in each iteration will be a report,
such as:

::

	 #           2 of           20,  n  =           38
	test =   0.19101,  Entropy =    3.4668138
	 SQRT((Chi^2)/n):  target =  12.14204774     % off =    26.2698
	        f-vector:     sum =   0.00452313  % change =     0.4663

The report is saying that for the second iteration out of 20, where
there are 38 intensity data points, the difference between the entropy
and ChiSquared gradients is 0.19101, the entropy of the distribution
just plotted is 3.4668138 (whose units are exact), the target for the
SQRT((ChiSquared)/n), where ChiSquared is derived from intensity
calculated from the distribution just plotted, of 12.14204774 was
missed by 26.2698%, and the total volume fraction of scatterers in the
distribution just plotted is 0.452313 %, which did not change much
from the previous iteration.  The iterations will continue.

Here is the screen output from the last iteration::

	 LOG (ChiSq) vs. iteration number
	          1 point(s) per column
	 0.405168291335942      units per row
	 ----------- 
	|           |
	|OO         |
	|  O        |
	|   O       |
	|           |
	|    O      |
	|           |
	|     O     |
	|      O    |
	|       O   |
	|           |
	|        O  |
	|         O |
	|           |
	|==========O|
	 ----------- 
	
	 Entropy vs. iteration number
	          1 point(s) per column
	 5.229508533621462E-002 units per row
	 ----------- 
	|===========|
	|           |
	|           |
	|O          |
	|           |
	|           |
	| O         |
	|           |
	|           |
	|           |
	|  O        |
	|           |
	|        OOO|
	|      OO   |
	|   OOO     |
	 ----------- 
	
	 Residuals
	          1 point(s) per column
	 0.315439585860872      standard deviations per row
	 -------------------------------------- 
	|                                    O |
	|                                      |
	|  OO                          O       |
	|                                      |
	|                                     O|
	|====O==O======O====O=====O============|
	|                       O   O      O   |
	|          OO              O           |
	|     O            O  O         O      |
	|        OO   O   O      O   OO        |
	|OO             OO   O O               |
	|======O=========================OO====|
	|                                      |
	|                                      |
	|            O                      O  |
	 -------------------------------------- 
	
	 Distribution
	          1 point(s) per column
	 6.975972282206744E-005 units per row
	 ------------------------------------------ 
	|      O                                   |
	|                                          |
	|                                          |
	|       O                                  |
	|                                          |
	|                                          |
	|                                          |
	|                                          |
	|        O                                 |
	|                  OO                      |
	|     O   OOOO    O  O  O                  |
	|             O  O    OO OO                |
	|              OO                          |
	|                          OO              |
	|OOOOO=======================OOOOOOOOOOOOOO|
	 ------------------------------------------ 
	 #          11 of           20,  n  =           38
	test =   0.01587,  Entropy =    3.1391598
	 SQRT((Chi^2)/n):  target =   1.00000000     % off =     0.0504
	        f-vector:     sum =   0.00809492  % change =     4.0445


The problem has been solved in 11 iterations of the MaxEnt routine.
The two criteria for solution are that TEST <= 0.05 (5%) and that the
SQRT((Chi^2)/n) target be met within 0.5%.  Observe that the volume
fraction has not changed very much from the previous iteration.

Here is the summary screen output of the analysis::

	 Input file: Bimodal.Sas                             
	 Volume weighted size dist.: V(r)N(r) versus r
	  2.63513513513514      units per column
	 1.395194456441349E-005 units per row
	 --------------------------------------------------------------------------- 
	|           O                                                               |
	|                                                                           |
	|                                                                           |
	|             O                                                             |
	|                                                                           |
	|                                                                           |
	|                                                                           |
	|                                                                           |
	|               O                                                           |
	|                                  O O                                      |
	|         O       OO O O         O    O     O                               |
	|                        O     O        O O   O O                           |
	|                          O O                                              |
	|                                                 O O                       |
	|OO O O O                                             O OO O O O O O O O O O|
	 --------------------------------------------------------------------------- 
	 standardized residuals vs. point number
	          1 point(s) per column
	 0.315439585860872      standard deviations per row
	 -------------------------------------- 
	|                                    O |
	|                                      |
	|  OO                          O       |
	|                                      |
	|                                     O|
	|====O==O======O====O=====O============|
	|                       O   O      O   |
	|          OO              O           |
	|     O            O  O         O      |
	|        OO   O   O      O   OO        |
	|OO             OO   O O               |
	|======O=========================OO====|
	|                                      |
	|                                      |
	|            O                      O  |
	 -------------------------------------- 
	
	Input data: Bimodal.Sas                             
	Contrast =       1.0000000 x 10^28 m^-4.
	 spheroid: D x D x D*   1.00000000000000     
	Data conversion factor to 1/cm =  1.00000E+00
	Error scaling factor =  1.00000E+00
	Histogram bins are distributed in an increasing algebraic series.
	Minimum particle dimension D =        10.00 A.
	Maximum particle dimension D =       400.00 A.
	Number of histogram bins =   40.
	Maximum number of iterations allowed =   20.
	Program left MaxEnt routine after   11 iterations.
	Target chi-squared (# data points) =    38.
	Best value of chi-squared achieved =    38.038279.
	Entropy of the final distribution =    3.1389744.
	Entropy of    a flat distribution =    3.6888795.
	Total particles  =     1.57329E+16 per cubic cm.
	Total volume fraction of all scatterers =     0.008094741.
	Part of distribution smaller than        10.00 A =     0.00000000%.
	Part of distribution  larger than       400.00 A =     0.00215387%.
	Volume-weighted   mode D value =     70.00000 A.
	Volume-weighted   mean D value =    153.25044 A.
	Volume-weighted std. deviation =     72.92106 A.
	Number-weighted   mode D value =     70.00000 A.
	Number-weighted   mean D value =     83.89447 A.
	Number-weighted std. deviation =     33.47500 A.
	Minimum Q-vector =   7.4935700E-03 1/A.
	Maximum Q-vector =   9.9937470E-02 1/A.
	   User-specified background =        0.000000000 input data units
	        Suggested background =        0.000064250 input data units
	StDev of shift in background =        0.000330552 input data units
	 New background should give ChiSq =    36.6534794792953   


By now, the tabular data of the size distribution and the intensity
fit have been written to the data file BIMODAL.OUT.  The summary
analysis of the distribution has also been written to the output file.

The program has found about 0.8% by volume of scatterers.  The
ChiSquared matches the number of intensity points (a requirement for
solution) and the background suggested is not far from the background
used (with respect to the sigma of the last intensity of 0.000123
1/cm).  Observe that the standard deviation of the suggested shift
in the background is much larger than the suggested shift.

It appears that we have a reasonable solution in hand.  Don't be too
hasty to believe it yet.  In order to test the stability of the
answer, analyze the BIMODAL.SAS data again, using all the same
parameters as before (just take the defaults) except use the
background (0.000064251 1/cm) suggested by MaxSas.

Therefore, you should respond affirmatively to the Stability Check
question. Note the default answer on the Stability Check question is
"N".

::

	  The change in ChiSquared should be < 5%.
	  Run the Stability Check? (Y/<N>)
	Y{CR}
	  Setting BASE constant at    1.000000000000000E-12
	  MaxEnt routine beginning ...

For this documentation, the results of the stability check will not be
shown.  There is not much change in the volume fraction after the
stability check (about 1% or so) for the supplied test distribution.

After the Stability Check, you should get this question again.
This time, respond like the following to quit the program.

::

	  The change in ChiSquared should be < 5%.
	  Run the Stability Check? (Y/<N>)
	{CR}
	
	 The program is finished.
	 The output file is: BIMODAL.MAX
	
	Size distributions from SAS data using the maximum entropy criterion
	   version: 3.1 (PRJ)                ,   7 February 1990          
	  Input file? <Quit>
	{CR}
	
	$

