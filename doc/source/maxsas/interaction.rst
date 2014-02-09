User Interaction with the program
======================================================


A FEW WORDS ABOUT NUMERICAL RESPONSES BY THE USER

If you respond to a numerical question with a "zero", the default
answer will be used.  That is the way this program works to give you
default answers.  If you want to set a parameter to be "zero", use an
infinitesimal value such as 1.0E-25.

All floating point responses should include a decimal point somewhere
in the mantissa of the response, otherwise the results are
unpredictable and very system dependent!



EXPLANATION OF QUESTIONS ASKED BY THE PROGRAM

Q: Input file? <Quit>

The input file contains the SAS intensity data as ordered triples of
Q-vector (in 1/A units), Intensity (in arbitrary units), and
statistical error of intensity (same as intensity units).  Note that
there are no initial header lines in the input file.  No more than the
first 300 data points (ordered triples) will be read from the input
file.

If you were to press {CR} without typing in a file name, the program
would quit (as indicated by the default).

If the input file does not exist, the program will happily proceed to
ask you all the remaining questions it has.  Then and only then will
it find out that the file you named does not exist.  This will
generate a program crash.

A suggestion for input file name extensions is ".SAS" but this only a
suggestion.  The input file name may be up to 80 characters long.


Q: Output file?

This is the only question which has no default answer.  You must
answer this question with something.  If your answer is the same as
the input file name, the program will start over asking you for the
input file name.  This may be used as an easy exit if you specified
the wrong name.

The program does not check to see if the named output file already
exists.  On some systems (Macintosh and MS-DOS), the old file will be
erased and a new file created.  On other systems (VAX), a file with
the same name but a new version number will be created.  Forewarned is
forearmed.

My suggestion for MaxSas solution file name extensions is ".MAX" or
".DIS".  The output file name may be up to 80 characters long.


Q: Minimum (Maximum) q-vector? [1/A]

Use Q-vector (actually Q-vector magnitude) in units of 1/Angstrom.
The user is allowed to exclude data points from the ends of the input
data. Only those data points satisfying qMin <= Q-vector <= qMax will
be analyzed.  The program is designed to only handle positive
Q-vectors.

Initially, qMin is set ridiculously low so that even the lowest data
point will be used.  Correspondingly, qMax is set high enough to
include all typical SAS Q-vectors.

The user should generally cut off the data when the signal-to-noise
ratio becomes poor.  Truncating earlier than this will lose
information about the smallest particles present in the sample.  Users
might note that it is not neccessary for all the intensity values to
be positive, although it is probably inadvisable to include more than
five negative ones.


Q: Scattering contrast? [10^28 m^-4]

The scattering contrast is the squared difference between the
scattering length density of particle and matrix.  If the contrast is
1.27E30 1/m**4, then enter the value 127.0.
By the way, 1.E28 1/m**4 = 1.E20 1/cm**4.

The user can either enter the true contrast here or reply {CR}, in
which case the final "volume fractions" obtained will have to be
divided by the contrast (in units of 1.E28 1/m**4) in order to obtain
genuine volume fractions.  The program is coded to accept scattering
contrast values no larger than one million units of 1.0E28 (1.0e34)
1/m**4.


Q: Factor to convert data to 1/cm?

If the intensity values in the input file were not in units of 1/cm,
enter the constant to convert them into such units.  If they were
already in 1/cm units, good for you, so just press {CR} to accept the
default.  The program is coded to accept conversion values no larger
than 1000.0.


Q: Error scaling factor?

Here is an opportunity for you to try analyzing your data with
different ratios of signal to noise.  If you think that the errors in
the input file were underspecified, you may multiply them by this
constant.  More on this later as this will have a major influence upon
the analysis.


Q: Background?

This program has left you the opportunity to subtract a constant
intensity value.  A good initial approximation will put you on the
road to a good analysis of the data.  Remarkably, the background may
take any value, positive or negative.  If you want to set the
background back to zero, use infinitesimal (such as 1.E-25) rather
instead.  More on background later.


Q: Spheroids: D x D x vD,  Aspect ratio (v)?

The scattering form factor currently implemented is that discussed by
[Roess, 1947].  A special case of this ellipsoid of revolution,
whose outside dimensions are D x D x vD, is the sphere whose form
factor is described in [Culverwell, 1986] and [Potton, 1988a].

It is possible to select any aspect ratio (within reason) using this
model and the program only checks to see that you have entered a
positive value.  Special care has been taken to ensure that the volume
fractions determined by this model are correct.

For a full explanation of the coding of this model (from eq. 4, 5, & 6
of [Roess, 1947]), see the source code listing.  Look for the routine
named "Spheroid."

Remember that the distributions that are output are in terms of the
dimension "D".  The volume of this type of spheroid is (4Pi/3) v r**3.


Q: Bin step scale? (1=Linear, 2=Log)

"Linear" binning means that the diametral bins will increase in size
according to an algebraic series (e.g. 1, 2, 3, ...).  The other
method currently available is "logarithmic" binning where the increase
is according to a geometric series (e.g. 1.0, 1.05, 1.1025, ...).  Use
whichever method gives you a sufficient number of points over all the
peaks in the distribution.  Be aware that the calculated volume
fractions and number densities for the first few bins on the "log
scale" are likely to be artificially high because of the small bin
width and small particle volume corresponding to that bin (both these
terms divide the quantity that MaxSas derives to give you the volume
fraction").

The width of each bin indexed by "i" is dD(i) = D(i+1) - D(i) so that
the number density of scatterers whose size is between D and D + dD is
truly N(D) dD.  The bin width appears in the output file as "dD."


Q: Number of histogram bins?

This is an integer between 2 and 100, limited by computer memory and
execution CPU time.  Use as few bins as you think you need to
adequately describe the distribution or as many bins as you want, up
to the maximum of 100.


Q: Maximum (Minimum) value of D? [A]

Use Angstrom units.  Because each intensity is a statistical
representation of ALL dimensions D in the sample, weighted by a
particular form factor (model function), the choice of maximum and
minimum D is left to the user. You may specify values that are beyond
the "peripheral vision" of your data to see if there is any
statistical support for such sizes in your data. Usually, one knows
something about the size distribution to be solved and a maximum
particle diameter can be estimated.  Ideally Dmax should be an
over-estimate; if too small a diameter range (Dmin to Dmax) is
specified, the program will likely fail.

The largest value for Dmax is something unreasonable for most SAS data
(1 million Angstroms).  If you try to exceed this limit, the program
will patiently ask you again for the maximum D value.  The smallest
Dmin value you may enter is 1.0 Angstrom.  The program will always
suggest Dmin = Dmax / (number of bins).

If Dmin >= Dmax, the program will start asking you questions all over.
You can use this as an easy way to correct a bad input prior to this
question, without having to stop and restart the program.


Q: Maximum number of iterations?

The number of iterations is best estimated by experience.  Skilling
and Bryan [Skilling, 1984] suggest that one should re-consider the
model if more than about 20 iterations are required for convergence
within the Maximum Entropy routine (MaxEnt).  The largest allowed
number of iterations is 200 but if you require this, your model is
probably not representing the data well.  The MaxEnt routine may not
require as many iterations as you specify. That just means the job was
easier than you "thought".

If, while the MaxEnt routine is iterating, you see that a few more
iterations will be required to achieve a satisfying solution than you
have specified here, all is not lost.  If the limit specified is
reached with no satisfying maximum entropy solution yet in hand, the
program will ask you if you want to iterate more.  You can then extend
the process.  For this reason, it is suggested that you specify a
lower value (rather than higher) so that you may check the program's
progress.  A low limit allows the MaxEnt routine to escape should the
fitting process fail to converge.  In such an event, one or more of
the input parameters should be adjusted to achieve a more harmonius
solution.

A good general suggestion for the number of iterations is the maximum
number that you are willing to see the MaxEnt routine perform and not
converge.  If the MaxEnt routine needs more iterations, it will ask
you for permission.


Q: The change in ChiSquared should be < 5%. ::

   Run the Stability Check? (Y/<N>)

The Stability Check will perform the same analysis on the data set
with all the same parameters except that the suggested background will
be used. If the answer is stable, then all the results should be the
same.  If the answer is unsteady, then things will look different in
some way.  The prompt for a stability check will not appear unless the
program calculates that the shift should produce less than a 5%
change in the ChiSquared.


Q: Maximum iterations have been reached.::

   How many more iterations? <none>

This question occurs inside the MaxEnt routine when the maximum number
of iterations that you specified have been reached.  If you want the
MaxEnt routine to keep trying, specify a positive integer, otherwise
take the default which will generate the following output::

        No convergence! # iter. = "IterMax"
        File was: "InFile"

The program will then start over at the first question.
