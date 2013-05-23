	PROGRAM Roess
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*25 ProgVers, EditDate
	PARAMETER (ProgVers = '1.0 (PRJ)')
	PARAMETER (EditDate = '18 August 1992')
C	Analysis of small-angle scattering data using the technique of
C	entropy maximization.

C	Designed for BATCH execution on the VAX

C	Uses scatterer shape factor of Roess and Shull.
C	This is an ellipsoid of revolution with principal radii
C	R x R x vR  where v is the aspect ratio.

C	Credits:
C	G.J. Daniell, Dept. of Physics, Southampton University, UK
C	J.A. Potton, UKAEA Harwell Laboratory, UK
C	I.D. Culverwell, UKAEA Harwell Laboratory, UK
C	G.P. Clarke, UKAEA Harwell Laboratory, UK
C	A.J. Allen, UKAEA Harwell Laboratory, UK
C	P.R. Jemian, Northwestern University, USA

C	References:
C	1. J Skilling and RK Bryan; MON NOT R ASTR SOC
C		211 (1984) 111 - 124.
C	2. JA Potton, GJ Daniell, and BD Rainford; Proc. Workshop
C		Neutron Scattering Data Analysis, Rutherford
C		Appleton Laboratory, UK, 1986; ed. MW Johnson,
C		IOP Conference Series 81 (1986) 81 - 86, Institute
C		of Physics, Bristol, UK.
C	3. ID Culverwell and GP Clarke; Ibid. 87 - 96.
C	4. JA Potton, GK Daniell, & BD Rainford,
C		J APPL CRYST 21 (1988) 663 - 668.
C	5. JA Potton, GJ Daniell, & BD Rainford,
C		J APPL CRYST 21 (1988) 891 - 897.
C	6. Roess and Shull; J Appl Phys 18 (1947) 308-313.

C	This is a modification of the program MAXE described below:
C	This progam (MAXE) was written in BASIC by GJ Daniell and later
C	  translated into FORTRAN and adapted for SANS analysis.  It
C	  has been further by PR Jemian to allow portability between
C	  the Digital Equipment Corporation VAX and Apple Macintosh
C	  computers. The form factor will now model homogeneous 
C	  spheroids of random orientation and arbitrary aspect ratio.
C	The input data file format is three columns of "Q I dI" which
C	  are separated by spaces or tabs.  There is no header line
C	  in the input data file.

	PARAMETER (cm2m = 0.01)	! convert cm to m units, but why?
	PARAMETER (MaxPts = 300, MaxBin = 102)
	PARAMETER (isLin = 1, isLog = 2, ioUnit = 1)

C  point-by-point mapping between reciprocal and real space
	COMMON /space1/ grid
	DIMENSION grid(MaxBin,MaxPts)

C  terms used in entropy maximization
	COMMON /space5/ chisq, chtarg, chizer, fSum, blank
	COMMON /space2/ beta, c1, c2, s1, s2
	DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)

C  terms used only by subroutine MaxEnt, allocated here to make memory tidy
	COMMON /space3/ ox, z, cgrad, sgrad, xi, eta
	DIMENSION ox(MaxPts), z(MaxPts)
	DIMENSION cgrad(MaxBin), sgrad(MaxBin)
	DIMENSION xi(MaxBin,3), eta(MaxPts,3)

C  space for the plotting frame, allocated here to make memory tidy
C    note the limits: MaxCol <= 100, MaxRow <= 150 (really large screens!)
	PARAMETER (MaxCol = 75, MaxRow = 15)
	PARAMETER (MxC2 = MaxCol+2, MxR2 = MaxRow+2)
	COMMON /space4/ screen, nCol, nRow, nCol2, nRow2
	CHARACTER*1 screen(100, 150)

C  space for main segment arrays
	DIMENSION q(MaxPts), datum(MaxPts), sigma(MaxPts)
	DIMENSION r(MaxBin), f(MaxBin), base(MaxBin), aNbr(MaxBin)
	DIMENSION fit(MaxPts), BinWid(MaxPts)
	CHARACTER*40 InFile, OutFil
	LOGICAL Yes
	CHARACTER*1 YN, aTab

	DATA one, zero /1.0, 0.0/	! for compiler-independence!

C ... Define (initially) the default responses
	DATA Aspect	/1.0/		! particle aspect ratio
	DATA LinLog	/isLin/		! linear binning scale
	DATA n		/40/		! number of bins
	DATA Dmin, Dmax	/8.00, 400.0/	! particle diameters
	DATA IterMax	/20/		! maximum number of iterations to try
	DATA RhoSq	/1.0/		! scattering contrast, x10**28 1/m**4
	DATA fac, err	/1.0, 1.0/	! scalars for intensity and errors
	DATA qMin, qMax	/1.e-8, 100./	! range to accept
	DATA Bkg	/0.0/		! intensity to subtract
	DATA sLengt	/1.0e-10/	! rectangular slit-length, 1/A
	DATA iPlot	/10/		! interim plots every iPlot-th iteration
	DATA dWave	/0.25/		! frac. FWHM of wavelength dispersion

C  Next line for MPW/Language Systems compiler, Macintosh only
C  Comment this out for other compilers
C  This is the only compiler-dependent line in this source code!!!!!!
C	CALL OutWindowScroll (1000)  ! for 1-line advance screen

	Pi = 4. * ATAN(1.)
	aTab = CHAR (9)

C  screen dimension variables for plots, in COMMON /space4/
	nCol = MaxCol
	nRow = MaxRow
	nCol2 = MxC2
	nRow2 = MxR2

    1	WRITE (*,*)
	WRITE (*,*) 'Size distributions from SAS data using the',
     >              ' maximum entropy criterion'
	WRITE (*,*) '   version: ', ProgVers, ',   ', EditDate

	CALL GetInf (InFile, OutFil, Aspect, LinLog,
     >		n, Dmin, Dmax, IterMax, RhoSq, fac, err, qMin,
     >		qMax, Bkg, sLengt, iPlot, dWave)
	IF (InFile .EQ. ' ') STOP

C   Read in the SAS data from the file "InFile"
	WRITE (*,*)  ' Reading from file: ', InFile
        OPEN (UNIT = ioUnit, FILE = InFile, STATUS = 'old')
	DO 94  j = 1, MaxPts
          READ (ioUnit, *, END = 95) q(j), datum(j), sigma(j)
   94   CONTINUE
   95   npt=j-1		! ignore any lines without an explicit EOL mark
        CLOSE (UNIT = ioUnit, STATUS = 'keep') 
	WRITE (*,*) npt, ' points were read from the file'

C   Subtract background, convert to 1/m units and extract data needed
     	i = 0			! count the points that are kept
        DO 2  j = 1, npt
	  IF (q(j) .GE. Qmin .AND. q(j) .LE. Qmax) THEN		!got one
	    i = i + 1
	    q(i) = q(j)
	    datum(i) = fac * (datum(j)-Bkg) / cm2m
	    sigma(i) = fac * err * sigma(j) / cm2m
	  END IF
    2	CONTINUE
	npt = i
	WRITE (*,*) npt, ' points were selected from the data'

C  PRJ: 24 May 1989
C	BinWid:	actual radial width of the indexed bin number
C	Step:	radial increment factor (for geometric series)
C	rWid:	radial width (for algebraic series)
	IF (LinLog .EQ. isLog) THEN	! geometric series
	  Step = (Dmax/Dmin)**(1. / FLOAT(n-1)) - 1.
	  rWid = 0.
	ELSE				! algebraic series
	  Step = 0.
	  rWid = 0.5*(Dmax - Dmin) / FLOAT(n-1)
	END IF
	r(1) = 0.5 * Dmin
	BinWid(1) = r(1) * Step + rWid
	DO 48  i = 2, n
	  r(i) = r(i-1) + BinWid(i-1)
	  BinWid(i) = r(i) * Step + rWid
   48	CONTINUE

! PRJ: 6 Feb 1990
!	Use Roess&Shull spheroid form factor
!	Note that sphere is a special case of this where Aspect = 1
	WRITE (*,*) ' Preparation of the GRID function...'
	WaveMean = 1.0 			! unit wavelength by default
	WaveMin = WaveMean - dWave	! for a triangular distribution
	WaveMax = WaveMean + dWave
	CALL Spheroid(r, N, q, npt, grid, aspect,
     >		WaveMean, dWave, WaveMin, WaveMax)  !figure the whole grid
	Pi43 = Pi * 4. / 3.
	DO i = 1, n				!multiply grid by other terms
	  pVol = Pi43 * Aspect * r(i)**3	!particle volume
	  DO j = 1, npt
	    grid(i,j) = cm2m * rhosq * pVol * grid(i,j)
	  END DO
	END DO

C	Attempt to account for scattering from very large and very small 
C	  particles using the limiting forms of grid(i,j).
	smallest = cm2m * rhosq * Pi43 * Aspect * r(1)**3		!V
	v = Aspect			!limit case for oblate spheroids
	IF (v .LT. 1.0) v = v**2	!limit case for prolate spheroids
	biggest = cm2m * rhosq * 3.0*Pi * (1.0 + v) / r(n) / Aspect	!S/V
	DO 227  j = 1, npt
	  thisQ4 = q(j)**3 * SQRT(q(j)**2 + sLengt**2)  !slit-length?
	  grid(n+1,j) = smallest
	  grid(n+2,j) = biggest / thisQ4
  227	CONTINUE

C   Submit the problem
  228	basis = 1.0e-6 / RhoSq		! Originally was 1.0e-6
	CALL MaxEnt (n+2, npt, f, datum, sigma, basis, base, 
     >		max, IterMax, iPlot)

C	"Max" counts the number of iterations inside MAXENT.
C	If Max < IterMax, then the problem has been solved.
	IF (Max .GE. IterMax) THEN
	  WRITE (*,*) ' No convergence! # iter. = ', Max
     	  WRITE (*,*) ' File was: ', InFile
C	  GO TO 1
	  STOP ' >>>>>>>>> NO CONVERGENCE <<<<<<<<<<<'
	END IF

C   Otherwise, SUCCESS!... so calculate the volume distribution
C	from the model SAS data
	CALL opus (n+2, npt, f, fit)

C   Calculate the absolute number distribution from the volume distribution
	modeV = 1	! peak of volume distribution
	modeN  = 1	! peak of number distribution
	DO 1915  i = 1, n
     	  pVol = Pi43 * Aspect * (r(i) * 1.e-8)**3	! V, cm**3
	  aNbr(i) = f(i) / pVol		! number / cm**3, in this bin
	  IF (f(i) .GT. f(modeV)) modeV = i
	  IF (aNbr(i) .GT. aNbr(modeN)) modeN = i
 1915	CONTINUE
	CALL Stats(r,aNbr,n, SumN, DnMean, DnSDev)	!from # dist.
	DnMean = 2.0 * DnMean		! put into diametral terms
	DnSDev = 2.0 * DnSDev		! put into diametral terms
	CALL Stats(r,f,n, SumV, DvMean, DvSDev)		!from vol. dist.
	DvMean = 2.0 * DvMean		! put into diametral terms
	DvSDev = 2.0 * DvSDev		! put into diametral terms
	VfTot = f(n+1) + SumV + f(n+2)	! don't forget the extra Vf

C   Figure the configurational entropy & convert to differential dist.s
	Entropy = zero
	DO 1919  i = 1, n
	  frac = f(i) / SumV			! normalized probability
	  Entropy = Entropy - frac * LOG (frac)
	  f(i) = f(i) / BinWid(i)
	  aNbr(i) = aNbr(i) / BinWid(i)
 1919	CONTINUE

C  Plot the differential volume-fraction distribution.
     	WRITE (*,*)
	WRITE (*,*) ' Input file: ', InFile
	WRITE (*,*) ' Volume weighted size dist.: V(r)N(r) versus r'
	CALL Plot (n, r, f)

	DO 918 j = 1, npt	! estimate residual background & its stdev
	  z(j) = one / (sigma(j)**2)	! statistical weighting
	  ox(j) = fit(j) -  datum(j)	! intensity differences
  918	CONTINUE
	CALL Stats(ox,z,nPt, xtra, shift, shiftDev) !from # dist.

C  Scale the data back to 1/cm units and calculate Chi-squared
	ChiSq  = zero
	Chi2Bk  = zero
	DO 919 j = 1, npt
	  z(j)  = (datum(j) - fit(j)) / sigma(j) 
	  ChiSq   = ChiSq + z(j)**2   
	  Chi2Bk  = Chi2Bk + (z(j) + shift/ sigma(j))**2
	  datum(j) = cm2m * datum(j)
	  sigma(j) = cm2m * sigma(j)
	  fit(j) = cm2m * fit(j)
  919	CONTINUE
	shift = cm2m * shift / fac
	shiftDev = cm2m * shiftDev / fac

	WRITE (*,*) ' standardized residuals vs. point number'
	CALL ResPlt (npt, z)

C  Let the file output begin!

	OPEN (UNIT = ioUnit, FILE=OutFil, STATUS='new')
	WRITE (ioUnit,*) ' Results of maximum entropy analysis of SAS'
	WRITE (ioUnit,*) '    version ',ProgVers, ', edited:', EditDate
	WRITE (ioUnit,*)
	WRITE (ioUnit,*) ' input file: ',  aTab, InFile
	WRITE (ioUnit,*) ' output file: ', aTab, OutFil
	WRITE (ioUnit,*) ' --------------------------------------------'
	WRITE (ioUnit,*)
	WRITE (ioUnit,*) ' N(D) dD is number of particles/cm**3'
	WRITE (ioUnit,*) '    whose size is between D and D + dD'
	WRITE (ioUnit,*) ' f(D) = V(D) * N(D)'
	WRITE (ioUnit,*)
	WRITE (ioUnit, 35591) 'D, A', aTab, 'f(D), 1/A',
     >			aTab, 'N(D), 1/A/cm^3', aTab, 'dD, A'
C	WRITE (ioUnit, 35591) '----', aTab, '--------------',
C     >			aTab, '--------------', aTab, '-----'
35591	FORMAT (1X, A12, A1, 1X, A15, A1, 1X, A15, A1, 1X, A6)

	DO 1001 i = 1, n	!write in terms of "diametral" distribution
 1001	  WRITE (ioUnit,3559) 2.*r(i), aTab, 0.5*f(i),
     >			aTab, 0.5*aNbr(i), aTab, 2.0*BinWid(i)
 3559	FORMAT (1X,F12.2, A1,1X,1PE15.5, A1,1X,E15.5, A1,1X,0PF12.4)


	WRITE (ioUnit,'(///)')
	WRITE (ioUnit, 1011) 'Q 1/A', aTab, 'I 1/cm', aTab,
     >			'MaxEnt 1/cm', aTab, 'dI 1/cm', aTab, 'z'
C	WRITE (ioUnit, 1011) '-----', aTab, '------', aTab,
C     >			'-------', aTab, '-------', aTab, '----'
 1011	FORMAT (A12, 3(A1, A12), 1X, A1, A12, 1X, A1, A12)

	DO 101  j = 1, npt
  101	  WRITE (ioUnit,560) q(j), aTab, datum(j), aTab, fit(j),
     >		aTab, sigma(j), aTab, z(j)
  560	  FORMAT (1PE12.4, 3(A1, E12.4), 1X, A1, 0PF12.6, 1X, A1, F12.6)

	WRITE (ioUnit,3301) InFile
	WRITE (*,3301) InFile
 3301	FORMAT (//' Input data: ', A40)

	WRITE (ioUnit,3302) RhoSq
	WRITE (*,3302) RhoSq
 3302	FORMAT (' Contrast = ', F15.7,' x 10^28 m^-4.')

	WRITE (ioUnit,*) ' spheroid: D x D x D*', Aspect
	WRITE (*,*) ' spheroid: D x D x D*', Aspect

	WRITE (ioUnit,53303) fac
	WRITE (*,53303) fac
53303	FORMAT (' Data conversion factor to 1/cm = ', 1PE12.5)

	WRITE (ioUnit,63303) err
	WRITE (*,63303) err
63303	FORMAT (' Error scaling factor = ', 1PE12.5)

	WRITE (ioUnit,63304) dWave
	WRITE (*,63304) dWave
63304	FORMAT (' Wavelength distribution frac. FWHM = ', 
     >		F12.5)

	IF (LinLog .EQ. isLog) THEN
	  WRITE (ioUnit,13304) 'geometric'
	  WRITE (*,13304) 'geometric'
	ELSE
	  WRITE (ioUnit,13304) 'algebraic'
	  WRITE (*,13304) 'algebraic'
	END IF
13304	FORMAT (' Histogram bins are distributed in an increasing ',
     >		A9, ' series.')

	WRITE (ioUnit,3304) 'Minimum', Dmin
	WRITE (*,3304) 'Minimum', Dmin
	WRITE (ioUnit,3304) 'Maximum', Dmax
	WRITE (*,3304) 'Maximum', Dmax
 3304	FORMAT (1X, A7, ' particle dimension D = ',F12.2,' A.')

	WRITE (ioUnit,3306) n
	WRITE (*,3306) n
 3306	FORMAT (' Number of histogram bins = ',I4,'.')

	WRITE (ioUnit,3307) itermax
	WRITE (*,3307) itermax
 3307	FORMAT (' Maximum number of iterations allowed = ',I4,'.')

	WRITE (ioUnit,3314) max
	WRITE (*,3314) max
 3314	FORMAT (' Program left MaxEnt routine after ',
     *    I4,' iterations.')

	WRITE (ioUnit,3308) npt
	WRITE (*,3308) npt
 3308	FORMAT (' Target chi-squared (# data points) = ',I5,'.')

	WRITE (ioUnit,3309) ChiSq
	WRITE (*,3309) ChiSq
 3309	FORMAT (' Best value of chi-squared achieved = ',F12.6,'.')

	WRITE (ioUnit, 33091) 'the final', Entropy
	WRITE (*, 33091) 'the final', Entropy
	WRITE (ioUnit, 33091) 'a flat', LOG (FLOAT (n))
	WRITE (*, 33091) 'a flat', LOG (FLOAT (n))
33091	FORMAT (' Entropy of ', A9, ' distribution = ', F12.7,'.')

	WRITE (ioUnit,33101) SumN
	WRITE (*,33101) SumN
33101	FORMAT (' Total particles  = ', 1PE15.5,' per cubic cm.')

	WRITE (ioUnit,3310) SumV
	WRITE (*,3310) SumV
 3310	FORMAT (' Total volume fraction of all scatterers = ',
     *    F15.9,'.')

	WRITE (ioUnit,3311) 'smaller', Dmin, f(n+1)/VfTot
	WRITE (ioUnit,3311) 'larger',  Dmax, f(n+2)/VfTot
	WRITE (*,3311) 'smaller', Dmin, f(n+1)/VfTot
	WRITE (*,3311) 'larger',  Dmax, f(n+2)/VfTot
 3311	FORMAT (' Part of distribution ',A7,' than ', F12.2,
     *    ' A = ', 2PF14.8,'%.')

	WRITE (ioUnit,3312) 'Volume', 'mode D value', 2.0 * r(modeV)
	WRITE (*,3312) 'Volume', 'mode D value', 2.0 * r(modeV)
	WRITE (ioUnit,3312) 'Volume', 'mean D value', DvMean
	WRITE (*,3312) 'Volume', 'mean D value', DvMean
	WRITE (ioUnit,3312) 'Volume', 'std. deviation', DvSDev
	WRITE (*,3312) 'Volume', 'std. deviation', DvSDev
	WRITE (ioUnit,3312) 'Number', 'mode D value', 2.0 * r(modeN)
	WRITE (*,3312) 'Number', 'mode D value', 2.0 * r(modeN)
	WRITE (ioUnit,3312) 'Number', 'mean D value', DnMean
	WRITE (*,3312) 'Number', 'mean D value', DnMean
	WRITE (ioUnit,3312) 'Number', 'std. deviation', DnSDev
	WRITE (*,3312) 'Number', 'std. deviation', DnSDev
 3312	FORMAT (1X, A6, '-weighted ', A14, ' = ', F12.5, ' A.')

	WRITE (ioUnit,3313) 'Min', q(1)
	WRITE (*,3313) 'Min', q(1)
	WRITE (ioUnit,3313) 'Max', q(npt)
	WRITE (*,3313) 'Max', q(npt)
 3313	FORMAT (1X, A3,'imum Q-vector = ', 1PE15.7, ' 1/A.')

	WRITE (ioUnit,3315) 'User-specified', Bkg
	WRITE (*,3315) 'User-specified', Bkg
	WRITE (ioUnit,3315) 'Suggested', Bkg - shift
	WRITE (*,3315) 'Suggested', Bkg - shift
	WRITE (ioUnit,3315) 'StDev of shift in', shiftDev
	WRITE (*,3315) 'StDev of shift in', shiftDev
 3315	FORMAT (1X, A17, ' background = ', F18.9,' input data units')

	WRITE (ioUnit,*) ' New background should give ChiSq = ', Chi2Bk
	WRITE (*,*) ' New background should give ChiSq = ', Chi2Bk

	CLOSE (UNIT=ioUnit, STATUS='keep')

C  Adjust the background default setting
C  Shift the intensity data just in case the user wants a Stability Check
C  Remember: background shifts down, intensity shifts up
C  Don't forget to put the data back into 1/m units!
     	Bkg = Bkg - shift
	DO 4010  j = 1, npt
	  datum(j) = (datum(j) + shift) / cm2m
	  sigma(j) = sigma(j) / cm2m
 4010	CONTINUE

C  Can't do this next step in BATCH mode, just gotta go on.
C
C	IF (ABS ((Chi2Bk-ChiSq)/FLOAT (npt)) .LE. 0.05) THEN
C	  WRITE (*,*) ' The change in ChiSquared should be < 5%.'
C 4000	  WRITE (*,*) ' Run the Stability Check? (Y/<N>)'
C	  READ (*,'(A1)') YN
C	  IF (YN .EQ. 'y'  .OR.  YN .EQ. 'Y') GO TO 228
C	  IF (YN.NE.' ' .AND. YN.NE.'n' .AND. YN.NE.'N') GO TO 4000
C	END IF

	WRITE (*,3200) OutFil
 3200	FORMAT (/,' The program is finished.', /, 
     >    ' The output file is: ', A40)
C	GO TO 1

 3199	STOP 
	END


	SUBROUTINE GetInf (InFile, OutFil, Aspect, LinLog,
     >		nBin, Dmin, Dmax, IterMax, RhoSq, fac, err, qMin,
     >		qMax, Bkg, sLengt, iPlot, dWave)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*40 InFile, OutFil
	PARAMETER (Ro2Max = 1.e6, ItrLim = 2000, AbsMax = 1.e3)
	PARAMETER (DiaMin = 1., DiaMax = 1.e6, ErrMax = 1.e6)
	PARAMETER (MaxPts = 300, MaxBin = 102)
	PARAMETER (isLin = 1, isLog = 2)

    1	WRITE (*,*) ' Input file? <Quit>'
	  READ (*, 2) InFile
    2	  FORMAT (A40)
	  IF (InFile.EQ.' ') RETURN

    3	WRITE (*,*) ' Output file?'
	  READ (*, 2) OutFil
	  IF (OutFil .EQ. ' ') GO TO 3
	  IF (OutFil .EQ. InFile) GO TO 1

	suggest = qMin
   16	WRITE (*,'(A,F,A,$)') ' Minimum q-vector? [1/A] <',
     >		suggest, '> '
	READ (*, '(F10.0)') qMin
	IF (qMin .LT. 0) GO TO 16
	IF (qMin .EQ. 0) qMin = suggest

	suggest = qMax
   17	WRITE (*,'(A,F,A,$)') ' Maximum q-vector? [1/A] <',
     >		suggest, '> '
	READ (*, '(F10.0)') qMax
	IF (qMax .EQ. 0) qMax = suggest
	IF (qMax .LE. 0) GO TO 17
	IF (qMax .LE. qMin) GO TO 1

	suggest = RhoSq
   13	WRITE (*,'(A,1PE,A,$)')
     >		' Scattering contrast? [10^28 m^-4] <',
     >		suggest,'> '
	READ (*, '(F10.0)') RhoSq
	IF (RhoSq .EQ. 0) RhoSq = suggest
	IF (RhoSq .LT. 0 .OR. RhoSq .GT. Ro2Max) GO TO 13

	suggest = fac
   14	WRITE (*,'(A,1PE,A,$)')
     >		' Factor to convert data to 1/cm? <',
     >		suggest, '> '
	READ (*, '(F10.0)') fac
	IF (fac .EQ. 0) fac = suggest
	IF (fac .LE. 0 .OR. fac .GT. AbsMax) GO TO 14

	suggest = err
   15	WRITE (*,'(A,F,A,$)') ' Error scaling factor? <',
     >		suggest, '> '
	READ (*, '(F10.0)') err
	IF (err .EQ. 0) err = suggest
	IF (err .LE. 0 .OR. err .GT. ErrMax) GO TO 15

	suggest = Bkg
   18	WRITE (*,'(A,1PE,A,$)') ' Background? <', suggest, '> '
	READ (*, '(F10.0)') Bkg
	IF (Bkg .EQ. 0) Bkg = suggest

	suggest = Aspect
    6	WRITE (*,'(A,1PE,A,$)')
     >		' Spheroids: D x D x vD,  Aspect ratio (v)?  <',
     >			suggest, '> '
	  READ (*,'(F10.0)') Aspect
	  IF (Aspect .EQ. 0) Aspect = suggest
	  IF (Aspect .LT. 0) GO TO 6

C	Reserved for future use to include slit-smearing effects
C	suggest = sLengt
C   61	WRITE (*,'(A,F,A,$)') ' Slit-length [1/A]? <', suggest, '> '
C	  READ (*,'(F10.0)') sLengt
C	  IF (sLengt .EQ. 0) sLengt = suggest
C	  IF (sLengt .LT. 0) GO TO 61

	Last = LinLog
    7	WRITE (*,'(A,I,A,$)') ' Bin step scale? (1=Linear, 2=Log) <',
     >		Last, '> '
	READ (*, '(I4)') LinLog
	IF (LinLog .EQ. 0) LinLog = Last
	IF (LinLog .NE. isLin  .AND. LinLog .NE. isLog) GO TO 7

	Last = nBin
    8	WRITE (*,'(A,I,A,$)') ' Number of histogram bins? <',
     >		Last, '> '
	READ (*, '(I4)') nBin
	IF (nBin .EQ. 0) nBin = Last
	IF (nBin .LT. 2 .OR. nBin .GT. (MaxBin-2)) GO TO 8

	suggest = Dmax
    9	WRITE (*,'(A,F,A,$)') ' Maximum value of D? [A] <',
     >		suggest, '> '
	READ (*, '(F10.0)') Dmax
	IF (Dmax .EQ. 0) Dmax = suggest
	IF (Dmax .LT. nBin*DiaMin .OR. Dmax .GE. DiaMax) GO TO 9

	Suggest = Dmax / FLOAT (nBin)
   11	WRITE (*,'(A,F,A,$)') ' Minimum value of D? [A] <',
     >		suggest, '> '
	READ (*, '(F10.0)') Dmin
	IF (Dmin .EQ. 0) Dmin = suggest
	IF (Dmin .GE. DMax .OR. Dmin .LT. DiaMin) GO TO 1

	IF (IterMax .GT. ItrLim) IterMax = ItrLim
	Last = IterMax
   12	WRITE (*,'(A,I,A,$)') ' Maximum number of iterations? <',
     >		Last, '> '
	READ (*, '(I4)') IterMax
	IF (IterMax .EQ. 0) IterMax = Last
	IF (IterMax .LT. 0 .OR. IterMax .GT. ItrLim) GO TO 12

	Last = iPlot
   20	WRITE (*,'(A,I,A,$)') ' Iterations between interim plots? <',
     >		Last, '> '
	READ (*, '(I4)') iPlot
	IF (iPlot .EQ. 0) iPlot = Last

	Suggest = dWave
   22	WRITE (*,'(A,F,A,$)') ' frac. wavelength dist. FWHM? [A] <', 
     >		suggest, '> '
	READ (*, '(F10.0)') suggest
	IF (suggest .LT. 0.0  .OR.  suggest .GE. 1.0) GO TO 22
	IF (suggest .GT. 0.0) dWave = suggest

	RETURN
	END


	SUBROUTINE opus(n,npt,x,ox)	! solution-space -> data-space
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (MaxPts=300, MaxBin=102)
	COMMON /space1/ grid
	DIMENSION x(MaxBin), grid(MaxBin,MaxPts), ox(MaxPts)
	DO 3  j = 1, npt
	  sum = 0.
	  DO 4  i = 1, n
	   sum = sum + x(i) * grid(i,j)
    4	  CONTINUE
	  ox(j) = sum
    3	CONTINUE
	RETURN
	END


	SUBROUTINE tropus(n,npt,ox,x)	! data-space -> solution-space
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (MaxPts=300, MaxBin=102)
	COMMON /space1/ grid
	DIMENSION x(MaxBin), grid(MaxBin,MaxPts), ox(MaxPts)
	DO 5  i = 1, n
	  sum = 0.
	  DO 6  j = 1, npt
	    sum = sum + ox(j) * grid(i,j)
    6	  CONTINUE
	  x(i) = sum
    5	CONTINUE
	RETURN
	END


	SUBROUTINE MaxEnt (n, npt, f, datum, sigma,  flat, base,
     >		iter, IterMax, iPlot)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (MaxPts=300, MaxBin=102)
	DIMENSION f(MaxBin), datum(MaxPts), sigma(MaxPts)
	DIMENSION base(MaxBin)

	COMMON /space1/ grid
	DIMENSION grid(MaxBin,MaxPts)

	COMMON /space5/ chisq, chtarg, chizer, fSum, blank
	COMMON /space2/ beta, c1, c2, s1, s2
	DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)

	COMMON /space3/ ox, z, cgrad, sgrad, xi, eta
	DIMENSION ox(MaxPts), z(MaxPts)
	DIMENSION cgrad(MaxBin), sgrad(MaxBin)
	DIMENSION xi(MaxBin,3), eta(MaxPts,3)

     	DIMENSION Entropy(201), Convrg(201)
	PARAMETER (TstLim = 0.05)	! for convergence
	DATA one, zero /1.0, 0.0/	! compiler-independence!

	blank = flat
	exp1 = EXP(one)

	IF (blank .EQ. zero) THEN
	  DO 1004 i = 1, n
 1004	    blank = blank + base(i)
	  blank = blank / FLOAT(n)
	  WRITE (*,*) ' Average of BASE = ', blank
	ELSE
	  WRITE (*,*) ' Setting BASE constant at ', blank
	  DO 1003 i = 1, n
 1003	    base(i) = blank
	ENDIF

     	WRITE (*,*) ' MaxEnt routine beginning ...'

	chizer = FLOAT(npt)
	chtarg = chizer
	m = 3
	DO 8 i = 1, n
    8	  f(i) = base(i)	! initial distribution is featureless

	iter = 0
    6	iter = iter + 1		! The iteration loop begins here!
	CALL opus (n, npt, f, ox)	! calc. the model intensity from "f"
	chisq = zero
	DO 10 j = 1, npt
	  a = (ox(j) - datum(j)) / sigma(j)
	  chisq = chisq + a**2
   10	  ox(j) = 2. * a / sigma(j)
	CALL tropus(n,npt,ox,cgrad)	! cGradient = Grid * ox
	test = zero	! mismatch between entropy and ChiSquared gradients
	snorm = zero	! entropy term
	cnorm = zero	! ChiSqr term
	tnorm = zero	! norm for the gradient term TEST
	fSum = zero	! find the sum of the f-vector
	DO 12  i = 1, n
	  fSum = fSum + f(i)
	  sgrad(i) = -LOG(f(i)/base(i)) / (blank*exp1)
	  snorm = snorm + f(i) * sgrad(i)**2
	  cnorm = cnorm + f(i) * cgrad(i)**2
	  tnorm = tnorm + f(i) * sgrad(i) * cgrad(i)
   12	CONTINUE
	snorm = SQRT(snorm)
	cnorm = SQRT(cnorm)
	a = one
	b = one / cnorm
	IF (iter .GT. 1) THEN
	  test = SQRT(0.5*(one-tnorm/(snorm*cnorm)))
	  a = 0.5 / (snorm * test)
	  b = 0.5 * b / test
	ENDIF
	DO 13 i = 1, n
	  xi(i,1) = f(i) * cgrad(i) / cnorm
	  xi(i,2) = f(i) * (a * sgrad(i) - b * cgrad(i))
   13	CONTINUE
	CALL opus (n,npt,xi(1,1),eta(1,1))
	CALL opus (n,npt,xi(1,2),eta(1,2))
	DO 14 j = 1, npt
	  ox(j) = eta(j,2) / (sigma(j)**2)
   14	CONTINUE
	CALL tropus (n,npt,ox,xi(1,3))
	a = zero
	DO 15 i = 1, n
	  b = f(i) * xi(i,3)
	  a = a + b * xi(i,3)
	  xi(i,3) = b
   15	CONTINUE
	a = one / SQRT(a)
	DO 16 i = 1, n
	  xi(i,3) = a * xi(i,3)
   16	CONTINUE
	CALL opus (n,npt,xi(1,3),eta(1,3))
	DO 17 k = 1, m
	  s1(k) = zero
	  c1(k) = zero
	  DO 18 i = 1, n
	    s1(k) = s1(k) + xi(i,k) * sgrad(i)
	    c1(k) = c1(k) + xi(i,k) * cgrad(i)
   18	  CONTINUE
	  c1(k) = c1(k) / chisq
   17	CONTINUE
	DO 19 k = 1, m
	  DO 19 l = 1, k
	    s2(k,l) = zero
	    c2(k,l) = zero
	    DO 20 i = 1, n
	      s2(k,l) = s2(k,l) - xi(i,k) * xi(i,l) / f(i)
   20	    CONTINUE
	    DO 21 j = 1, npt
	      c2(k,l) = c2(k,l) + eta(j,k) * eta(j,l) / (sigma(j)**2)
   21	    CONTINUE
	    s2(k,l) = s2(k,l) / blank
	    c2(k,l) = 2. * c2(k,l) / chisq
   19	CONTINUE
	c2(1,2) = c2(2,1)
	c2(1,3) = c2(3,1)
	c2(2,3) = c2(3,2)
	s2(1,2) = s2(2,1)
	s2(1,3) = s2(3,1)
	s2(2,3) = s2(3,2)
	beta(1) = -0.5 * c1(1) / c2(1,1)
	beta(2) = zero
	beta(3) = zero
	IF (iter .GT. 1) CALL Move(3)

C  Modify the current distribution (f-vector)
	fSum = zero		! find the sum of the f-vector
	fChange = zero		! and how much did it change?
	DO 23 i = 1, n
	  df = beta(1)*xi(i,1)+beta(2)*xi(i,2)+beta(3)*xi(i,3)
	  IF (df .LT. -f(i)) df = 0.001 * base(i) - f(i)	! a patch
	  f(i) = f(i) + df		! adjust the f-vector
	  fSum = fSum + f(i)
	  fChange = fChange + df
   23	CONTINUE

	s = zero
	DO 24  i = 1, n
	  temp = f(i) / fSum		! fraction of f(i) in this bin
	  s = s - temp * LOG (temp)	! from Skilling and Bryan, eq. 1
   24	CONTINUE

     	CALL opus (n, nPt, f, z)	! model the data-space from f(*)
	ChiSq = zero			! get the new ChiSquared
     	DO 25  j = 1, nPt
     	  z(j) = (datum(j) - z(j)) / sigma(j)	! the residuals
	  ChiSq = ChiSq + z(j)**2	! report this ChiSq, not the one above
   25	CONTINUE

     	Entropy(iter) = s
     	Convrg(iter) = LOG (ChiSq)

  300	CONTINUE
	IF (MOD (iter, iPlot) .EQ. 0) THEN
	  WRITE (*,*)
	  WRITE (*,*) ' Residuals'
	  CALL ResPlt (npt, z)

	  WRITE (*,*)
	  WRITE (*,*) ' Distribution'
	  CALL BasPlt (n, f, base)
	END IF

	WRITE (*,*) ' #', iter, ' of ', itermax, ',  n  = ', npt
	WRITE (*,200) test, s
  200	FORMAT (' test = ', F9.5, ',  Entropy = ', F12.7)
	WRITE (*,201) 'SQRT((Chi^2)/n)', 'target', SQRT(chtarg/npt),
     >			'% off', SQRT(chisq/chtarg)-1.0
	WRITE (*,201) 'f-vector', 'sum', fSum,
     >			'% change', fChange/fSum
  201	FORMAT (A17, ':', A8, ' = ', F12.8, A10,' = ', 2PF10.4)

C  See if we have finished our task.
	IF (ABS(chisq/chizer-one) .LT. 0.01) THEN  ! hardest test first
	  IF (test .LT. TstLim) THEN		! same solution gradient?
C		We've solved it but now must check for a bizarre condition.
C		Calling routine says we failed if "iter = iterMax".
C		Let's increment (maybe) iterMax so this doesn't happen.
	    IF (iter .EQ. iterMax) iterMax = iterMax + 1
	    RETURN
	  END IF
	END IF
	IF (iter .LT. iterMax) GO TO 6

	RETURN
C  Can't ask this question in BATCH mode.  Gotta just return
C
C	C  Ask for more time to finish the job.
C		WRITE (*,*)
C		WRITE (*,*) ' Maximum iterations have been reached.'
C	 2001	WRITE (*,*) ' How many more iterations? <none>'
C		READ (*,'(I4)') more
C		IF (more .LT. 0) GO TO 2001
C		IF (more .EQ. 0) RETURN
C		iterMax = iterMax + more
C		GO TO 6
	END


	SUBROUTINE Move(m)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER ( MxLoop = 500 )	! for no solution
	PARAMETER ( Passes = 1.e-3 )	! convergence test
	COMMON /space5/ chisq, chtarg, chizer, fSum, blank
	COMMON /space2/ beta, c1, c2, s1, s2
	DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
	DATA one, zero /1.0, 0.0/	! compiler-independence!
	a1 = zero			! lower bracket  "a"
	a2 = one			! upper bracket of "a"
        cmin = ChiNow (a1, m)
	IF (cmin*chisq .GT. chizer) THEN
	  ctarg = 0.5*(one + cmin)
	ELSE
	  ctarg = chizer/chisq
	END IF
	f1 = cmin - ctarg
	f2 = ChiNow (a2,m) - ctarg
	DO 1  loop = 1, MxLoop
	  anew = 0.5 * (a1+a2)		! choose a new "a"
	  fx = ChiNow (anew,m) - ctarg
	  IF (f1*fx .GT. zero) THEN
	    a1 = anew
	    f1 = fx
	  END IF
	  IF (f2*fx .GT. zero) THEN
	    a2 = anew
	    f2 = fx
	  END IF
	  IF (abs(fx) .LT. Passes) GO TO 2
    1	CONTINUE

C  If the preceding loop finishes, then we do not seem to be converging.
C	Stop gracefully because not every computer uses control-C (etc.)
C	as an exit procedure.
	WRITE (*,*) ' Loop counter = ', MxLoop
	PAUSE ' No convergence in alpha chop (MOVE).  Press return ...'
	STOP ' Program cannot continue.'

    2	w = Dist (m)
	IF (w .LE. 0.1*fSum/blank) GO TO 1042
	DO 1044 k=1,m
	  beta(k) = beta(k) * SQRT(0.1 * fSum/(blank * w))
 1044	CONTINUE
 1042	chtarg = ctarg * chisq
	RETURN
	END


	DOUBLE PRECISION FUNCTION Dist (m)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	COMMON /space5/ chisq, chtarg, chizer, fSum, blank
	COMMON /space2/ beta, c1, c2, s1, s2
	DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
	DATA one, zero /1.0, 0.0/	! compiler-independence!
	w = zero
	DO 26  k = 1, m
	  z = zero
	  DO 27  l = 1, m
	    z = z - s2(k,l) * beta(l)
   27	  CONTINUE
	  w = w + beta(k) * z
   26	CONTINUE
	Dist = w
	RETURN
	END


	DOUBLE PRECISION FUNCTION ChiNow(ax,m)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	COMMON /space5/ chisq, chtarg, chizer, fSum, blank
	COMMON /space2/ beta, c1, c2, s1, s2
	DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
	DIMENSION a(3,3), b(3)
	DATA one, zero /1.0, 0.0/	! compiler-independence!
	bx = one - ax
	DO 28  k = 1, m
	  DO 29  l = 1, m
	    a(k,l) = bx * c2(k,l)  -  ax * s2(k,l)
   29	  CONTINUE
	  b(k) = -(bx * c1(k)  -  ax * s1(k))
   28	CONTINUE
	CALL ChoSol(a,b,m,beta)
	w = zero
	DO 31  k = 1, m
	  z = zero
	  DO 32  l = 1, m
	    z = z + c2(k,l) * beta(l)
   32	  CONTINUE
	  w = w + beta(k) * (c1(k) + 0.5 * z)
   31	CONTINUE
	ChiNow = one +  w
	RETURN
	END


	SUBROUTINE ChoSol(a, b, n, beta)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION fl(3,3), a(3,3), bl(3), b(3), beta(3)
	DATA one, zero /1.0, 0.0/	! compiler-independence!
	IF (a(1,1) .LE. zero) THEN
	  WRITE (*,*) ' Fatal error in CHOSOL: a(1,1) = ', a(1,1)
	  PAUSE ' Press <RETURN> to end program ...'
	  STOP ' Program cannot continue.'
	END IF
	fl(1,1) = SQRT(a(1,1))
	DO 35  i = 2, n
	  fl(i,1) = a(i,1) / fl(1,1)
	  DO 35  j = 2, i
	    z = zero
	    DO 36  k = 1, j-1
	      z = z + fl(i,k) * fl(j,k)
   36	    CONTINUE
	    z = a(i,j) - z
	    IF (j .EQ. i) THEN
	      fl(i,j) = SQRT(z)
	    ELSE
	      fl(i,j) = z / fl(j,j)
	    END IF
35      CONTINUE
	bl(1) = b(1) / fl(1,1)
	DO 37  i=2, n
	  z = zero
	  DO 38  k = 1, i-1
	    z = z + fl(i,k) * bl(k)
   38	  CONTINUE
	  bl(i) = (b(i) - z) / fl(i,i)
   37	CONTINUE
	beta(n) = bl(n) / fl(n,n)
	DO 39  i1 = 1, n-1
	  i = n - i1
	  z = zero
	  DO 40  k = i+1, n
	    z = z + fl(k,i) * beta(k)
   40	  CONTINUE
	  beta(i) = (bl(i) - z) / fl(i,i)
   39	CONTINUE
	RETURN
	END


	SUBROUTINE BlankScreen (nC, nR)
C	Blank the screen array and mark the borders.
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*1 Blank, hBordr, vBordr
	PARAMETER (Blank = ' ', hBordr = '-', vBordr = '|')
	COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
	CHARACTER*1 screen(100, 150)
	DO 1  j = 1, nC
	  DO 1  i = 1, nR
	    screen(i,j) = Blank
    1	CONTINUE
	DO 2  i = 2, nC - 1
	  screen(nR,i) = hBordr
    2	  screen(1,i) = hBordr
	DO 3  i = 2, nR - 1
	  screen(i,nC) = vBordr
    3	  screen(i,1) = vBordr
	RETURN
	END

	SUBROUTINE ResPlt (n, x)
C	Draw a plot of the standardized residuals on the screen.
C	Mark the rows of + and - one standard deviation.
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION x(1)
	CHARACTER*1 Symbol, resSym
	PARAMETER (Symbol = 'O', resSym = '=')
	COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
	CHARACTER*1 screen(100, 150)

     	IF (n .LT. 2) RETURN	! not enough data

C  Find out how many points to pack per column and how many columns
	nPack = 1 + INT(FLOAT (n) / MaxCol - 1./n)
	nCol = INT((n - 1./n)/nPack + 1)

	CALL BlankScreen (nCol + 2, MxR2)	! prep the "screen"

C  get the data limits
        xMax = 1.
	xMin = -1.
        DO 4 i = 1, n
	  IF (x(i) .GT. xMax) xMax = x(i)
	  IF (x(i) .LT. xMin) xMin = x(i)
    4	CONTINUE
	RowDel = (MaxRow - 1) / (xMax - xMin)

C  show the standard deviation bars
	mPlus = 1 + INT((1 - xMin)*RowDel + 1)
	mMinus = 1 + INT((-1 - xMin)*RowDel + 1)
	DO 5  i = 2, nCol + 1
	  screen(mMinus,i) = resSym
    5	  screen(mPlus,i) = resSym

C  draw the data (overdrawing the residuals bars if necessary)
	DO 6  i = 1, n
	  mCol = 1 + INT((i - 1./n)/nPack + 1)		! addressing function
	  mRow = 1 + INT((x(i) - xMin)*RowDel + 1)	! +1 for the plot frame
	  screen(mRow, mCol) = Symbol
    6	CONTINUE

C  convey the "screen" to the default output
	WRITE (*,*) nPack, ' point(s) per column'
	WRITE (*,*) 1./RowDel, ' standard deviations per row'
	DO 7  i = MxR2, 1, -1
    7	  WRITE (*,*) (screen(i,j), j = 1, nCol + 2)
	RETURN
	END


	SUBROUTINE BasPlt (n, x, basis)
C	Draw a plot of some data and indicate a basis line on the
C	the plot.  That is, indicate that line below which the data is
C	not meaningful.  The basis here is taken to be a constant.
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION x(1)
	CHARACTER*1 Symbol, BasSym
	PARAMETER (Symbol = 'O', BasSym = '=')
	COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
	CHARACTER*1 screen(100, 150)

     	IF (n .LT. 2) RETURN	! not enough data

C  Find out how many points to pack per column and how many columns
	nPack = 1 + INT(FLOAT (n) / MaxCol - 1./n)
	nCol = INT((n - 1./n)/nPack + 1)

	CALL BlankScreen (nCol + 2, MxR2)	! prep the "screen"

C  get the data limits
        xMax = basis
	xMin = basis
        DO 4 i = 1, n
	  IF (x(i) .GT. xMax) xMax = x(i)
	  IF (x(i) .LT. xMin) xMin = x(i)
    4	CONTINUE
	RowDel = (MaxRow - 1) / (xMax - xMin)

C  show the basis line
	mPlus = 1 + INT((basis - xMin)*RowDel + 1)
	DO 5  i = 2, nCol + 1
    5	  screen(mPlus,i) = basSym

C  draw the data (overdrawing the basis bars if necessary)
	DO 6  i = 1, n
	  mCol = 1 + INT((i - 1./n)/nPack + 1)		! addressing function
	  mRow = 1 + INT((x(i) - xMin)*RowDel + 1)	! +1 for the plot frame
	  screen(mRow, mCol) = Symbol
    6	CONTINUE

C  convey the "screen" to the default output
	WRITE (*,*) nPack, ' point(s) per column'
	WRITE (*,*) 1./RowDel, ' units per row'
	DO 7  i = MxR2, 1, -1
    7	  WRITE (*,*) (screen(i,j), j = 1, nCol + 2)
	RETURN
	END


        SUBROUTINE Plot (n,x,y)
C	Make a scatter plot on the default display device (UNIT=*).
C	MaxRow and MaxCol correspond to the display dimensions.
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION x(1), y(1)
	CHARACTER*1 Symbol
	PARAMETER (Symbol = 'O')
	COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
	CHARACTER*1 screen(100, 150)

     	IF (n .LT. 2) RETURN			! not enough data
	CALL BlankScreen (MxC2, MxR2)		! prep the "screen"

C  get the data limits
	xMin = x(1)
        xMax = x(1)
	yMin = y(1)
        yMax = y(1)
        DO 4 i = 2, n
	  IF (x(i).GT.xMax) xMax=x(i)
	  IF (x(i).LT.xMin) xMin=x(i)
	  IF (y(i).GT.yMax) yMax=y(i)
	  IF (y(i).LT.yMin) yMin=y(i)
    4	CONTINUE
	ColDel = (MaxCol - 1) / (xMax - xMin)
	RowDel = (MaxRow - 1) / (yMax - yMin)

C  data scaling functions are offset by +1 for plot frame
	DO 5  i = 1, n
	  mCol = 1 + INT((x(i) - xMin)*ColDel + 1)
	  mRow = 1 + INT((y(i) - yMin)*RowDel + 1)
    5	  screen(mRow, mCol) = Symbol

C  convey the "screen" to the default output
	WRITE (*,*) 1./ColDel, ' units per column'
	WRITE (*,*) 1./RowDel, ' units per row'
	DO 6  i = MaxRow + 2, 1, -1
    6	  WRITE (*,*) (screen(i,j), j = 1, MaxCol + 2)
        RETURN
        END


	SUBROUTINE Spheroid(r, nR, h, nH, G, v, 
     >		WaveMean, dWave, WaveMin, WaveMax)
C	Particle form factor for spheroids of R x R x vR, R = radius
C	  ref: Roess and Shull; J Appl Phys 18 (1947) 308-313.
C	  coded: PRJ on 7.2.90
C	  modified: PRJ on 4.May.92	included triangular wavelength dispersion
C	The method of eq.6 (second part) is implemented by direct evaluation
C	  using eq. (4) & (5) and numerical averaging.
C	Lookup tables are used to speed up the calculations.
C
C	h: array of scattering vectors numbered 1 to nH
C	r: array of radii numbered 1 to nR
C	G: grid of (form factors)**2 numbered 1 to nH & 1 to nR
C	v: aspect ratio
C	WaveMean: default (average) wavelength
C	dWave:    FWHM of wavelength
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION h(1), r(1)
	PARAMETER (MaxPts=300, MaxBin=102)
	DIMENSION G(MaxBin,MaxPts)
	PARAMETER (eta = 0.01)		!how close to a sphere?
	PARAMETER (nT = 500)		!# of table values
	DIMENSION t(nT), u(nT), hr(nT)	!table ordinates
	DIMENSION S(nT), SS(nT)		!table abcissae
	PARAMETER (nI = 100)		!# of integrands
	DIMENSION x(nI), y(nI), p(nI)	!integrands
	DOUBLE PRECISION FUNCTION Area
	Step(a,b,N) = EXP( LOG(b/a)/( N-1.0 ) ) + 1.0e-5	!step-size

C
C  Calculate the table of "spheres" form factors: S(u)
C
	IF (v .LT. 1.0) THEN
	  vMin = v
	  vMax = 1.0
	ELSE
	  vMin = 1.0
	  vMax = v
	END IF
	uMin =  h(1) *  r(1) * (WaveMean/WaveMax) * vMin
	uMax = h(nH) * r(nR) * (WaveMean/WaveMin) * vMax
	du = Step(uMin, uMax, nT)	! the step size
	uNow = uMin
	DO  i = 1, nT
	  u(i) = uNow			! keep for the lookup table
	  S(i) = SphereForm (uNow)	! exact calculation at u = uNow
	  uNow = uNow * du		! for the next u
	END DO
C	CALL WriteFile (u, S, nT, 'S_s.dat')


C
C  Integrate to get the ellipsoidal form factor
C
	IF (ABS(v-1) .GT. eta) THEN	! integrate only if not spherical
	  tMin = h( 1) * r( 1) * (WaveMean/WaveMax)
	  tMax = h(nH) * r(nR) * (WaveMean/WaveMin)
	  dt = Step(tMin, tMax, nT)			! the step size
	  DO j = 1, nI
	    x(j) = (j - 1.0) / (nI - 1.0)		! ordinate of integration
	  END DO
	  tNow = tMin					! first point
	  DO  i = 1, nT
	    t(i) = tNow					! keep for the lookup table
	    DO j = 1, nI
	      uNow = tNow * SQRT (1 + x(j)**2 * (v**2 - 1))
	      CALL LookUp (uNow, u, nT, index)		!get the index
	      y(j) = ThisForm(uNow, u,S,nT, index)	!get the value
	    END DO
	    SS(i) = Area(x, y, nI)			! integration
	    tNow = tNow * dt				! for the next t
	  END DO
	  DO i = 1, nT
	    u(i) = t(i)
	    S(i) = SS(i)
	  END DO
C	  CALL WriteFile (u, S, nT, 'S_e.dat')
	END IF


C
C  Integrate to include the wavelength distribution
C
	IF (dWave .GE. 1.0E-6) THEN	! integrate only if finite dispersion
C		define the wavelength distribution
	  dw = (WaveMax - WaveMin) / (nI - 1.0)
	  DO j = 1, nI
	    x(j) = WaveMin + (j - 1.0) * dw		! wavelengths for integration
	    p(j) = (1.0 - ABS((x(j)-WaveMean)/dWave))/dWave	! triangular dist.
C	    p(j) = Gaussian (x(j), WaveMean, dWave)	! dist. function
	  END DO
	  WaveNorm = Area(x, p, nI)		! to make unit area
C	WRITE (*,*) ' WaveMean = ', WaveMean
C	WRITE (*,*) ' dWave = ', dWave
C	WRITE (*,*) ' WaveNorm = ', WaveNorm
C	WRITE (*,*) ' Wave step = ', dw
C	CALL WriteFile (x, p, nI, 'P_lambda.dat')
	  tMin = h( 1) * r( 1)
	  tMax = h(nH) * r(nR)
	  dt = Step(tMin, tMax, nT)			! the step size
	  tNow = tMin					! first point
	  DO  i = 1, nT
	    t(i) = tNow					! keep for the lookup table
	    DO j = 1, nI
	      uNow = tNow * WaveMean / x(j)		!indexed into ellipsoids table
	      CALL LookUp (uNow, u,nT, index)		!get the index
	      y(j) = p(j) * ThisForm(uNow, u,S,nT, index)	!get the value
	    END DO
C	IF (i .EQ. nT-10) CALL WriteFile (x, y, nI, 'integ.dat')
	    SS(i) = Area(x, y, nI) / WaveNorm 		! integration
	    tNow = tNow * dt				! for the next t
	  END DO
	  DO i = 1, nT
	    u(i) = t(i)
	    S(i) = SS(i)
	  END DO
C	  CALL WriteFile (u, S, nT, 'S_we.dat')
	END IF


C
C  Calculate G( h(i), r(j) ) from the table
C
	DO  i = 1, nR		!at each scattering vector
	  DO  j = 1, nH		!at each radius
	    hrNow = r(i) * h(j)				!indexed into table
	    CALL LookUp (hrNow, u,nT, index)		!get the index
	    G(i,j) = ThisForm(hrNow, u,S,nT, index)	!get the value
	  END DO
	END DO
	RETURN
	END

	DOUBLE PRECISION FUNCTION SphereForm(u)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C	Standard form factor for spheres.  Square it for ease of use.
	SphereForm = (3.0 * (SIN(u) - u*COS(u)) / (u**3))**2
	RETURN
	END

	SUBROUTINE LookUp(u, x, n, j)
C	Find x[j-1] < u <= x[j] where 1 <= j <= n.  If u <= x[1], j = 1.
C	  If u > x[n], j = n + 1.  Returns the value of j.
C	u: value to search
C	x: array of table values, numbered 1 to n
C	j: returned index (also used as an input, possibly)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION x(1)
	IF (u .LE. x(1)) THEN	!below range?
	  j = 1
	  RETURN
	END IF
	IF (u .GT. x(n)) THEN	!above range?
	  j = n + 1
	  RETURN
	END IF
	IF (j .LT. 1  .OR.  j .GT. n) j = 1	! get j into the range
	IF (j .GT. 1) THEN	!Is it already bracketed?
	  IF (x(j-1) .LT. u  .AND.  u .LE. x(j)) RETURN
	END IF
	jlo = 1			!Nope, execute a binary search.
	jhi = n
10	j = (jhi + jlo + 1) / 2		!the "+1" forces rounding up
	IF (u .LE. x(j)) THEN
	  jhi = j			!it's in the in upper section
	ELSE
	  jlo = j			!it's in the in lower section
	END IF
	IF (jhi - jlo .GT. 1) GO TO 10	!not yet bracketed, squeeze limits
	j = jhi				!take the higher bracket
	RETURN
	END

	DOUBLE PRECISION FUNCTION ThisForm(u, x,S,n, j)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION x(1), S(1)
	IF (j .GT. n) THEN
	  sNow = SphereForm (u)		! no table data so get exactly
	ELSE				! get by 2-pt. linear interpolation
	  IF (j .EQ. 1) THEN		! this is a limiting case
	    x1 = 0.0
	    S1 = 1.0
	  ELSE				! the general case
	    x1 = x(j-1)
	    S1 = S(j-1)
	  END IF
	  x2 = x(j)			! always have this point
	  S2 = S(j)
	  sNow = S1 + (S2 - S1)*(u - x1)/(x2 - x1)
	END IF
	ThisForm = sNow
	RETURN
	END

	DOUBLE PRECISION FUNCTION Area(x, y, N)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION x(1), y(1)
	sum = 0.0
	DO i = 2, N
	  sum = sum + (x(i) - x(i-1)) * (y(i) + y(i-1))
	END DO
	Area = 0.5 * sum
	RETURN
	END

	DOUBLE PRECISION FUNCTION Gaussian(x, c, s)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	PARAMETER ( Root2Pi = 2.506628274631000502 )
	g = (x - c) / s
	Gaussian = (1 / s / Root2Pi) * EXP(-g*g/2)
	RETURN
	END

	SUBROUTINE Stats(x,w,n, total, average, stdev)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C	Determine statistics from weighted data
	DIMENSION x(1), w(1)
	total = 0.0
	average = 0.0
	stdev = 0.0
	DO 10  i = 1, n
	  total = total + w(i)
	  average = average + w(i) * x(i)
	  stdev = stdev + w(i) * x(i)**2
   10	CONTINUE
	average = average / total
	stdev = SQRT (stdev/total - average**2)
	RETURN
	END

	SUBROUTINE WriteFile (x, y, N, f)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION x(1), y(1)
	REAL*4 xx, yy
	CHARACTER*(*) f
	CHARACTER*1 aTab
	aTab = CHAR(9)
	OPEN (UNIT=1, FILE=f, STATUS='NEW')
	WRITE (*,*) 'file: ', f, ', points = ', N
	WRITE (1, *) 'x', aTab, 'y'
	DO i = 1, N
	  xx = x(i)
	  yy = y(i)
	  WRITE (1, *) xx, aTab, yy
	END DO
	CLOSE (UNIT = 22)
	RETURN
	END

