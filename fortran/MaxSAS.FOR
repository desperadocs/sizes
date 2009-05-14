        PROGRAM MaxSAS
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        CHARACTER*25 ProgVers, EditDate
        PARAMETER (ProgVers = '3.0 (PRJ)')
        PARAMETER (EditDate = '27 November 1989')
C       Analysis of small-angle scattering data using the technique of
C       entropy maximization.

C       Credits:
C       G.J. Daniell, Dept. of Physics, Southampton University, UK
C       J.A. Potton, UKAEA Harwell Laboratory, UK
C       I.D. Culverwell, UKAEA Harwell Laboratory, UK
C       G.P. Clarke, UKAEA Harwell Laboratory, UK
C       A.J. Allen, UKAEA Harwell Laboratory, UK
C       P.R. Jemian, Northwestern University, USA

C       References:
C       1. J Skilling and RK Bryan; MON NOT R ASTR SOC
C               211 (1984) 111 - 124.
C       2. JA Potton, GJ Daniell, and BD Rainford; Proc. Workshop
C               Neutron Scattering Data Analysis, Rutherford
C               Appleton Laboratory, UK, 1986; ed. MW Johnson,
C               IOP Conference Series 81 (1986) 81 - 86, Institute
C               of Physics, Bristol, UK.
C       3. ID Culverwell and GP Clarke; Ibid. 87 - 96.
C       4. JA Potton, GK Daniell, & BD Rainford,
C               J APPL CRYST 21 (1988) 663 - 668.
C       5. JA Potton, GJ Daniell, & BD Rainford,
C               J APPL CRYST 21 (1988) 891 - 897.

C       This progam was written in BASIC by GJ Daniell and later
C         translated into FORTRAN and adapted for SANS analysis.  It
C         has been further modified by AJ Allen to allow use with a
C         choice of particle form factors for different shapes.  It
C         was then modified by PR Jemian to allow portability between
C         the Digital Equipment Corporation VAX and Apple Macintosh
C         computers.
C       The input data file format is three columns of "Q I dI" which
C         are separated by spaces or tabs.  There is no header line
C         in the input data file.

        PARAMETER (cm2m = 0.01) ! convert cm to m units, but why?
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
        DIMENSION r(MaxBin), f(MaxBin), base(MaxBin), dNdr(MaxBin)
        DIMENSION fit(MaxPts), BinWid(MaxPts)
        CHARACTER*40 InFile, OutFil
        LOGICAL Yes
        CHARACTER*1 YN, aTab

        DATA one, zero /1.0, 0.0/       ! compiler-independence!
        DATA hrDamp /8.0/       ! model 7: sets transition range
        DATA htDamp /0.9/       ! model 7: amount of damping
C  The value "hrDamp" sets the range where the transistion occurs.
C  The value "htDamp" sets the maximum proportion of damping.

C ... Define (initially) the default responses
        DATA iOption    /4/     ! usual form factor for spheres
        DATA Aspect     /1.0/   ! particle aspect ratio
        DATA LinLog     /isLin/ ! linear binning scale
        DATA n          /40/    ! number of bins
        DATA Dmin, Dmax /8.00, 400.0/   ! particle diameters
        DATA IterMax    /20/    ! maximum number of iterations to try
        DATA RhoSq      /1.0/   ! scattering contrast, x10**28 1/m**4
        DATA fac, err   /1.0, 1.0/      ! scalars for intensity and errors
        DATA qMin, qMax /1.e-8, 100./   ! range to accept
        DATA Bkg        /0.0/   ! intensity to subtract
        DATA sLengt     /100.0/ ! rectangular slit-length, 1/A

C  Next line for MPW/Language Systems version 1.2.1, Macintosh only
C  Comment this out for other compilers
C  This is the only compiler-dependent line in this source code!!!!!!
C       CALL OutWindowScroll (1000)  ! for 1-line advance screen

        pi = 4. * ATAN(1.)
        aTab = CHAR (9)

C  screen dimension variables for plots, in COMMON /space4/
        nCol = MaxCol
        nRow = MaxRow
        nCol2 = MxC2
        nRow2 = MxR2

    1   WRITE (*,*)
        WRITE (*,*) 'Size distributions from SAS data using the',
     >              ' maximum entropy criterion'
        WRITE (*,*) '       version: ', ProgVers
        WRITE (*,*) '   Last edited: ', EditDate

        CALL GetInf (InFile, OutFil, iOption, Aspect, LinLog,
     >          n, Dmin, Dmax, IterMax, RhoSq, fac, err, qMin,
     >          qMax, Bkg, sLengt)
        IF (InFile .EQ. ' ') STOP

C   Read in the SAS data from the file "InFile"
        WRITE (*,*)  ' Reading from file: ', InFile
        OPEN (UNIT = ioUnit, FILE = InFile, STATUS = 'old')
        DO 94  j = 1, MaxPts
          READ (ioUnit, *, END = 95) q(j), datum(j), sigma(j)
   94   CONTINUE
   95   npt=j-1         ! ignore any lines without an explicit EOL mark
        CLOSE (UNIT = ioUnit, STATUS = 'keep') 
        WRITE (*,*) npt, ' points were read from the file'

C   Subtract background, convert to 1/m units and
C       shift for the selected data range
        i = 0
        DO 2  j = 1, npt
          IF (q(j) .GE. Qmin .AND. q(j) .LE. Qmax) THEN
            i = i + 1
            q(i) = q(j)
            datum(i) = fac * (datum(j)-Bkg) / cm2m
            sigma(i) = fac * err * sigma(j) / cm2m
          END IF
    2   CONTINUE
        npt = i
        WRITE (*,*) npt, ' points were selected from the data'

C  PRJ: 24 May 1989
C       BinWid: actual radial width of the indexed bin number
C       Step:   radial increment factor (for geometric series)
C       rWid:   radial width (for algebraic series)
        IF (LinLog .EQ. isLog) THEN     ! geometric series
          Step = (Dmax/Dmin)**(1. / FLOAT(n-1)) - 1.
          rWid = 0.
        ELSE                            ! algebraic series
          Step = 0.
          rWid = 0.5*(Dmax - Dmin) / FLOAT(n-1)
        END IF
        r(1) = 0.5 * Dmin
        BinWid(1) = r(1) * Step + rWid
        DO 48  i = 2, n
          r(i) = r(i-1) + BinWid(i-1)
          BinWid(i) = r(i) * Step + rWid
   48   CONTINUE

        WRITE (*,*) ' Preparation of the GRID function...'
C  Calculate the form-factor pre-terms
  111   IF (iOption .EQ. 1) THEN        ! Rods, using model of AJ Allen
          alphan1 = 2. * pi * Aspect
          alphan2 = 4. * pi
          preform = alphan1
          sLengt = 0.                   ! "pinhole" collimation
        ELSE IF (iOption .EQ. 2) THEN   ! Disks, using model of AJ Allen
          alphan1 = 2. * pi / (Aspect**2)
          alphan2 = 2. * pi
          preform = alphan1
          sLengt = zero
        ELSE IF (iOption .EQ. 3) THEN   ! Globules, using model of AJ Allen
          alphan1 = 4. * pi * Aspect / 3.
          IF (Aspect .LT. 0.99) THEN    ! hamburger-shaped
            sqqt = SQRT (one - Aspect**2)
            argument = (2. - Aspect**2 + 2. * sqqt) / (Aspect**2)
            surchi = (one + Aspect**2 * LOG(argument) / (2.*sqqt) )
     >          / (2. * Aspect)
          ELSE IF (Aspect .GT. 1.01) THEN       ! peanut shaped
            sqqt = SQRT(Aspect**2 - one)
            argument = sqqt / Aspect
            surchi = (one + Aspect**2 * ASIN(argument) / sqqt)
     >          / (2. * Aspect)
          ELSE                          ! spheroidal
            surchi = one
          END IF
          alphan2 = 6. * pi * surchi
          preform = alphan1
          sLengt = zero
        ELSE IF (iOption .EQ. 4) THEN   ! Spheres, delta-function
          alphan1 = 4. * pi / 3.
          alphan2 = 6. * pi
          preform = 9. * alphan1
          sLengt = zero
        ELSE IF (iOption .EQ. 5) THEN   ! Spheres, box-distribution
          alphan1 = 4. * pi / 3.        ! This model by PRJ
          alphan2 = 6. * pi
          preform = 48. * pi
          sLengt = zero
        ELSE IF (iOption .EQ. 6) THEN   ! smeared, spheroidal globs
          preform = 4. * Pi / 3.        ! This model by PRJ
          alphan1 = preform
          alphan2 = 6. * Pi
          Cgs = SQRT (3. * Pi)          ! for low-Q region
          Cps = 9. * Pi / 4.            ! for med. high-Q region
          Cp = 9. / 2.                  ! for high-Q region
        ELSE IF (iOption .EQ. 7) THEN   ! spheroidal globs, no smearing
          preform = 4. * Pi / 3.        ! This model by PRJ
          alphan1 = preform
          alphan2 = 6. * Pi
          sLengt = zero
        END IF

C  alphaN1 is RhoSq * the particle volume
C  alphaN2 is RhoSq * the particle surface area / the particle volume
C       ... and later divided by q**4
        alphan1 = cm2m * alphan1 * rhosq * r(1)**3
        alphan2 = cm2m * alphan2 * rhosq / r(n)
        preform = cm2m * preform * rhosq

        DO 226  i = 1, n
          rCubed = r(i)**3
          DO 226  j = 1, npt
            Qr = q(j) * r(i)
            IF (iOption .EQ. 1) THEN
              QH = q(j) * Aspect * r(i)         ! rod 1/2 - length
              topp = one + 2.*Pi* QH**3 * Qr / (9 * (4 + Qr**2))
     >                   + (QH**3 * Qr**4) / 8.
              bott = one + QH**2 * (one + QH**2 * Qr)/9
     >                   + (QH**4 * Qr**7) / 16
            ELSE IF (iOption .EQ. 2) THEN
              h = r(i)                  ! disk 1/2 - thickness
              Rd = h / Aspect           ! disk radius
              Qh = q(j) * h
              QRd = q(j) * Rd
              topp = one + QRd**3 / (3. + Qh**2)
     >                   + (Qh**2 * QRd / 3.)**2
              bott = one + QRd**2 * (one + Qh * QRd**2) / 16
     >                   + (Qh**3 * QRd**2 / 3.)**2
            ELSE IF (iOption .EQ. 3) THEN
              topp = one
              bott = one + Qr**2 * (2. + Aspect**2) / 15.
     >                   + 2. * Aspect * Qr**4 / (9. * surchi)
            ELSE IF (iOption .EQ. 4) THEN
              topp = (SIN(Qr) - Qr * COS(Qr))**2
              bott = Qr**6
            ELSE IF (iOption .EQ. 5) THEN
              Qj = q(j)
              rP = r(i) + BinWid(i)
              rM = r(i)
              bP = 0.5*rP + (Qj**2)*(rP**3)/6.
     >          + (0.25*(Qj * rP**2) - 0.625/Qj) * SIN (2.*Qj*rP)
     >          + 0.75 * rP * COS (2.*Qj*rP)
              bM = 0.5*rM + (Qj**2)*(rM**3)/6.
     >          + (0.25*(Qj * rM**2) - 0.625/Qj) * SIN (2.*Qj*rM)
     >          + 0.75 * rM * COS (2.*Qj*rM)
              topp = bP - bM
              bott = Qj**6 * (rP**4 - rM**4) * rCubed
            ELSE IF (iOption .EQ. 6) THEN
              rL = r(i) * sLengt
              topp = Cgs
              bott = rL*(one + (Qr**2)/5. + Cgs/Cps * Qr**3)
     >                  + Cgs/Cp * Qr**4
            ELSE IF (iOption .EQ. 7) THEN
C  The value "hrDamp" sets the range where the transistion occurs.
C  The value "htDamp" sets the maximum proportion of damping.
C  The weight is a "step" function with a broad edge.
              weight = htDamp * EXP (-((Qr/hrDamp)**2)) + (one - htDamp)
              topp = 3. * (SIN(Qr) - Qr * COS(Qr)) / Qr**3
              bott = 4.5 / Qr**4        ! bott=<topp**2> for large Qr
              topp = weight * topp**2 + (one-weight) / (one + one/bott)
              bott = one
            END IF
            grid(i,j) = preform * rCubed * topp / bott
C               factors of 4Pi/3 are already included in "preform"
  226   CONTINUE

C       Attempt to account for scattering from very large and very small 
C       particles by use of the limiting forms of grid(i,j).
        DO 227  j = 1, npt
          grid(n+1,j) = alphan1 ! next line accounts for a slit-length
          grid(n+2,j) = alphan2 / (q(j)**3 * SQRT(q(j)**2 + sLengt**2))
  227   CONTINUE

C  Try to solve the problem
  228   basis = 1.0e-12 / RhoSq         ! Originally was 1.0e-6
        CALL MaxEnt (n+2,npt, f,datum,sigma, basis,base, max,itermax)

C       "Max" counts the number of iterations inside MAXENT.
C       If Max < IterMax, then the problem has been solved.
        IF (max .GE. itermax) THEN
          WRITE (*,*) ' No convergence! # iter. = ', max
          WRITE (*,*) ' File was: ', InFile
          GO TO 1
        END IF

C   Otherwise, SUCCESS!... so calculate the volume distribution
C       from the model SAS data
        CALL opus (n+2, npt, f, fit)

C   ... and remove the bin width effect.
C   Also, calculate the total volume fraction, the mode, mean, and
C       standard deviations of the volume and number distributions.
        SumV   = zero
        SumVR  = zero
        SumVR2 = zero
        SumN   = zero
        SumNR  = zero
        SumNR2 = zero
        modeV = 1
        modeN  = 1
        DO 1919  i = 1, n
          size = r(i)
          frac = f(i)
          pVol = 4*Pi/3 * (size * 1.e-8)**3     ! particle volume, cm**3
          IF (iOption .EQ. 1)  pVol = pVol * Aspect     ! rods
          IF (iOption .EQ. 2)  pVol = pVol / Aspect     ! disks
          IF (iOption .EQ. 3)  pVol = pVol * Aspect     ! globs
          amount = frac / pVol                  ! number / cm**3
          f(i) = frac / BinWid(i)
          dNdr(i) = amount / BinWid(i)
          IF (i .GT. 3) THEN                    ! ignore 1st few bins
            SumN   = SumN   + amount
            SumNR  = SumNR  + amount * size
            SumNR2 = SumNR2 + amount * size**2
          END IF
          IF (dNdr(i) .GT. dNdr(modeN)) modeN = i       ! get the mode
          SumV   = SumV   + frac
          SumVR  = SumVR  + frac * size
          SumVR2 = SumVR2 + frac * size**2
          IF (f(i) .GT. f(modeV)) modeV = i             ! get the mode
 1919   CONTINUE
        DnMean = 2.0 * SumNR / SumN
        DnSDev = 2.0 * SQRT((SumNR2 / SumN) - (SumNR / SumN)**2)
        DvMean = 2.0 * SumVR / SumV
        DvSDev = 2.0 * SQRT((SumVR2 / SumV) - (SumVR / SumV)**2)

        Entropy = zero
        DO 1920  i = 1, n
          frac = BinWid(i) * f(i) / SumV        ! Skilling & Bryan, eq. 1
          Entropy = Entropy - frac * LOG (frac)
 1920   CONTINUE

C  Show the final distribution, corrected for bin width.

        WRITE (*,*)
        WRITE (*,*) ' Input file: ', InFile
        WRITE (*,*) ' Volume weighted size dist.: V(r)N(r) versus r'
        CALL Plot (n, r, f)

C   Estimate a residual background that remains in the data.
        Sum1 = zero
        Sum2 = zero
        DO 918 j = 1, npt
          weight = one / (sigma(j)**2)
          Sum1 = Sum1 + weight * (fit(j) -  datum(j))
          Sum2 = Sum2 + weight
  918   CONTINUE
        shift = Sum1 / Sum2

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
  919   CONTINUE
        shift = cm2m * shift / fac

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
        WRITE (ioUnit,*) '    of size between D and D + dD'
        WRITE (ioUnit,*)
        WRITE (ioUnit, 35591) 'D, A', aTab, 'V(D)*N(D), 1/A',
     >                  aTab, 'N(D), 1/A/cm^3'
        WRITE (ioUnit, 35591) '----', aTab, '--------------',
     >                  aTab, '--------------'
35591   FORMAT (1X, A12, A1, 1X, A15, A1, 1X, A15)

        DO 1001 i = 1, n
 1001     WRITE (ioUnit,3559) 2.*r(i), aTab, 0.5*f(i), aTab, 0.5*dNdr(i)
 3559   FORMAT (1X, F12.2, A1, 1X, 1PE15.5, A1, 1X, E15.5)


        WRITE (ioUnit,'(///)')
        WRITE (ioUnit, 1011) 'Q 1/A', aTab, 'I 1/cm', aTab,
     >                  '^I 1/cm', aTab, 'dI 1/cm', aTab, 'z'
        WRITE (ioUnit, 1011) '-----', aTab, '------', aTab,
     >                  '-------', aTab, '-------', aTab, '----'
 1011   FORMAT (A12, 3(A1, A12), 1X, A1, A12, 1X, A1, A12)

        DO 101  j = 1, npt
  101     WRITE (ioUnit,560) q(j), aTab, datum(j), aTab, fit(j),
     >          aTab, sigma(j), aTab, z(j)
  560     FORMAT (1PE12.4, 3(A1, E12.4), 1X, A1, 0PF12.6, 1X, A1, F12.6)

        WRITE (ioUnit,3301) InFile
        WRITE (*,3301) InFile
 3301   FORMAT (//' Input data: ', A40)

        WRITE (ioUnit,3302) RhoSq
        WRITE (*,3302) RhoSq
 3302   FORMAT (' Contrast = ', F15.7,' x 10^28 m^-4.')

        IF (iOption .EQ. 1) THEN
          WRITE (ioUnit,*) ' rods: dia=D, length=D*', Aspect
          WRITE (*,*) ' rods: dia=D, length=D*', Aspect
        ELSE IF (iOption .EQ. 2) THEN
          WRITE (ioUnit,*) ' disks: thickness=D, dia=D/', Aspect
          WRITE (*,*) ' disks: thickness=D, dia=D/', Aspect
        ELSE IF (iOption .EQ. 3) THEN
          WRITE (ioUnit,*) ' globs: D x D x D*', Aspect
          WRITE (*,*) ' globs: D x D x D*', Aspect
        ELSE IF (iOption .EQ. 4) THEN
          WRITE (ioUnit,*) ' delta-function Spheres: diameter=D'
          WRITE (*,*) ' delta-function Spheres: diameter=D'
        ELSE IF (iOption .EQ. 5) THEN
          WRITE (ioUnit,*) ' box-function Spheres: diameter=D'
          WRITE (*,*) ' box-function Spheres: diameter=D'
        ELSE IF (iOption .EQ. 6) THEN
          WRITE (ioUnit,*) ' slit-smeared spheroidal globs: diameter=D'
          WRITE (*,*) ' slit-smeared spheroidal globs: diameter=D'
          WRITE (ioUnit,*) ' slit-length (1/A) = ', sLengt
          WRITE (*,*) ' slit-length (1/A) = ', sLengt
        ELSE IF (iOption .EQ. 7) THEN
          WRITE (ioUnit,*) ' spheroidal globs: diameter=D'
          WRITE (*,*) ' spheroidal globs: diameter=D'
        END IF

        WRITE (ioUnit,53303) fac
        WRITE (*,53303) fac
53303   FORMAT (' Data conversion factor to 1/cm = ', 1PE12.5)

        WRITE (ioUnit,63303) err
        WRITE (*,63303) err
63303   FORMAT (' Error scaling factor = ', 1PE12.5)

        IF (LinLog .EQ. isLog) THEN
          WRITE (ioUnit,13304) 'geometric'
          WRITE (*,13304) 'geometric'
        ELSE
          WRITE (ioUnit,13304) 'algebraic'
          WRITE (*,13304) 'algebraic'
        END IF
13304   FORMAT (' Histogram bins are distributed in an increasing ',
     >          A9, ' series.')

        WRITE (ioUnit,3304) 'Minimum', Dmin
        WRITE (*,3304) 'Minimum', Dmin
        WRITE (ioUnit,3304) 'Maximum', Dmax
        WRITE (*,3304) 'Maximum', Dmax
 3304   FORMAT (1X, A7, ' particle dimension D = ',F12.2,' A.')

        WRITE (ioUnit,3306) n
        WRITE (*,3306) n
 3306   FORMAT (' Number of histogram bins = ',I4,'.')

        WRITE (ioUnit,3307) itermax
        WRITE (*,3307) itermax
 3307   FORMAT (' Maximum number of iterations allowed = ',I4,'.')

        WRITE (ioUnit,3314) max
        WRITE (*,3314) max
 3314   FORMAT (' Program left MaxEnt routine after ',
     *    I4,' iterations.')

        WRITE (ioUnit,3308) npt
        WRITE (*,3308) npt
 3308   FORMAT (' Target chi-squared (# data points) = ',I5,'.')

        WRITE (ioUnit,3309) ChiSq
        WRITE (*,3309) ChiSq
 3309   FORMAT (' Best value of chi-squared achieved = ',F12.6,'.')

        WRITE (ioUnit, 33091) 'the final', Entropy
        WRITE (*, 33091) 'the final', Entropy
        WRITE (ioUnit, 33091) 'a flat', LOG (FLOAT (n))
        WRITE (*, 33091) 'a flat', LOG (FLOAT (n))
33091   FORMAT (' Entropy of ', A9, ' distribution = ', F12.7,'.')

        WRITE (ioUnit,33101) SumN
        WRITE (*,33101) SumN
33101   FORMAT (' Total particles  = ', 1PE15.5,' per cubic cm.')

        WRITE (ioUnit,3310) SumV
        WRITE (*,3310) SumV
 3310   FORMAT (' Total volume fraction of all scatterers = ',
     *    F15.9,'.')

        WRITE (ioUnit,3311) 'smaller', Dmin, f(n+1)
        WRITE (ioUnit,3311) 'larger',  Dmax, f(n+2)
        WRITE (*,3311) 'smaller', Dmin, f(n+1)
        WRITE (*,3311) 'larger',  Dmax, f(n+2)
 3311   FORMAT (' Volume fraction ',A7,' than ', F12.2,
     *    ' A = ', 1PE13.5,'.')

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
 3312   FORMAT (1X, A6, '-weighted ', A14, ' = ', F12.5, ' A.')

        WRITE (ioUnit,3313) 'Min', q(1)
        WRITE (*,3313) 'Min', q(1)
        WRITE (ioUnit,3313) 'Max', q(npt)
        WRITE (*,3313) 'Max', q(npt)
 3313   FORMAT (1X, A3,'imum Q-vector = ', 1PE15.7, ' 1/A.')

        WRITE (ioUnit,3315) 'User-specified', Bkg
        WRITE (*,3315) 'User-specified', Bkg
        WRITE (ioUnit,3315) 'Suggested', Bkg - shift
        WRITE (*,3315) 'Suggested', Bkg - shift
 3315   FORMAT (1X, A14, ' background = ', F18.9,' input data units')

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
 4010   CONTINUE

        IF (ABS ((Chi2Bk-ChiSq)/FLOAT (npt)) .LE. 0.05) THEN
          WRITE (*,*) ' The change in ChiSquared should be < 5%.'
 4000     WRITE (*,*) ' Run the Stability Check? (Y/<N>)'
          READ (*,'(A1)') YN
          IF (YN .EQ. 'y'  .OR.  YN .EQ. 'Y') GO TO 228
          IF (YN.NE.' ' .AND. YN.NE.'n' .AND. YN.NE.'N') GO TO 4000
        END IF

        WRITE (*,3200) OutFil
 3200   FORMAT (/,' The program is finished.', /, 
     1    ' The output file is: ', A40)
        GO TO 1

 3199   STOP 
        END


        SUBROUTINE GetInf (InFile, OutFil, iOption, Aspect, LinLog,
     >          nBin, Dmin, Dmax, IterMax, RhoSq, fac, err, qMin,
     >          qMax, Bkg, sLengt)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        CHARACTER*40 InFile, OutFil
        PARAMETER (Ro2Max = 1.e6, ItrLim = 200, AbsMax = 1.e3)
        PARAMETER (DiaMin = 1., DiaMax = 1.e6, ErrMax = 1.e6)
        PARAMETER (MaxPts = 300, MaxBin = 102)
        PARAMETER (isLin = 1, isLog = 2)

    1   WRITE (*,*) ' Input file? <Quit>'
          READ (*, 2) InFile
    2     FORMAT (A40)
          IF (InFile.EQ.' ') RETURN

    3   WRITE (*,*) ' Output file?'
          READ (*, 2) OutFil
          IF (OutFil .EQ. ' ') GO TO 3
          IF (OutFil .EQ. InFile) GO TO 1

        suggest = qMin
   16   WRITE (*,*) ' Minimum q-vector? [1/A] <', suggest, '>'
        READ (*, '(F10.0)') qMin
        IF (qMin .LT. 0) GO TO 16
        IF (qMin .EQ. 0) qMin = suggest

        suggest = qMax
   17   WRITE (*,*) ' Maximum q-vector? [1/A] <', suggest, '>'
        READ (*, '(F10.0)') qMax
        IF (qMax .EQ. 0) qMax = suggest
        IF (qMax .LE. 0) GO TO 17
        IF (qMax .LE. qMin) GO TO 1

        suggest = RhoSq
   13   WRITE (*,*) ' Scattering contrast? [10^28 m^-4] <',suggest,'>'
        READ (*, '(F10.0)') RhoSq
        IF (RhoSq .EQ. 0) RhoSq = suggest
        IF (RhoSq .LT. 0 .OR. RhoSq .GT. Ro2Max) GO TO 13

        suggest = fac
   14   WRITE (*,*) ' Factor to convert data to 1/cm? <', suggest, '>'
        READ (*, '(F10.0)') fac
        IF (fac .EQ. 0) fac = suggest
        IF (fac .LE. 0 .OR. fac .GT. AbsMax) GO TO 14

        suggest = err
   15   WRITE (*,*) ' Error scaling factor? <', suggest, '>'
        READ (*, '(F10.0)') err
        IF (err .EQ. 0) err = suggest
        IF (err .LE. 0 .OR. err .GT. ErrMax) GO TO 15

        suggest = Bkg
   18   WRITE (*,*) ' Background? <', suggest, '>'
        READ (*, '(F10.0)') Bkg
        IF (Bkg .EQ. 0) Bkg = suggest

        Last = iOption
    4     WRITE (*,*) '  Select a form model for the scatterer:'
          WRITE (*,*) '  (See the User Guide for complete explanations)'
          WRITE (*,*) ' 1: rods        2: disks       3: globules'
          WRITE (*,*) ' 4: spheres (usual form)       ',
     >                  '5: spheres (integrated)'
          WRITE (*,*) ' 6: spheroids (slit-smeared)   ',
     >                  '7: spheroids (not smeared)'
          WRITE (*,*) '  Which option number?  <', Last, '>'
          READ (*, '(I4)') iOption
          IF (iOption .EQ. 0) iOption = Last
          IF (iOption .LT. 1  .OR.  iOption .GT. 7) GO TO 4

        suggest = Aspect
    6   IF (iOption .GE. 1  .AND.  iOption .LE. 3) THEN
          WRITE (*,*) ' AR = Aspect Ratio, useful ranges are indicated'
          IF (iOption .EQ. 1) THEN
            WRITE (*,*) ' diameter D, length D * AR, AR > 5'
          ELSE IF (iOption .EQ. 2) THEN
            WRITE (*,*) ' thickness D, diameter D / AR, AR < 0.2'
          ELSE IF (iOption .EQ. 3) THEN
            WRITE (*,*) ' D x D x D * AR, 0.3 < AR < 3'
          END IF
          WRITE (*,*) ' Aspect ratio?  <', suggest, '>'
          READ (*,'(F10.0)') Aspect
          IF (Aspect .EQ. 0) Aspect = suggest
          IF (Aspect .LT. 0) GO TO 6
        END IF

        suggest = sLengt
   61   IF (iOption .EQ. 6) THEN
          WRITE (*,*) ' Slit-smeared globs.  ',
     >          'Slit-length [1/A]? <', suggest, '>'
          READ (*,'(F10.0)') sLengt
          IF (sLengt .EQ. 0) sLengt = suggest
          IF (sLengt .LT. 0) GO TO 61
        END IF

        Last = LinLog
    7   WRITE (*,*) ' Bin step scale? (1=Linear, 2=Log) <', Last, '>'
        READ (*, '(I4)') LinLog
        IF (LinLog .EQ. 0) LinLog = Last
        IF (LinLog .NE. isLin  .AND. LinLog .NE. isLog) GO TO 7

        Last = nBin
    8   WRITE (*,*) ' Number of histogram bins? <', Last, '>'
        READ (*, '(I4)') nBin
        IF (nBin .EQ. 0) nBin = Last
        IF (nBin .LT. 2 .OR. nBin .GT. (MaxBin-2)) GO TO 8

        suggest = Dmax
    9   WRITE (*,*) ' Maximum value of D? [A] <', suggest, '>'
        READ (*, '(F10.0)') Dmax
        IF (Dmax .EQ. 0) Dmax = suggest
        IF (Dmax .LT. nBin*DiaMin .OR. Dmax .GE. DiaMax) GO TO 9

        Suggest = Dmax / FLOAT (nBin)
   11   WRITE (*,*) ' Minimum value of D? [A] <', suggest, '>'
        READ (*, '(F10.0)') Dmin
        IF (Dmin .EQ. 0) Dmin = suggest
        IF (Dmin .GE. DMax .OR. Dmin .LT. DiaMin) GO TO 1

        IF (IterMax .GT. ItrLim) IterMax = ItrLim
        Last = IterMax
   12   WRITE (*,*) ' Maximum number of iterations? <', Last, '>'
        READ (*, '(I4)') IterMax
        IF (IterMax .EQ. 0) IterMax = Last
        IF (IterMax .LT. 0 .OR. IterMax .GT. ItrLim) GO TO 12

        RETURN
        END


        SUBROUTINE opus(n,npt,x,ox)     ! solution-space -> data-space
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        PARAMETER (MaxPts=300, MaxBin=102)
        COMMON /space1/ grid
        DIMENSION x(MaxBin), grid(MaxBin,MaxPts), ox(MaxPts)
        DO 3  j = 1, npt
          sum = 0.
          DO 4  i = 1, n
           sum = sum + x(i) * grid(i,j)
    4     CONTINUE
          ox(j) = sum
    3   CONTINUE
        RETURN
        END


        SUBROUTINE tropus(n,npt,ox,x)   ! data-space -> solution-space
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        PARAMETER (MaxPts=300, MaxBin=102)
        COMMON /space1/ grid
        DIMENSION x(MaxBin), grid(MaxBin,MaxPts), ox(MaxPts)
        DO 5  i = 1, n
          sum = 0.
          DO 6  j = 1, npt
            sum = sum + ox(j) * grid(i,j)
    6     CONTINUE
          x(i) = sum
    5   CONTINUE
        RETURN
        END


        SUBROUTINE MaxEnt(n,npt, f,datum,sigma, flat,base,iter,itermax)
        IMPLICIT REAL*8 (A-H,O-Z)
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
        PARAMETER (TstLim = 0.05)       ! for convergence
        DATA one, zero /1.0, 0.0/       ! compiler-independence!

        blank = flat
        exp1 = EXP(one)

        IF (blank .EQ. zero) THEN
          DO 1004 i = 1, n
 1004       blank = blank + base(i)
          blank = blank / FLOAT(n)
          WRITE (*,*) ' Average of BASE = ', blank
        ELSE
          WRITE (*,*) ' Setting BASE constant at ', blank
          DO 1003 i = 1, n
 1003       base(i) = blank
        ENDIF

        WRITE (*,*) ' MaxEnt routine beginning ...'

        chizer = FLOAT(npt)
        chtarg = chizer
        m = 3
        DO 8 i = 1, n
    8     f(i) = base(i)        ! initial distribution is featureless

        iter = 0
    6   iter = iter + 1         ! The iteration loop begins here!
        CALL opus (n, npt, f, ox)       ! calc. the model intensity from "f"
        chisq = zero
        DO 10 j = 1, npt
          a = (ox(j) - datum(j)) / sigma(j)
          chisq = chisq + a**2
   10     ox(j) = 2. * a / sigma(j)
        CALL tropus(n,npt,ox,cgrad)     ! cGradient = Grid * ox
        test = zero     ! mismatch between entropy and ChiSquared gradients
        snorm = zero    ! entropy term
        cnorm = zero    ! ChiSqr term
        tnorm = zero    ! norm for the gradient term TEST
        fSum = zero     ! find the sum of the f-vector
        DO 12  i = 1, n
          fSum = fSum + f(i)
          sgrad(i) = -LOG(f(i)/base(i)) / (blank*exp1)
          snorm = snorm + f(i) * sgrad(i)**2
          cnorm = cnorm + f(i) * cgrad(i)**2
          tnorm = tnorm + f(i) * sgrad(i) * cgrad(i)
   12   CONTINUE
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
   13   CONTINUE
        CALL opus (n,npt,xi(1,1),eta(1,1))
        CALL opus (n,npt,xi(1,2),eta(1,2))
        DO 14 j = 1, npt
          ox(j) = eta(j,2) / (sigma(j)**2)
   14   CONTINUE
        CALL tropus (n,npt,ox,xi(1,3))
        a = zero
        DO 15 i = 1, n
          b = f(i) * xi(i,3)
          a = a + b * xi(i,3)
          xi(i,3) = b
   15   CONTINUE
        a = one / SQRT(a)
        DO 16 i = 1, n
          xi(i,3) = a * xi(i,3)
   16   CONTINUE
        CALL opus (n,npt,xi(1,3),eta(1,3))
        DO 17 k = 1, m
          s1(k) = zero
          c1(k) = zero
          DO 18 i = 1, n
            s1(k) = s1(k) + xi(i,k) * sgrad(i)
            c1(k) = c1(k) + xi(i,k) * cgrad(i)
   18     CONTINUE
          c1(k) = c1(k) / chisq
   17   CONTINUE
        DO 19 k = 1, m
          DO 19 l = 1, k
            s2(k,l) = zero
            c2(k,l) = zero
            DO 20 i = 1, n
              s2(k,l) = s2(k,l) - xi(i,k) * xi(i,l) / f(i)
   20       CONTINUE
            DO 21 j = 1, npt
              c2(k,l) = c2(k,l) + eta(j,k) * eta(j,l) / (sigma(j)**2)
   21       CONTINUE
            s2(k,l) = s2(k,l) / blank
            c2(k,l) = 2. * c2(k,l) / chisq
   19   CONTINUE
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
        fSum = zero             ! find the sum of the f-vector
        fChange = zero          ! and how much did it change?
        DO 23 i = 1, n
          df = beta(1)*xi(i,1)+beta(2)*xi(i,2)+beta(3)*xi(i,3)
          IF (df .LT. -f(i)) df = 0.001 * base(i) - f(i)        ! a patch
          f(i) = f(i) + df              ! adjust the f-vector
          fSum = fSum + f(i)
          fChange = fChange + df
   23   CONTINUE

        s = zero
        DO 24  i = 1, n
          temp = f(i) / fSum            ! fraction of f(i) in this bin
          s = s - temp * LOG (temp)     ! from Skilling and Bryan, eq. 1
   24   CONTINUE

        CALL opus (n, nPt, f, z)        ! model the data-space from f(*)
        ChiSq = zero                    ! get the new ChiSquared
        DO 25  j = 1, nPt
          z(j) = (datum(j) - z(j)) / sigma(j)   ! the residuals
          ChiSq = ChiSq + z(j)**2       ! report this ChiSq, not the one above
   25   CONTINUE

        Entropy(iter) = s
        Convrg(iter) = LOG (ChiSq)
        IF (iter .GT. 2) THEN   ! show our progress
          temp = (Convrg(iter) + Convrg(iter-1) + Convrg (iter-2))/3.
          IF (ABS (one - Convrg(iter)/temp) .GT. 0.02) THEN
            WRITE (*,*)
            WRITE (*,*) ' LOG (ChiSq) vs. iteration number'
            CALL BasPlt (iter, Convrg, LOG (ChiZer))
          END IF
          WRITE (*,*)
          WRITE (*,*) ' Entropy vs. iteration number'
          temp = LOG (FLOAT (n))        ! the maximum entropy possible
          CALL BasPlt (iter, Entropy, temp)
        END IF

  300   WRITE (*,*)
        WRITE (*,*) ' Residuals'
        CALL ResPlt (npt, z)

        WRITE (*,*)
        WRITE (*,*) ' Distribution'
        CALL BasPlt (n, f, blank)

        WRITE (*,*) ' #', iter, ' of ', itermax, ',  n  = ', npt
        WRITE (*,200) test, s
        WRITE (*,201) 'target',SQRT(chtarg/npt), 'now',SQRT(chisq/npt)
        WRITE (*,202) 'sum', fSum, ' % change', 100.*fChange/fSum
  200   FORMAT (' test = ', F9.5, ',  Entropy = ', F12.7)
  201   FORMAT (' SQRT((Chi^2)/n):', A8,' = ', F12.8,A10,' = ', F12.8)
  202   FORMAT ('        f-vector:', A8,' = ', F12.8,A10,' = ', F12.8)

C  See if we have finished our task.
        IF (ABS(chisq/chizer-one) .LT. 0.01) THEN  ! hardest test first
          IF (test .LT. TstLim) THEN            ! same solution gradient?
C               We've solved it but now must check for a bizarre condition.
C               Calling routine says we failed if "iter = iterMax".
C               Let's increment (maybe) iterMax so this doesn't happen.
            IF (iter .EQ. iterMax) iterMax = iterMax + 1
            RETURN
          END IF
        END IF
        IF (iter .LT. iterMax) GO TO 6

C  Ask for more time to finish the job.
        WRITE (*,*)
        WRITE (*,*) ' Maximum iterations have been reached.'
 2001   WRITE (*,*) ' How many more iterations? <none>'
        READ (*,'(I4)') more
        IF (more .LT. 0) GO TO 2001
        IF (more .EQ. 0) RETURN
        iterMax = iterMax + more
        GO TO 6
        END


        SUBROUTINE Move(m)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        PARAMETER ( MxLoop = 500 )      ! for no solution
        PARAMETER ( Passes = 1.e-3 )    ! convergence test
        COMMON /space5/ chisq, chtarg, chizer, fSum, blank
        COMMON /space2/ beta, c1, c2, s1, s2
        DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
        DATA one, zero /1.0, 0.0/       ! compiler-independence!
        a1 = zero                       ! lower bracket  "a"
        a2 = one                        ! upper bracket of "a"
        cmin = ChiNow (a1, m)
        IF (cmin*chisq .GT. chizer) ctarg = 0.5*(one + cmin)
        IF (cmin*chisq .LE. chizer) ctarg = chizer/chisq
        f1 = cmin - ctarg
        f2 = ChiNow (a2,m) - ctarg
        DO 1  loop = 1, MxLoop
          anew = 0.5 * (a1+a2)          ! choose a new "a"
          fx = ChiNow (anew,m) - ctarg
          IF (f1*fx .GT. zero) a1 = anew
          IF (f1*fx .GT. zero) f1 = fx
          IF (f2*fx .GT. zero) a2 = anew
          IF (f2*fx .GT. zero) f2 = fx
          IF (abs(fx) .LT. Passes) GO TO 2
    1   CONTINUE

C  If the preceding loop finishes, then we do not seem to be converging.
C       Stop gracefully because not every computer uses control-C (etc.)
C       as an exit procedure.
        WRITE (*,*) ' Loop counter = ', MxLoop
        PAUSE ' No convergence in alpha chop (MOVE).  Press return ...'
        STOP ' Program cannot continue.'

    2   w = Dist (m)
        IF (w .LE. 0.1*fSum/blank) GO TO 1042
        DO 1044 k=1,m
          beta(k) = beta(k) * SQRT(0.1 * fSum/(blank * w))
 1044   CONTINUE
 1042   chtarg = ctarg * chisq
        RETURN
        END


        REAL*8 FUNCTION Dist (m)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        COMMON /space5/ chisq, chtarg, chizer, fSum, blank
        COMMON /space2/ beta, c1, c2, s1, s2
        DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
        DATA one, zero /1.0, 0.0/       ! compiler-independence!
        w = zero
        DO 26  k = 1, m
          z = zero
          DO 27  l = 1, m
            z = z - s2(k,l) * beta(l)
   27     CONTINUE
          w = w + beta(k) * z
   26   CONTINUE
        Dist = w
        RETURN
        END


        REAL*8 FUNCTION ChiNow(ax,m)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        COMMON /space5/ chisq, chtarg, chizer, fSum, blank
        COMMON /space2/ beta, c1, c2, s1, s2
        DIMENSION beta(3), c1(3), c2(3,3), s1(3), s2(3,3)
        DIMENSION a(3,3), b(3)
        DATA one, zero /1.0, 0.0/       ! compiler-independence!
        bx = one - ax
        DO 28  k = 1, m
          DO 29  l = 1, m
            a(k,l) = bx * c2(k,l)  -  ax * s2(k,l)
   29     CONTINUE
          b(k) = -(bx * c1(k)  -  ax * s1(k))
   28   CONTINUE
        CALL ChoSol(a,b,m,beta)
        w = zero
        DO 31  k = 1, m
          z = zero
          DO 32  l = 1, m
            z = z + c2(k,l) * beta(l)
   32     CONTINUE
          w = w + beta(k) * (c1(k) + 0.5 * z)
   31   CONTINUE
        ChiNow = one +  w
        RETURN
        END


        SUBROUTINE ChoSol(a, b, n, beta)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION fl(3,3), a(3,3), bl(3), b(3), beta(3)
        DATA one, zero /1.0, 0.0/       ! compiler-independence!
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
   36       CONTINUE
            z = a(i,j) - z
            IF (j .EQ. i) fl(i,j) = SQRT(z)
            IF (j .NE. i) fl(i,j) = z / fl(j,j)
35      CONTINUE
        bl(1) = b(1) / fl(1,1)
        DO 37  i=2, n
          z = zero
          DO 38  k = 1, i-1
            z = z + fl(i,k) * bl(k)
   38     CONTINUE
          bl(i) = (b(i) - z) / fl(i,i)
   37   CONTINUE
        beta(n) = bl(n) / fl(n,n)
        DO 39  i1 = 1, n-1
          i = n - i1
          z = zero
          DO 40  k = i+1, n
            z = z + fl(k,i) * beta(k)
   40     CONTINUE
          beta(i) = (bl(i) - z) / fl(i,i)
   39   CONTINUE
        RETURN
        END


        SUBROUTINE ResPlt (n, x)
C       Draw a plot of the standardized residuals on the screen.
C       Mark the rows of + and - one standard deviation.
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION x(1)
        CHARACTER*1 Blank, Symbol, hBordr, vBordr, resSym
        PARAMETER (Blank = ' ', Symbol = 'O', resSym = '=')
        PARAMETER (hBordr = '-', vBordr = '|')
        COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
        CHARACTER*1 screen(100, 150)
        IF (n .LT. 2) RETURN    ! not enough data

C  Find out how many points to pack per column and how many columns
        nPack = 1 + INT(FLOAT (n) / MaxCol - 1./n)
        nCol = INT((n - 1./n)/nPack + 1)

C  prepare the "screen" for drawing
        DO 1  j = 1, nCol + 2
          DO 1  i = 1, MxR2
            screen(i,j) = Blank
    1   CONTINUE
        DO 2  i = 2, nCol + 1
          screen(MxR2,i) = hBordr
    2     screen(1,i) = hBordr
        DO 3  i = 2, MaxRow + 1
          screen(i,nCol+2) = vBordr
    3     screen(i,1) = vBordr

C  get the data limits
        xMax = 1.
        xMin = -1.
        DO 4 i = 1, n
          IF (x(i) .GT. xMax) xMax = x(i)
          IF (x(i) .LT. xMin) xMin = x(i)
    4   CONTINUE
        RowDel = (MaxRow - 1) / (xMax - xMin)

C  show the standard deviation bars
        mPlus = 1 + INT((1 - xMin)*RowDel + 1)
        mMinus = 1 + INT((-1 - xMin)*RowDel + 1)
        DO 5  i = 2, nCol + 1
          screen(mMinus,i) = resSym
    5     screen(mPlus,i) = resSym

C  draw the data (overdrawing the residuals bars if necessary)
        DO 6  i = 1, n
          mCol = 1 + INT((i - 1./n)/nPack + 1)          ! addressing function
          mRow = 1 + INT((x(i) - xMin)*RowDel + 1)      ! +1 for the plot frame
          screen(mRow, mCol) = Symbol
    6   CONTINUE

C  convey the "screen" to the default output
        WRITE (*,*) nPack, ' point(s) per column'
        WRITE (*,*) 1./RowDel, ' standard deviations per row'
        DO 7  i = MxR2, 1, -1
    7     WRITE (*,*) (screen(i,j), j = 1, nCol + 2)

        RETURN
        END


        SUBROUTINE BasPlt (n, x, basis)
C       Draw a plot of some data and indicate a basis line on the
C       the plot.  That is, that line below which the data is
C       not meaningful.  The basis here is taken to be a constant.
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION x(1)
        CHARACTER*1 Blank, Symbol, hBordr, vBordr, BasSym
        PARAMETER (Blank = ' ', Symbol = 'O', BasSym = '=')
        PARAMETER (hBordr = '-', vBordr = '|')

        COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
        CHARACTER*1 screen(100, 150)

        IF (n .LT. 2) RETURN    ! not enough data

C  Find out how many points to pack per column and how many columns
        nPack = 1 + INT(FLOAT (n) / MaxCol - 1./n)
        nCol = INT((n - 1./n)/nPack + 1)

C  prepare the "screen" for drawing
        DO 1  j = 1, nCol + 2
          DO 1  i = 1, MxR2
            screen(i,j) = Blank
    1   CONTINUE
        DO 2  i = 2, nCol + 1
          screen(MxR2,i) = hBordr
    2     screen(1,i) = hBordr
        DO 3  i = 2, MaxRow + 1
          screen(i,nCol+2) = vBordr
    3     screen(i,1) = vBordr

C  get the data limits
        xMax = basis
        xMin = basis
        DO 4 i = 1, n
          IF (x(i) .GT. xMax) xMax = x(i)
          IF (x(i) .LT. xMin) xMin = x(i)
    4   CONTINUE
        RowDel = (MaxRow - 1) / (xMax - xMin)

C  show the basis line
        mPlus = 1 + INT((basis - xMin)*RowDel + 1)
        DO 5  i = 2, nCol + 1
    5     screen(mPlus,i) = basSym

C  draw the data (overdrawing the basis bars if necessary)
        DO 6  i = 1, n
          mCol = 1 + INT((i - 1./n)/nPack + 1)          ! addressing function
          mRow = 1 + INT((x(i) - xMin)*RowDel + 1)      ! +1 for the plot frame
          screen(mRow, mCol) = Symbol
    6   CONTINUE

C  convey the "screen" to the default output
        WRITE (*,*) nPack, ' point(s) per column'
        WRITE (*,*) 1./RowDel, ' units per row'
        DO 7  i = MxR2, 1, -1
    7     WRITE (*,*) (screen(i,j), j = 1, nCol + 2)

        RETURN
        END


        SUBROUTINE Plot (n,x,y)
C       Make a scatter plot on the default display device (UNIT=*).
C       MaxRow and MaxCol correspond to the display dimensions.
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION x(1), y(1)
        CHARACTER*1 Blank, Symbol, hBordr, vBordr
        PARAMETER (Blank = ' ', Symbol = 'O')
        PARAMETER (hBordr = '-', vBordr = '|')

        COMMON /space4/ screen, MaxCol, MaxRow, MxC2, MxR2
        CHARACTER*1 screen(100, 150)

        IF (n .LT. 2) RETURN    ! not enough data

C  prepare the "screen" for drawing
        DO 1  j = 1, MxC2
          DO 1  i = 1, MxR2
            screen(i,j) = Blank
    1   CONTINUE
        DO 2  i = 2, MaxCol+1
          screen(MxR2,i) = hBordr
    2     screen(1,i) = hBordr
        DO 3  i = 2, MaxRow+1
          screen(i,MxC2) = vBordr
    3     screen(i,1)  = vBordr

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
    4   CONTINUE
        ColDel = (MaxCol - 1) / (xMax - xMin)
        RowDel = (MaxRow - 1) / (yMax - yMin)

C  data scaling functions are offset by +1 for plot frame
        DO 5  i = 1, n
          mCol = 1 + INT((x(i) - xMin)*ColDel + 1)
          mRow = 1 + INT((y(i) - yMin)*RowDel + 1)
    5     screen(mRow, mCol) = Symbol

C  convey the "screen" to the default output
        WRITE (*,*) 1./ColDel, ' units per column'
        WRITE (*,*) 1./RowDel, ' units per row'
        DO 6  i = MaxRow + 2, 1, -1
    6     WRITE (*,*) (screen(i,j), j = 1, MaxCol + 2)
        RETURN
        END
