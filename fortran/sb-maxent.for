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


        END


        SUBROUTINE opus(n,npt,x,ox)     ! solution-space -> data-space
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        PARAMETER (MaxPts=300, MaxBin=102)
        COMMON /space1/ grid
        DIMENSION x(MaxBin), grid(MaxBin,MaxPts), ox(MaxPts)
C		ox = grid^T * x
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
C		x = grid * ox
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
          chisq = chisq + a**2			! get ChiSqr
   10     ox(j) = 2. * a / sigma(j)		! prepare to get cgrad (ChiSqr gradient)
        CALL tropus(n,npt,ox,cgrad)     ! cgrad = Grid * ox		SB1984: Eq. (8)
        test = zero     ! mismatch between entropy and ChiSquared gradients
        snorm = zero    ! norm (magnitude) of entropy vector: S(f)
        cnorm = zero    ! norm (magnitude) of ChiSqr vector: C(f)
        tnorm = zero    ! norm for the gradient term TEST
        fSum = zero     ! find the sum of the f-vector
        DO 12  i = 1, n
          fSum = fSum + f(i)							! volume fraction (sum of image vector)
          sgrad(i) = -LOG(f(i)/base(i)) / (blank*exp1)	! entropy gradient
          snorm = snorm + f(i) * sgrad(i)**2
          cnorm = cnorm + f(i) * cgrad(i)**2
          tnorm = tnorm + f(i) * sgrad(i) * cgrad(i)
   12   CONTINUE
        snorm = SQRT(snorm)			! SB1984: Eq. (22)
        cnorm = SQRT(cnorm)
        a = one
        b = one / cnorm
        IF (iter .GT. 1) THEN
          test = SQRT(0.5*(one-tnorm/(snorm*cnorm)))
          a = 0.5 / (snorm * test)
          b = 0.5 * b / test
        ENDIF
        DO 13 i = 1, n			! "xi" here is called "eta" in the SB1984 manuscript
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
        DO 17 k = 1, m		! prepare the search directions for the conjugate gradient technique
          s1(k) = zero		! S_mu
          c1(k) = zero		! C_mu
          DO 18 i = 1, n
            s1(k) = s1(k) + xi(i,k) * sgrad(i)
            c1(k) = c1(k) + xi(i,k) * cgrad(i)
   18     CONTINUE
          c1(k) = c1(k) / chisq
   17   CONTINUE
        DO 19 k = 1, m
          DO 19 l = 1, k
            s2(k,l) = zero		!  g_(mu,nu)
            c2(k,l) = zero		!  M_(mu,nu)
            DO 20 i = 1, n
              s2(k,l) = s2(k,l) - xi(i,k) * xi(i,l) / f(i)
   20       CONTINUE
            DO 21 j = 1, npt
              c2(k,l) = c2(k,l) + eta(j,k) * eta(j,l) / (sigma(j)**2)
   21       CONTINUE
            s2(k,l) = s2(k,l) / blank
            c2(k,l) = 2. * c2(k,l) / chisq
   19   CONTINUE
        c2(1,2) = c2(2,1)		! reflect on body-diagonal
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
C          CALL BasPlt (iter, Entropy, temp)
        END IF

  300   WRITE (*,*)
        WRITE (*,*) ' Residuals'
C        CALL ResPlt (npt, z)

        WRITE (*,*)
        WRITE (*,*) ' Distribution'
C        CALL BasPlt (n, f, blank)

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
        DO   k = 1, m
          DO   l = 1, m
            a(k,l) = bx * c2(k,l)  -  ax * s2(k,l)
          ENDDO
          b(k) = -(bx * c1(k)  -  ax * s1(k))
        ENDDO
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
