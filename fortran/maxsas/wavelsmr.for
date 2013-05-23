	SUBROUTINE Spheroid(r, nR, h, nH, G, v)
C	Particle form factor for spheroids of R x R x vR, R = radius
C	  ref: Roess and Shull; J Appl Phys 18 (1947) 308-313.
C	  coded: PRJ on 7.2.90
C	  modified: PRJ on 4.May.92
C	The method of eq.6 (second part) is implemented by direct evaluation
C	  using eq. (4) & (5) and numerical averaging.
C	Lookup tables are used to speed up the calculations.

C	h: array of scattering vectors numbered 1 to nH
C	r: array of radii numbered 1 to nR
C	G: grid of (form factors)**2 numbered 1 to nH & 1 to nR
C	v: aspect ratio
C	WaveMean: default (average) wavelength
C	WaveMin:  minimum wavelength for integration
C	WaveMax:  maximum wavelength for integration
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION h(1), r(1)
	PARAMETER (MaxPts=300, MaxBin=102)
	DIMENSION G(MaxBin,MaxPts)
C	PARAMETER (HalfPi = 1.570796326794896619)
C	PARAMETER (eta = 0.01)		!how close to a sphere?
	PARAMETER (nT = 500)		!# of table values
	DIMENSION t(nT), u(nT)		!table ordinates
	DIMENSION S(nT), Se(nT)		!table abcissae
	PARAMETER (nI = 100)		!# of integrands
	DIMENSION x(nI), y(nI)		!integrands
	Step(a,b,N) = EXP(LOG(b/a)/(FLOAT(N)-1.0)) + 1.0e-5	!step-size

C
C  Calculate the table of "spheres" form factors: S(u)
C
	uMin = h( 1) * r( 1) * (WaveMean/WaveMax) * MIN (v, 1.0)
	uMax = h(nH) * r(nR) * (WaveMean/WaveMin) * MAX (v, 1.0)
	du = Step (uMin, uMax, nU)	! the step size
	uNow = uMin
	DO  i = 1, nU
	  u(i) = uNow			! keep for the lookup table
	  S(i) = SphereNow (uNow)	! exact calculation at u = uNow
	  uNow = uNow * du		! for the next u
	END DO


C
C  Integrate to get the ellipsoidal form factor: Se(t)
C
	tMin = h( 1) * r( 1) * (WaveMean/WaveMax)
	tMax = h(nH) * r(nR) * (WaveMean/WaveMin)
	du = Step (tMin, tMax, nU)	! the step size
	tNow = tMin
	DO i = 1, nI
	  x(i) = (i - 1.0) / (nI - 1.0)	! ordinate of integration
	END DO
	DO  i = 1, nU
	  t(i) = tNow			! keep for the lookup table
	  DO j = 1, nI
	    uNow = tNow * SQRT (1 + x(j)**2 * (v**2 - 1))
	    CALL LookUp (uNow, u, nU, index)		!get the index
	    y(j) = ThisForm(uNow, u,S,nU, index)	!get the value
	  END DO
	  Se(i) = Area (x, y, nI)	! integration
	  tNow = tNow * dt		! for the next t
	END DO


C
C  Integrate to get the ellipsoidal form factor: Se(t)
C
	tMin = h( 1) * r( 1) * (WaveMean/WaveMax)
	tMax = h(nH) * r(nR) * (WaveMean/WaveMin)
	du = Step (tMin, tMax, nU)	! the step size
	tNow = tMin
	DO i = 1, nI
	  x(i) = (i - 1.0) / (nI - 1.0)	! ordinate of integration
	END DO
	DO  i = 1, nU
	  t(i) = tNow			! keep for the lookup table
	  DO j = 1, nI
	    uNow = tNow * SQRT (1 + x(j)**2 * (v**2 - 1))
	    CALL LookUp (uNow, u, nU, index)		!get the index
	    y(j) = ThisForm(uNow, u,S,nU, index)	!get the value
	  END DO
	  Se(i) = Area (x, y, nI)	! integration
	  tNow = tNow * dt		! for the next t
	END DO


C	Calculate the grid at each h(i) and r(j) from the modified E(*)
	DO 50  i = 1, nR		!at each scattering vector
	  DO 50  j = 1, nH		!at each radius
	    hr = r(i) * h(j)		!indexed into spheroids table
	    CALL LookUp (hr, u, nU, index)	!get the index
	    G(i,j) = ThisForm(hr, u,S,nU, index)	!get the value
   50	CONTINUE
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
	  sum = sum + (x(i) - x(i-1)) * (y(i) + y(i+1))
	END DO
	Area = 0.5 * sum
	RETURN
	END

