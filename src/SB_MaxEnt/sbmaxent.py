#!/usr/bin/env python

'''
sbmaxent

Entropy maximization routine as described in the article
J Skilling and RK Bryan; Mon Not R Astr Soc 211 (1984) 111 - 124.
("Mon Not R Astr Soc" is the 
"Monthly Notices of the Royal Astronomical Society")

This program is a translation from the C which is from the FORTRAN 
which is from the BASIC version.

:license: Copyright (c) 2013, UChicago Argonne, LLC
:license: This file is distributed subject to a Software License Agreement found
     in the file LICENSE that is included with this distribution. 
'''

########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################


import os
import sys
import math
import numpy
import rst_table


TEST_LIMIT        = 0.05                    # for convergence
CHI_SQR_LIMIT     = 0.01                    # maximum difference in ChiSqr for a solution
SEARCH_DIRECTIONS = 3                       # <10.  This code requires value = 3
SD1               = SEARCH_DIRECTIONS+1     # for convenience below
RESET_STRAYS      = 1                       # was 0.001, correction of stray negative values
DISTANCE_LIMIT_FACTOR = 0.1                 # limitation on df to constrain runaways

MAX_MOVE_LOOPS    = 500                     # for no solution in routine: move, 
MOVE_PASSES       = 0.001                   # convergence test in routine: move


chisq = 0.0
chizer = 0.0
beta = None
c1 = None
c2 = None
s1 = None
s2 = None


class MaxEntException(Exception): 
    '''Any exception from this module'''
    pass


class MaxEnt(object):
    '''manage the iterative method described by Skilling and Bryan'''
    # TODO: Consider making this a class, due to global variable use 
    
    def __init__(datum, sigma, base, IterMax, image_to_data, data_to_image, G):
        pass


def opus (image, G):
    '''
    opus: transform solution-space -> data-space:  [G]^tr * image
    
    default definition, caller can use this definition or provide an alternative
    
    :param float[N] image: solution
    :param float[M][N] G: transformation matrix
    :returns float[M]: calculated data
    '''
    return G.transpose().dot(image)


def tropus (data, G):
    '''
    tropus: transform data-space -> solution-space:  [G] * data
    
    default definition, caller can use this definition or provide an alternative
    
    :param float[M] data: observations
    :param float[M][N] G: transformation matrix
    :returns float[N]: calculated image
    '''
    return G.dot(data)


def MaxEnt_SB(datum, sigma, base, IterMax, image_to_data, data_to_image, G, report=True):
    '''
    do the complete Maximum Entropy algorithm of Skilling and Bryan
    
    :param float datum[]:
    :param float sigma[]:
    :param float base[]:
    :param int IterMax:
    :param obj image_to_data: opus function 
    :param obj data_to_image: tropus function
    :param float[][] G: transformation matrix
    
    :returns float[]: f (:math:`f(r) dr`)
    '''
    global beta
    
    n   = len(base)
    npt = len(datum)
    ox      = 0*numpy.ndarray((npt))
    z       = 0*numpy.ndarray((npt))
    cgrad   = 0*numpy.ndarray((n))
    sgrad   = 0*numpy.ndarray((n))
    # Note that the order of subscripts for
    # "xi" and "eta" has been reversed from
    # the convention used in the FORTRAN version
    # to enable parts of them to be passed as
    # as vectors to "image_to_data" and "data_to_image".
    xi      = 0*numpy.ndarray((SEARCH_DIRECTIONS, n))
    eta     = 0*numpy.ndarray((SEARCH_DIRECTIONS, npt))
    beta    = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    s1      = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    s2      = 0*numpy.ndarray((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))
    c1      = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    c2      = 0*numpy.ndarray((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))

    exp1 = math.exp (1.0)
    blank = sum(base) / len(base)   # use the average value of base

    chizer, chtarg = npt*1.0, npt*1.0
    f = base * 1.0                              # starting distribution is base

    fSum  = sum(f)                              # find the sum of the f-vector
    z = (datum - image_to_data (f, G)) / sigma  # standardized residuals, SB eq. 3
    chisq = sum(z*z)                            # Chi^2, SB eq. 4

    for iter in range(IterMax):
        ox = -2 * z / sigma                        # gradient of Chi^2

        cgrad = data_to_image (ox, G)              # cgrad[i] = del(C)/del(f[i]), SB eq. 8
        sgrad = -numpy.log(f/base) / (blank*exp1)  # sgrad[i] = del(S)/del(f[i])
        snorm = math.sqrt(sum(f * sgrad*sgrad))    # entropy term
        cnorm = math.sqrt(sum(f * cgrad*cgrad))    # ChiSqr term 
        tnorm = sum(f * sgrad * cgrad)             # norm for gradient term TEST 

        a = 1.0
        b = 1.0 / cnorm
        if iter == 0:
            test = 0.0     # mismatch between entropy and ChiSquared gradients
        else:
            test = math.sqrt( ( 1.0 - tnorm/(snorm*cnorm) )/2 ) # SB eq. 37?
            a = 0.5 / (snorm * test)
            b *= 0.5 / test
        xi[0] = f * cgrad / cnorm
        xi[1] = f * (a * sgrad - b * cgrad)

        eta[0] = image_to_data (xi[0], G);          # image --> data
        eta[1] = image_to_data (xi[1], G);          # image --> data
        ox = eta[1] / (sigma * sigma)
        xi[2] = data_to_image (ox, G);              # data --> image
        a = 1.0 / math.sqrt(sum(f * xi[2]*xi[2]))
        xi[2] = f * xi[2] * a
        eta[2] = image_to_data (xi[2], G)           # image --> data
        
#         print_arr("MaxEnt: eta.transpose()", eta.transpose())
#         print_arr("MaxEnt: xi.transpose()", xi.transpose())

        c1 = xi.dot(cgrad) / chisq
        s1 = xi.dot(sgrad)
#         print_vec("MaxEnt: c1", c1)
#         print_vec("MaxEnt: s1", s1)

        for k in range(SEARCH_DIRECTIONS):
            for l in range(k+1):
                c2[k][l] = 2 * sum(eta[k] * eta[l] / sigma/sigma) / chisq
                s2[k][l] = -sum(xi[k] * xi[l] / f) / blank
#         print_arr("MaxEnt: c2", c2)
#         print_arr("MaxEnt: s2", s2)

        # reflect across the body diagonal
        for k, l in ((0,1), (0,2), (1,2)):
            c2[k][l] = c2[l][k]
            s2[k][l] = s2[l][k]
 
        beta[0] = -0.5 * c1[0] / c2[0][0]
        beta[1] = 0.0
        beta[2] = 0.0
        if (iter > 0):
            MaxEntMove(fSum, blank, chisq, chizer, c1, c2, s1, s2)

        f_old = f.copy()    # preserve the last image
        f += xi.transpose().dot(beta)   # move the image towards the solution
        
        # As mentioned at the top of p.119,
        # need to protect against stray negative values.
        # In this case, set them to RESET_STRAYS * base[i]
        #f = f.clip(RESET_STRAYS * blank, f.max())
        for i in range(n):
            if f[i] <= 0.0:
                f[i] = RESET_STRAYS * blank
        df = f - f_old
        fSum = sum(f)
        fChange = sum(df)

        # calculate the normalized entropy
        S = -sum((f/fSum) * numpy.log(f/fSum))      # normalized entropy, S&B eq. 1
        z = (datum - image_to_data (f, G)) / sigma  # standardized residuals
        chisq = sum(z*z)                            # report this ChiSq

        if report:
            print "%3d/%3d" % ((iter+1), IterMax)
            print " %5.2lf%% %8lg" % (100*test, S)
            if iter > 0:
                value = 100*( math.sqrt(chisq/chtarg)-1)
            else:
                value = 0
            print " %12.5lg %10.4lf" % ( math.sqrt (chtarg/npt), value )
            print "%12.6lg %8.2lf\n" % (fSum, 100*fChange/fSum)

        # See if we have finished our task.
        # do the hardest test first
        if (abs(chisq/chizer-1.0) < CHI_SQR_LIMIT) and  (test < TEST_LIMIT):
            return f     # solution FOUND returns here
    
    return None              # no solution after IterMax iterations


def MaxEntMove(fSum, blank, chisq, chizer, c1, c2, s1, s2):
    '''
    move beta one step closer towards the solution
    '''
    global beta
    
    # TODO: integrate NumPy
    a_lower, a_upper = 0., 1.          # bracket  "a"
    cmin = ChiNow (a_lower, c1, c2, s1, s2)
    #print "MaxEntMove: cmin = %g" % cmin
    if cmin*chisq > chizer:
        ctarg = (1.0 + cmin)/2
    else:
        ctarg = chizer/chisq
    f_lower = cmin - ctarg
    f_upper = ChiNow (a_upper, c1, c2, s1, s2) - ctarg

    fx = 2*MOVE_PASSES      # just to start off
    loop = 1
    while abs(fx) >= MOVE_PASSES and loop <= MAX_MOVE_LOOPS:
        a_new = (a_lower + a_upper) * 0.5           # search by bisection
        fx = ChiNow (a_new, c1, c2, s1, s2) - ctarg
        # tighten the search range for the next pass
        if f_lower*fx > 0:
            a_lower, f_lower = a_new, fx
        if f_upper*fx > 0:
            a_upper, f_upper = a_new, fx
        loop += 1

    if abs(fx) >= MOVE_PASSES or loop > MAX_MOVE_LOOPS:
        msg = "MaxEntMove: Loop counter = " 
        msg += str(MAX_MOVE_LOOPS)
        msg += '  No convergence in alpha chop'
        raise MaxEntException(msg)

    w = Dist (s2, beta);
    m = SEARCH_DIRECTIONS
    if (w > DISTANCE_LIMIT_FACTOR*fSum/blank):        # invoke the distance penalty, SB eq. 17
        for k in range(m):
            beta[k] *= math.sqrt (fSum/(blank*w))
    chtarg = ctarg * chisq
    return w, chtarg, loop, a_new, fx, beta


def Dist(s2, beta):
    '''measure the distance of this possible solution'''
    w = 0
    n = beta.shape[0]
    for k in range(n):
        z = -sum(s2[k] * beta)
        w += beta[k] * z
    return w


def ChiNow(ax, c1, c2, s1, s2):
    '''ChiNow'''
    global beta
    
    bx = 1 - ax
    a =   bx * c2  -  ax * s2
    b = -(bx * c1  -  ax * s1)

    beta = ChoSol(a, b)
    w = 1.0
    for k in range(SEARCH_DIRECTIONS):
        z = sum(c2[k] * beta)
        w += beta[k] * (c1[k] + z/2)
    return w


def ChoSol(a, b):
    '''
    ChoSol: ? chop the solution vectors ?
    
    :returns: new vector beta
    '''
    n = b.shape[0]
    fl = numpy.ndarray((n, n))*0
    bl = numpy.ndarray((n))*0
    
    #print_arr("ChoSol: a", a)
    #print_vec("ChoSol: b", b)

    if (a[0][0] <= 0):
        msg = "ChoSol: a[0][0] = " 
        msg += str(a[0][0])
        msg += '  Value must be positive'
        raise MaxEntException(msg)

    # first, compute fl from a
    # note fl is a lower triangular matrix
    fl[0][0] = math.sqrt (a[0][0])
    for i in (1, 2):
        fl[i][0] = a[i][0] / fl[0][0]
        for j in range(1, i+1):
            z = 0.0
            for k in range(j):
                z += fl[i][k] * fl[j][k]
                #print "ChoSol: %d %d %d  z = %lg" % ( i, j, k, z)
            z = a[i][j] - z
            if j == i:
                y = math.sqrt(z)
            else:
                y = z / fl[j][j]
            fl[i][j] = y
    #print_arr("ChoSol: fl", fl)

    # next, compute bl from fl and b
    bl[0] = b[0] / fl[0][0]
    for i in (1, 2):
        z = 0.0
        for k in range(i):
            z += fl[i][k] * bl[k]
            #print "\t", i, k, z
        bl[i] = (b[i] - z) / fl[i][i]
    #print_vec("ChoSol: bl", bl)

    # last, compute beta from bl and fl
    beta = numpy.ndarray((n))
    beta[-1] = bl[-1] / fl[-1][-1]
    for i in (1, 0):
        z = 0.0
        for k in range(i+1, n):
            z += fl[k][i] * beta[k]
            #print "\t\t", i, k, 'z=', z
        beta[i] = (bl[i] - z) / fl[i][i]
    #print_vec("ChoSol: beta", beta)

    return beta


def G_term(q, r, V, rhosq):
    '''Calculates the response matrix :math:`G(Q,r)` 
    
    :param double q: :math:`Q`
    :param double r: :math:`r`
    :param double rhosq: :math:`|\\Delta\\rho|^2`, the scattering contrast
    :returns double: G(Q,r)
    '''
    qr = q*r
    F = 3 * (math.sin(qr) - qr*math.cos(qr)) / (qr*qr*qr)
    shape = F*F
    value = rhosq*1e20 * V * shape
    return value


def print_vec(text, a):
    '''print the contents of a vector to the console'''
    n = a.shape[0]
    print "%s[ = (" % text,
    for i in range(n):
        s = " %g, " % a[i]
        print s,
    print ")"


def print_arr(text, a):
    '''print the contents of an array to the console'''
    n, m = a.shape
    print "%s[][] = (" % text
    for i in range(n):
        print " (",
        for j in range(m):
            print " %g, " % a[i][j],
        print "),"
    print ")"


def readTextData(filename):
    '''return q, I, dI from a 3-column text file'''
    if not os.path.exists(filename):
        raise Exception("file not found: " + test_data_file)
    buf = [line.split() for line in open(filename, 'r').readlines()]
    M = len(buf)
    buf = zip(*buf)         # transpose rows and columns
    q  = numpy.array(buf[0], dtype=numpy.float64)
    I  = numpy.array(buf[1], dtype=numpy.float64)
    dI = numpy.array(buf[2], dtype=numpy.float64)
    return q, I, dI


def test_opus_tropus():
    '''test routine for Dist()
    
    test values are trivial checks
    '''
    b = numpy.array((1,2,3))
    R = numpy.array(((1,0,0),(0,1,0),(0,0,1)))      # identity matrix
    print "tropus:", tropus(b, R), b
    print "opus:", opus(b, R), b
    m = numpy.array( ((1,1,0), (0,1,1)) )
    print "tropus:", tropus(b, m), [3, 5]


def test_Dist():
    '''test routine for Dist()
    
    test values obtained from C version of sizes program
    '''
    beta    = numpy.array((-0.0287502, 0.0705577, -0.0257728))
    s2      = numpy.array(( (-1e+06, 177480, -248889), 
                            (177480, -1e+06, -641283),  
                            (-248889, -641283, -1e+06) ))
    print "Dist:", Dist(s2, beta), 5225.79


def test_ChoSol():
    '''test routine for ChoSol()
    
    test values obtained from C version of sizes program
    '''
    # test ChoSol()
    a = numpy.array(( ( 89.2816, 16.6402, 55.369,),
                      ( 16.6402, 42.875, 66.8581,),
                      ( 55.369, 66.8581, 124.94,),
                     ))
    b = numpy.array((-9.72508, 1.72601, -2.42046))
    print "ChoSol:", ChoSol(a, b), ( -0.076233, 0.286152,  -0.138715)

    a = numpy.array(( ( 1e+06, -177480, 248889,),
                      ( -177480, 1e+06, 641283,),
                      ( 248889, 641283, 1e+06,),
                     ))
    b = numpy.array((165268, 31304, 84048.4))
    print "ChoSol:", ChoSol(a, b), ( 0.17638, 0.0626079,  6.42578e-17)


def test_ChiNow():
    '''test routine for ChiNow()
    
    test values obtained from C version of sizes program
    '''
    ax = 0.0
    c1 = numpy.array( ( 9.72508,  -1.72601,  2.42046, ) )
    c2 = numpy.array( (
      ( 89.2816,  16.6402,  55.369, ),
      ( 16.6402,  42.875,  66.8581, ),
      ( 55.369,  66.8581,  124.94, ),
    ))
    s1 = numpy.array( ( 165268,  31304,  84048.4, ) )
    s2 = numpy.array( (
      ( -1e+06,  177480,  -248889, ),
      ( 177480,  -1e+06,  -641283, ),
      ( -248889,  -641283,  -1e+06, ),
    ))
    # ChiNow: result=0.214487
    print "ChiNow:", ChiNow(ax, c1, c2, s1, s2), 0.214487


def test_MaxEntMove():
    '''test routine for MaxEntMove()
    
    test values obtained from C version of sizes program
    '''
    chisq = 4038.84
    chizer = 92.0
    fSum = 0.00768981
    blank = 1e-6
    c1 = numpy.array(( 9.72508,  -1.72601,  2.42046, ))
    c2 = numpy.array( (
          ( 89.2816,  16.6402,  55.369, ),
          ( 16.6402,  42.875,  66.8581, ),
          ( 55.369,  66.8581,  124.94, ),
        ))
    s1 = numpy.array( ( 165268,  31304,  84048.4, ) )
    s2 = numpy.array((
          ( -1e+06,  177480,  -248889, ),
          ( 177480,  -1e+06,  -641283, ),
          ( -248889,  -641283,  -1e+06, ),
        ) )
    w, chtarg, loop, a_new, fx, beta = MaxEntMove(fSum, blank, chisq, chizer, c1, c2, s1, s2)
    print "MaxEntMove: w", w, 5225.79
    print "MaxEntMove: chtarg", chtarg, 2452.56
    print "MaxEntMove: loop", loop, 482
    print "MaxEntMove: a_new", a_new, 3.24249e-05
    print "MaxEntMove: fx", fx,  -0.000196529
    print "MaxEntMove: beta", beta, (-0.0348757,  0.0855907,  -0.0312639)


def test_MaxEnt_SB(report=True):
    print "MaxEnt_SB: "
    test_data_file = os.path.join('..', '..', 'data', 'test.sas')
    rhosq = 100     # scattering contrast, 10^20 1/cm^-4
    bkg   = 0.1     #   I = I - bkg
    dMin, dMax, nRadii = 25, 9000, 40
    defaultDistLevel = 1.0e-6
    IterMax = 40
    errFac = 1.05
    
    r    = numpy.exp(numpy.linspace(math.log(dMin), math.log(dMax), nRadii))/2
    dr   = r * (r[1]/r[0] - 1)          # step size
    f_dr = numpy.ndarray((nRadii)) * 0  # volume fraction histogram
    b    = numpy.ndarray((nRadii)) * 0 + defaultDistLevel  # MaxEnt "sky background"
    
    qVec, I, dI = readTextData(test_data_file)
    G = numpy.ndarray((nRadii, len(qVec)))
    V = (4./3.) * math.pi * r*r*r * 1e-24
    for i, rVal in enumerate(r):
        for j, qVal in enumerate(qVec):
            G[i][j] = G_term(qVal, rVal, V[i], rhosq)
    
    f_dr = MaxEnt_SB(I - bkg, dI*errFac, b, IterMax, opus, tropus, G, report=report)
    if f_dr is None:
        print "no solution"
        return
    
    print "solution reached"
    t = rst_table.Table()
    t.labels = ('r', 'dr', 'f(r) dr', )
    t.rows.extend(zip(*(r.tolist(), dr.tolist(), f_dr.tolist())))
    print t.reST()


def tests():
    test_opus_tropus()
    test_Dist()
    test_ChoSol()
    test_ChiNow()
    test_MaxEntMove()
    test_MaxEnt_SB(report=False)


if __name__ == '__main__':
    tests()
