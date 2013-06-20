             test-fm : Project Name  (only 1st item is read)
    ../data/test.sas : SAS file, contains columns: Q  i  esd
     1e-08       100 : qMin qMax, 1/A  (1.0e-8 to 100 means all data)
                 100 : rhosq       : scattering contrast, 10^20 1/cm^-4
                   1 : fac         :   I = fac * ( i - bkg )
                   3 : err         : ESD = fac * err * esd
                 0.1 : bkg         :   I = fac * ( i - bkg )
                   1 : shapeModel  (1=spheroids, no others yet)
                   1 : Aspect Ratio
                   0 : Bin Type    (1=Lin, 0=Log)
                  40 : nRadii
        25      9000 : dMin dMax, A
                   1 : n, in N(D)*V^n, 0=N(D), 1=f(D), 2=i(D)
              1.0e-6 : defaultDistLevel  (MaxEnt only)
                  91 : IterMax
                   0 : slitLength, 1/A
              0.0002 : dLambda/Lambda
                   1 : method (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)
