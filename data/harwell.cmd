             harwell : ProjectName
 harwell-bimodal.sas : SAS data file
    1e-08      0.098 : qMin  qMax, 1/A
                  10 : scattering contrast (10^20 cm^-4)
                   1 : scaling factor
                   1 : error scaling factor
              0.0001 : background
                   1 : shape model
                   1 : aspect ratio, v, as in r * r * v*r
                   0 : bin type (1=Lin, 0=Log)
                  40 : number of bins
       30        800 : Dmin & Dmax
                   1 : n, in N(D)*V^n
               1e-06 : default distribution level (MaxEnt)
                  32 : maximum number of MaxEnt iterations
                   0 : slit length, 1/A
              0.0002 : wavelength distribution, dLambda/Lambda
                   1 : analysisType (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)
