             reverse : ProjectName
         reverse.sas : SAS data file
   0.0005      0.025 : qMin  qMax, 1/A
                  10 : scattering contrast (10^20 cm^-4)
                   1 : scaling factor
                   1 : error scaling factor
                4.78 : background
                   1 : shape model
                   1 : aspect ratio, v, as in r * r * v*r
                   0 : bin type (1=Lin, 0=Log)
                 100 : number of bins
       80       8000 : Dmin & Dmax
                   1 : n, in N(D)*V^n
               1e-06 : default distribution level (MaxEnt)
                  50 : maximum number of MaxEnt iterations
                   0 : slit length, 1/A
              0.0002 : wavelength distribution, dLambda/Lambda
                   1 : analysisType (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)
