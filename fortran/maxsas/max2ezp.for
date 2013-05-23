C	Max2Ezp.For
C	Read the output file from MaxSas.For
C	and extract the size distribution part, writing it
C	as input for the EASYPLOT program on the IPNS VAX.
C	The user must know how many data points there
C	are AND the file must have exactly "LinesSkipped" lines
C	of text before the first data point.

     	PROGRAM Max2Ezp
     	CHARACTER*50 InFile, OutFile
     	CHARACTER*150 aLine
     	PARAMETER ( MaxArr = 500 )
     	PARAMETER ( io = 11 )
     	PARAMETER (LinesSkipped = 12)
     	REAL X(MaxArr), Y(MaxArr)

     	WRITE (*, '(X,A,$)') 'What is the input file name? ==> '
     	READ (*, '(A)') InFile
     	WRITE (*, '(X,A,$)') 'What is the output file name? ==> '
     	READ (*, '(A)') OutFile
     	WRITE (*, '(X,A,$)') 'How many data points? ==> '
     	READ (*, *) NumPts
     	OPEN (io, FILE=InFile, READONLY, STATUS='OLD')
     	DO i = 1, LinesSkipped
     	  READ (io, '(A)'), aLine
     	END DO
     	DO i = 1, NumPts
     	  READ (io, FMT=*) D, f, xN, Sky
     	  X(i) = D
     	  Y(i) = f
     	END DO
  999	CLOSE (io)
     	OPEN (io, FILE=OutFile, STATUS='NEW')
     	WRITE (io, 1010) 
 1010	FORMAT (' ''O''  0.0 0.0')
     	DO i = 1, NumPts
     	  WRITE (io, *), X(i), Y(i)
     	END DO
     	WRITE (io, *) -999.0, -999.0
     	CLOSE (io)
     	END

