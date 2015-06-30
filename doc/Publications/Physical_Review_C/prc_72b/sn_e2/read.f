        program dill
        IMPLICIT NONE
         integer  n, i, a
         real*8 b, c
          n=0
 1         READ(5,*,end=2)
          n=n+1
          GOTO 1
 2        CONTINUE
          REWIND 5

         DO i=1,n
         read(5,*) a, b
         c = b* 23.596433203
           write(6,'(I4,2x,4E12.5)') a, b, c, c*5/10000.0, 
     :                b*5*((1.0031773*a**(1./6.))**4)/10000.0
         ENDDO

         end
