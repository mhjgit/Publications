      PROGRAM betaeq
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/densities/adens
      PARAMETER (e=15.8, s=32., gamma=0.6)
      PARAMETER (delta=0.2, rho0=0.16)

      a0=2.D0/3.D0/acos(-1.)/acos(-1.)
      DO i=1,1000
         rho=0.01*i
         x=rho/rho0
         adens=a0*(4.*s/197.327*(x**gamma))**3/rho
         CALL root(chi)         
         proton=(1.-chi)/2.
         erg=s*(x**gamma)*(chi**2)
         WRITE(6,'(E12.6,2x,e12.6,2x,e12.6)')rho,proton
      ENDDO

      END


      SUBROUTINE root(rtbis)
      IMPLICIT REAL*8 (a-h,o-z)      
      INTEGER jmax, j
      PARAMETER (jmax=1000)
      PARAMETER ( x1=0.5D0, x2=1.D0, xacc= 0.000001)
      REAL*8 dx,f,fmid,xmid

      fmid=beta_eq(x2)
      f=beta_eq(x1)
      IF (f*fmid.GE.0.d0) pause 'root must be bracketed in rtbis'
      IF (f.LT.0.d0) THEN
        rtbis=x1
        dx=x2-x1
      ELSE
        rtbis=x2
        dx=x1-x2
      ENDIF
      DO j=1,JMAX
         dx=dx*.5d0
         xmid=rtbis+dx
         fmid=beta_eq(xmid)
         IF(fmid.LE.0.d0) rtbis=xmid
         IF (DABS(dx).LT.xacc .OR. fmid.EQ.0.d0) RETURN
      ENDDO
      pause 'too many bisections in rtbis'

      END

      REAL FUNCTION beta_eq*8(x)
      IMPLICIT REAL*8 (a-h,o-z)      
      COMMON/densities/adens

      beta_eq=x**3-1./adens+x/adens

      END
