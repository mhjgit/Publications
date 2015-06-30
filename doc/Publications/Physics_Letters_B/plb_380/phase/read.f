
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(maxd=1000,mamax=30)
      DIMENSION xn(maxd),rhop(maxd),rhon(maxd),rhoe(maxd)
      DIMENSION yrhop(maxd),yrhon(maxd),yrhoe(maxd)
      DIMENSION radius(maxd), cden(maxd)


      OPEN(UNIT=5, FILE='electron.dat')
      OPEN(UNIT=7, FILE='out1')


      ndata=39
      nstar=105
      DO i=1,ndata
         READ (5,*) xn(i),rhoe(i)
      ENDDO

      DO i=1,nstar
         READ (7,*) cden(i), radius(i)
      ENDDO


c      CALL spline(xn,rhop,ndata,yp1,ypn,yrhop)
c      CALL spline(xn,rhon,ndata,yp1,ypn,yrhon)
      CALL spline(xn,rhoe,ndata,yp1,ypn,yrhoe)

      DO i=1,nstar
         x=cden(i)
         IF(x.gt.0.31) then
c         CALL splint(xn,rhop,yrhop,ndata,x,yp)
c         CALL splint(xn,rhon,yrhon,ndata,x,yn)
         CALL splint(xn,rhoe,yrhoe,ndata,x,ye)
         WRITE(6,1000) radius(i), ye
         endif
      ENDDO

 1000 FORMAT(5E12.4)

      END



      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=500)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END



      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
