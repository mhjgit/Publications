C
C
C
C     *************************************************************
C            Program to calculate the mass (SM), 
C            metric function and 
C            radius (r) for a neutron star, with a  
C            user supplied EOS. You should have     
C            energy density rho as function of the  
C            particle density n 
C            The dimensions are:                    
C            [rho] = MeV/fm^3  energy density given by 
C            user, like energy/baryon
C            [n]   = 1/fm^3    particle density                     
C            [P]   = MeV/fm^3  pressure                       
C            [sm]  = Kg        mass of neutron star                     
C            [r]   = m         radius of neutron star                     
C
C            Recall also that the constants defined 
C            below have the following dim           
C             [hbarc]=197.315 MeVfm
C            [g] = m^3Kg^-1s^-2 gravitational const 
C                = 6.67259*10^-11 m^3Kg^-1s^-2      
C                  or 6.67259 MeV^-2\hbarc/c^2
C            [sm_solar]=1.989D+30  kg
C            The solar mass is converted into MeV. Then
C            one has the relation
C                   1 kg=1.D+30/1.78266270D0  1 MeV/c^2
C            All taken from the Particle Data Group 
C            (1988)                                 
C            The program runs with dim-less variables
C            P=P/eps_c
C            rho=rho/eps_c
C            with eps_c the central energy-density of the star
C            r=r/r0 with r0=1.D+19 fm= 10km
C            m=m/m_solar
C            The dimension less constants used are
C            
C            const_1a=r0**3/sm_solar/speedoflight**2/eps_c
C            const_2=sm_solar*g/r0/speedoflight**2=const_2
C
C
C            The program uses runge-kutta till 4th order
C            The EOS is parametrized in routine EOSFIT
C            with routines taken from Numerical recepies
C
C
C            It is rather easy to plug in EOS for various
C            densities. You must just go carefully through
C            routine EOSFIT.
C
C
C            Coded by : M. Hjorth-Jensen            
C                       Dept. of Physics, Univ. of  
C                       Oslo, Norway, june/1991     
C                       based on white dwarf prog   
C                       whited.c of march 91                   
C            Last update : 13/feb/1996              
C            c++ version star.c is also available at
C            mhjensen@ect.unitn.it or 
C            http://www.ect.unitn.it/~mhjensen
C
C
C
C            Present version includes also the metric
C            function, though the calculation of
C            the moment of inertia is not yet 
C            implemented.
C
C            Observables given by print_output:
C          Prints out central density in fm^-3, radius in km,
C          neutron star mass in terms of solar masses, redshift
C          and reltivistic keplerian angular velocity in units of
C          1./(10^4 s^-1)= 10^-4 s
C
C
C     *************************************************************
C
C
C

      PROGRAM  nstar
      IMPLICIT REAL*8 (a-h,o-z)      
      PARAMETER(int_max=20000)
      PARAMETER( nync=30, n_diff_eq=3, maxd=1000)    
      COMMON/pressure/pp,x1,x2,tol
      COMMON/limits/xlim1,xlim2
      COMMON/parameters/const_1,const_2,cc,eps_c
      DIMENSION ynew(n_diff_eq), y(n_diff_eq), ders(n_diff_eq)
      COMMON/results/cen_dens(nync),radius(nync),smass(nync)

      EXTERNAL press
      LOGICAL log_div


C        *************************************
C             Read in parametrized EOS      *
C        *************************************

      CALL  eosfit                         ! read in user supplied EOS

C        *************************************
C        *        Set up of constants        *
C        *************************************

      r0=1.D+19                           ! 10^19 fm = 10 km,typical value of R
      xkg=1.D+30/1.78266270D0             ! 1 MeV/c^2=1.78266270D-30 kg
      hbarc=197.315D0                     ! \hbar  x c
      g=6.67259D-45*hbarc                 ! gravitational constant in fmMeV^-1
      fourpi=4.D0*DACOS(-1.D0)            ! 4*\pi 
      sm_solar=1.989D+30*xkg              ! solar mass in MeV/c^2
      delta=(3.D0-0.01D0)/DFLOAT(nync)    ! steps for central density
      c1=fourpi*(r0**3)/sm_solar          
      const_2=sm_solar*g/r0               ! dimless constant in TOV equation         
      xlim1=0.1096
      xlim2=0.31045
      x1=0.001
      x2=2.D0
      tol=1.D-8

C        ***************************************************
C          Calc mass, metric function  & radius of ns-star  
C        ***************************************************

      DO i=1,30                          ! loop over central densities
         ync=0.1D0+DFLOAT(i-1)*delta       ! central density in fm^-3
         pp=0.D0
         pp=press(ync)                     ! central pressure
         log_div=.FALSE.              
         eps_c=rho(ync,log_div)            ! central density in units MeV/fm^3
         const_1=c1*eps_c                  ! dimless const_1
         y(1)=pp/eps_c                     ! Dimensionless pressure
         y(2)=0.D0                         ! Dimensionless mass of the star
         x=0.001D0                         ! Dimensionless start radius
         j=0                               ! number of RK4 steps
         DO WHILE((pp.GT.0.D0).AND.(j.LE.int_max))
            j=j+1
            h=0.001                           ! Dimensionless step
            CALL derivatives(x,y,ders)        ! start values for derivatives
            CALL runge_kutta_4(x,y,ynew,ders,h) ! RK4 procedure
            x=x+h                             ! new dimless radius
            DO l=1,n_diff_eq                  ! y(1)= dimless pressure
               y(l)=ynew(l)                   ! y(2)= dimless mass
            ENDDO                             ! y(3)= metric function
            pp=y(1)*eps_c                     ! new pressure in MeV/fm^3
         ENDDO
         cen_dens(i)=ync                      ! central density
         radius(i)=x*10.                      ! radius in km
         smass(i)=y(2)                        ! star mass
      ENDDO          
      CALL print_output

      END

C
C     ************************************************************
C     *        4th-Runge-Kutta solution of coupled equations     *
C     *        See any textbook on numerical methods for details *
C     *        Taken from numerical recipies, Press et al.
C     ************************************************************
C

      SUBROUTINE runge_kutta_4(x,y,yout,dydx,h)
      IMPLICIT REAL*8 (a-h,o-z)      
      PARAMETER (n_diff_eq=3, nmax=50)
      DIMENSION y(n_diff_eq),dydx(n_diff_eq),yout(n_diff_eq),
     &          yt(nmax),dyt(nmax),dym(nmax)
      hh=h*0.5
      h6=h/6.
      xh=x+hh
      DO i=1,n_diff_eq
         yt(i)=y(i)+hh*dydx(i)
      ENDDO
      CALL derivatives(xh,yt,dyt)
      DO i=1,n_diff_eq
         yt(i)=y(i)+hh*dyt(i)
      ENDDO
      CALL derivatives(xh,yt,dym)
      DO i=1,n_diff_eq
         yt(i)=y(i)+h*dym(i)
         dym(i)=dyt(i)+dym(i)
      ENDDO
      CALL derivatives(x+h,yt,dyt)
      DO i=1,n_diff_eq
         yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      ENDDO

      END

C
C     ************************************************************
C                  dm/dr=ders(2)=const_1*r*r*energydensity  
C                  TOV equation=dp/dr= ders(1)
C                  metric function=ders(3)
C     ************************************************************
C

      SUBROUTINE derivatives(r,y,ders)
      IMPLICIT REAL*8 (a-h,o-z)      
      COMMON/pressure/pp,x1,x2,tol
      COMMON/parameters/const_1,const_2,cc,eps_c
      EXTERNAL press
      PARAMETER(n_diff_eq=3)
      DIMENSION y(n_diff_eq), ders(n_diff_eq)
      LOGICAL log_div

      log_div=.TRUE.
      pp=y(1)*eps_c                    ! pressure in MeV/fm^3
      IF(y(1).LE.0.D0) RETURN
      x=zbrent(press,x1,x2,tol)        ! the corresponding density
      ders(1)=-const_2*(rho(x,log_div)+y(1))*(y(2)+
     &         const_1*(r**3)*y(1))
     &         /(r*r-2*const_2*r*y(2))
      ders(2)=const_1*(r**2)*rho(x,log_div)
      ders(3)=2*const_2*(y(2)+const_1*y(1))/(r*r-2*const_2*r*y(2))
      END


C
C     ************************************************************
C               Energy in units of MeV/fm^3 if log_div=.FALSE.
C               dimless if log_div =.TRUE., divide by central
C               energy density eps_c
C     ************************************************************
C

      REAL FUNCTION rho*8(x,log_div)
      IMPLICIT REAL*8 (a-h,o-z)     
      COMMON/parameters/const_1,const_2,cc,eps_c 
      PARAMETER(mamax=30)
      COMMON/limits/xlim1,xlim2
      COMMON/set1/e1(mamax),ma1

      LOGICAL log_div
      rho=0.
      IF(x.LE.xlim1) THEN
        rho=x*940.6 + 247.5*(x**2.48)
      ELSEIF(x.GT.xlim1) THEN
        DO i=1,ma1
           rho=rho+e1(i)*(x**(DFLOAT(i-1)/3.d0))*197.315
        ENDDO
      ENDIF
      IF(log_div) THEN
         rho=rho/eps_c
      ENDIF
      END



C
C     ************************************************************
C               Pressure in units of MeV/fm^3 
C     ************************************************************
C

      REAL FUNCTION press*8(x)
      IMPLICIT REAL*8 (a-h,o-z)      
      PARAMETER(mamax=30)
      COMMON/set1/e1(mamax),ma1
      COMMON/limits/xlim1,xlim2
      COMMON/pressure/pp,x1,x2,tol      
      
      press=-pp
      IF(x.LE.xlim1) THEN
         press=press+366.3*(x**2.48)
      ELSEIF(x.GT.xlim1)THEN
         DO i=1,ma1
            press=press+(DFLOAT(i-1)/3.D0*e1(i)*
     &           (x**(DFLOAT(i-1)/3.D0))
     &           -e1(i)*(x**(DFLOAT(i-1)/3.d0)))*197.315
         ENDDO
      ENDIF

      END



C
C     ************************************************************
C          Prints out central density in fm^-3, radius in km,
C          neutron star mass in terms of solar masses, redshift
C          and reltivistic keplerian angular velocity in units of
C          1./(10^4 s^-1)= 10^-4 s
C     ************************************************************
C



      SUBROUTINE  print_output
      IMPLICIT REAL*8 (a-h,o-z)      
      PARAMETER( nync=30)    
      COMMON/results/cen_dens(nync),radius(nync),smass(nync)

      WRITE(6,1000)      
      DO i=1,30                          ! loop over central densities
         redshift=1./DSQRT(1.D0-2*smass(i)/radius(i))-1.D0
         rel_kep_vel=0.72*DSQRT(1.D+3*smass(i)/(radius(i)**3))
         WRITE(6,1001) cen_dens(i),radius(i),smass(i),redshift,
     &              rel_kep_vel
      ENDDO          
 1000 FORMAT('density in fm^3:', 'Radius in km:', 'Mass:','Redshift:',
     &        'Relativistic Keplerian angular velocity:')
 1001 FORMAT(5D12.3)                        
      
      END



C
C     ************************************************************
C     *      Reads in the parametrization of the energy density  
C     *      and returns the energy density as function of xn    
C     *      as well as the pressure as function of xn
C     ************************************************************
C

      SUBROUTINE  eosfit
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(maxd=1000,mamax=30)
      DIMENSION xn(maxd),p(maxd)
      COMMON/set1/e1(mamax),ma1
      DIMENSION e(maxd),sig(maxd)
      DIMENSION cvm(mamax,mamax)
      DIMENSION v(mamax,mamax),u(maxd,mamax),w(mamax)
      EXTERNAL func_e

      OPEN(UNIT=5, FILE='super.dat')


      READ (5,*) ndata, ma1                    ! read in data-points and
      IF (ndata.GT.maxd) THEN                 ! max polynomial x^ma
	 WRITE (6,*) 'bad input: maxd'
	 STOP
      ENDIF
      IF (ma1.GT.mamax) THEN
         WRITE (6,*) 'bad input: ma1'
         STOP
      ENDIF
      DO i=1,ndata
         READ (5,*) xn(i),e(i),xy,zy                ! read in energy as func of
         sig(I)=1.D0                          
      ENDDO

      CALL svdfit(xn,e,sig,ndata,e1,ma1,u,v,w,maxd,mamax,chisq,func_e)
      CALL svdvar(v,ma1,mamax,w,cvm,mamax)

      WRITE(6,*) ' Quality of energy fit'
      WRITE(6,*)  ma1
      DO i=1,ma1
         WRITE (6,*) e1(i)
      ENDDO
      WRITE (6,'(''CHISQ = '',F12.5)') chisq
      IF(chsq.GT.0.001D0) THEN
        WRITE(6,*) 'chi-square tolerance too high, try again'
        STOP
      ENDIF
      DO j=1,ndata
         sum=0.d0                              
         sump=0.d0
         DO i=1,ma1
            sum=sum+e1(i)*(xn(j)**(DFLOAT(i-1)/3.d0))
            sump=sump+DFLOAT(i-1)/3.D0*e1(i)*
     &           (xn(j)**(DFLOAT(i-1)/3.D0))
     &           -e1(i)*(xn(j)**(DFLOAT(i-1)/3.d0))
         ENDDO
         p(j)=sump
         WRITE(6,'(5d12.4)') xn(j),sum*197.315, sump*197.315
      ENDDO


      END


      SUBROUTINE func_e(x,afunc,ma)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION afunc(1)
      afunc(1)=1.D0
      DO I=2,ma
         afunc(i)=afunc(i-1)*x**(1./3.)
      ENDDO
      END



C     *   programs from numerical recepies     


      SUBROUTINE SVDFIT(X,Y,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ,FUNCS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=1000,MMAX=50,TOL=1.d-10)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *    U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)
      DO 12 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        TMP=1.d0/SIG(I)
        DO 11 J=1,MA
           U(I,J)=AFUNC(J)*TMP
11      CONTINUE
        B(I)=Y(I)*TMP
12    CONTINUE
      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.d0
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.d0
14    CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.d0
      DO 16 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        SUM=0.d0
        DO 15 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
15      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
16    CONTINUE
      END

      SUBROUTINE SVDVAR(V,MA,NP,W,CVM,NCVM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MMAX=30)
      DIMENSION V(NP,NP),W(NP),CVM(NCVM,NCVM),WTI(MMAX)
      DO 11 I=1,MA
        WTI(I)=0.d0
        IF(W(I).NE.0.d0) WTI(I)=1.d0/(W(I)*W(I))
11    CONTINUE
      DO 14 I=1,MA
        DO 13 J=1,I
          SUM=0.d0
          DO 12 K=1,MA
            SUM=SUM+V(I,K)*V(J,K)*WTI(K)
12        CONTINUE
          CVM(I,J)=SUM
          CVM(J,I)=SUM
13      CONTINUE
14    CONTINUE
      END

      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0d0
      SCALE=0.0d0
      ANORM=0.0d0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+dABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0d0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0d0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(dABS(W(I))+dABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0d0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0d0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0d0
            V(J,I)=0.0d0
31        CONTINUE
        ENDIF
        V(I,I)=1.0d0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0d0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0d0) THEN
          G=1.0d0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0d0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0d0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0d0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((dABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((dABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0d0
          S=1.0d0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((dABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0d0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0d0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0d0*H*Y)
          G=SQRT(F*F+1.0d0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0d0
          S=1.0d0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0d0) THEN
              Z=1.0d0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0d0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      END

      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.d0
        IF(W(J).NE.0.d0)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.d0
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      END





      FUNCTION ZBRENT(FUNC,X1,X2,TOL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.D-8)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
c      IF(FB*FA.GT.0.) PAUSE 'Root must be bracketed for ZBRENT.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(DABS(FC).LT.DABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*DABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(DABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(DABS(E).GE.TOL1 .AND. DABS(FA).GT.DABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=DABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-DABS(TOL1*Q),DABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(DABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+DSIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      PAUSE 'ZBRENT exceeding maximum iterations.'
      ZBRENT=B
      RETURN
      END







