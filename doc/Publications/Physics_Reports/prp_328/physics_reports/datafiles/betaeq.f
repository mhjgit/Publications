From lars.engvik@fys.uio.no Wed Sep  9 18:55:01 1998
Return-Path: <lars.engvik@fys.uio.no>
Delivery-Date: Wed, 9 Sep 1998 18:55:00 +0200
Received: from isak.uio.no (actually isak.uio.no [129.240.84.99]) by pat.uio.no 
          with SMTP (PP); Wed, 9 Sep 1998 18:54:45 +0200
Received: by isak.uio.no ; Wed, 9 Sep 1998 18:54:44 +0200 (METDST)
Date: Wed, 9 Sep 1998 18:54:44 +0200 (METDST)
From: lars.engvik@fys.uio.no
Message-Id: <199809091654.SAA19520@isak.uio.no>
To: m.h.jensen@fys.uio.no
Status: R
Content-Length: 15426

Hei morten,
nedenfor har de programmet som
regner ut betalikevekt ved hjelp av energi per partikkel
i nøytronmaterie og symmetrisk materie.
Det er  sparsomt med kommentarer i programmet.
Subroutinen eos_input leser inn disse fra fil. 
I tillegg setter den opp punkter på ny grid.
Tilstandslikningen som blir brukt til betalikevekt blir lagret i
common/eos/.
Her kan du mate inn de ønskede verdier direkte uavhengig av eos_input.

Programmet regner i tillegg ut trykk. Jeg er ikke helt trygg på
den numeriske nøyaktigheten for trykket.  Det er nok følsomt for 
hvilken tetthets grid  tilstandlikningen er beregnet på.

 
input filen ser slik ut
21
.0346    .0492    .0675    .0899    .1167    
.1484    .1854    .2280    .2767    .3319    
.3939    .4633    .5404    .6256    .7192    
.8218    .9338   1.0554   1.3295   1.4828   1.6474
1.0
  5.58     6.80     8.17     9.83    11.99    
 15.07    19.64    26.55    36.77    51.73    
 72.25    99.54   133.64   175.37   224.32   
281.02   344.80   417.02   618.23   762.07   894.63
0.0
  8.49    -8.46    -9.68   -10.99   -12.21   
-12.97   -12.71   -10.95    -6.78     1.05    
 13.36    31.72    56.19    88.28   127.31   
173.75   227.88   288.36   441.34   550.26   670.48


      subroutine eos_input
      IMPLICIT REAL*8 (A-H,O-Z)
      common/eos/sym_e(101),e_nuc(101),dens(101),ndens
      PARAMETER(maxd=1000)
      PARAMETER (ldi=101)
      dimension e_neu1(maxd),e_nuc1(maxd)
      dimension e_neu(101)
      dimension dens1(maxd)
      open(5,file='ipt.dat')
      read(5,*)ndens1
      read(5,*)(dens1(i),i=1,ndens1)
      ndens=80
      xxd=0.04d0
      xd=0.02d0
      do i=1,ndens
         dens(i)=xxd
         xxd=xxd+xd
      enddo
      read(5,*)frac
      read(5,*)(e_neu1(j),j=1,ndens1)
      read(5,*)frac
      read(5,*)(e_nuc1(j),j=1,ndens1)
      close(5)
      do j=1,ndens
         xd=dens(j)
         npol=4
         if(xd.le.dens1(1).or.xd.ge.dens1(ndens1-1))npol=3 
         call interp1d(dens1,e_nuc1,ndens1,npol,xd,xb)
         e_nuc(j)=xb
         call interp1d(dens1,e_neu1,ndens1,npol,xd,xb)
         e_neu(j)=xb
         sym_e(j)=e_neu(j)-e_nuc(j)
      enddo
 122  format(3f12.3)
      open(7,file='sym.dat')
      do i=1,ndens
         write(7,122)dens(i),sym_e(i)
      enddo
      close(7)
      close(17)
 101  format(10f12.4)
      return
      end


      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON/MEV/HC,HC3,wn 
      common/leptons/wmu,wel
      common/eos/sym_e(101),e_nuc(101),dens1(101),ndens
      DIMENSION CHEMP(100),CHEMN(100)
      DIMENSION YNEU(100),CHEM(100),EN(100)
     :         ,EM(100),EE(100),ETOT(100)
     :         ,y_el(100),dmass(100),dnum(100),edens(100) 
      dimension x1(1),fx1(1),sqq(120),sau(120)
      character*9 cc
      data pi/3.141592653589793d0/
      open(6,file='out')
      open(8,file='edens.dat')
      CALL eos_input
      PI2=PI*PI
      NRHO=ndens
      IMAX=800    
      ONE=1.D0              
      HALF=0.5D0     
      ZER=0.d0
      HC=197.328D0
      HC3=HC*HC*HC
      WN=938.926D0     
      wnkg=1.6738015d-27
      WMU=105.7D0
      wel=0.511d0
      WMU2=WMU*WMU
      xmevg=1000.d0*wnkg/wn
      xfmcm3=1.d39  
      dgcm=xmevg*xfmcm3/1.d14  !mass density given 10^-14 g/cm^3 
      X=1.D0                       ! x ->neutron fraction
      write(6,*)
      write(6,101)
      DO IB=1,ndens
         dens=dens1(ib)
         baryon=dens*hc3
         sym=sym_e(ib)
         call neutr_frac(baryon,sym,dmuon,pfmu,delec,u,x)
         XMUON=dmuon/baryon
         XELEC=delec/baryon
         ENUC=e_nuc(ib)+sym*x**2 ! EOS
         T=PFMU/WMU  
         T2=T*T
         TERM2=DLOG(T+(T2+1.d0)**.5D0)       
         TERM1=(2.D0*T2+1.D0)*T* (T2+1)**.5D0
         emu=0.d0
         if(dmuon.ne.0.d0)emu=WMU**4*(TERM1-TERM2)/(8.D0*PI2*dmuon)
         elec=3.d0*u/4.d0
         EOS=ENUC+EMU*xmuon+ELEC*xelec ! total energy per baryon
         asym=2.d0*x-1.d0
         uplus=2.d0*dedr(X,BARYON)-asym*u !u_p+u_n         
         un=0.5d0*(uplus+u)
         up=0.5d0*(uplus-u)
         rmass=1.67323*10.**15*baryon/hc3
         xkfe=(3.d0*pi*pi*delec)**.333333
         xkfm=(3.d0*pi*pi*dmuon)**.333333
         xkfp=(3.d0*pi*pi*(delec+dmuon))**.333333
         xkfn=(3.d0*pi*pi*(baryon-delec-dmuon))**.333333
         xe=xkfn-xkfp-xkfe
         xm=xkfn-xkfp-xkfm
         IF(I.EQ.IMAX)WRITE(6,*)'NO SELFCONSISTENCY'
       
C***** STORE DATA ****************            
         dnum(IB)=baryon/HC3
         YNEU(IB)=X
         Y_EL(IB)=delec/baryon
         CHEM(IB)=U
         CHEMP(IB)=up
         CHEMN(IB)=un
         EN(IB)=ENUC            !average energy per nucleon
         EM(IB)=EMU             !average energy per muon
         EE(IB)=ELEC            !average energy per electron
         ETOT(IB)=EOS           !total energy per nucleon
         edens(ib)=eos*baryon/hc3
         dmass(ib)=(wn*dnum(ib) +edens(ib))*dgcm
         WRITE(6,1)BARYON/HC3,dmass(ib),1.d0-X,xelec,un,
     :        U,ENUC,emu,emu*dmuon/hc3,elec,elec*delec/hc3,
     :        xe,xm,edens(ib)
C*********** NEXT DENSITY    ********************
      enddo
      WRITE(6,*)
      WRITE(6,*)'DENSITY:'                     
      WRITE(6,2)(dnum(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'n_p:'                     
      WRITE(6,2)(1.d0-YNEU(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'n_e:'                     
      WRITE(6,2)(Y_EL(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'CHEMICAL POTENTIAL neutrons :'                     
      WRITE(6,2)(CHEMN(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'CHEMICAL POTENTIAL protons  :'                     
      WRITE(6,2)(CHEMP(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'CHEMICAL POTENTIAL electrons :'                     
      WRITE(6,2)(CHEM(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'MUON ENERGY:'                     
      WRITE(6,2)(EM(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'MUON ENERGY DENSITY:'                     
      WRITE(6,2)(EM(N)*dnum(n)*(1-yneu(n)-y_el(n)),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'ELECTRON ENERGY:'                     
      WRITE(6,2)(EE(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'ELECTRON ENERGY DENSITY :'                     
      WRITE(6,2)(EE(N)*dnum(n)*y_el(n),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'TOTAL ENERGY DENSITY :'                     
      WRITE(6,2)(edens(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'TOTAL ENERGY:'                     
      WRITE(6,2)(ETOT(N),N=1,NRHO )
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'density   pressure'
      cc='       '
      WRITE(8,*)' '
      WRITE(8,*)'total energy density' 
      WRITE(6,102)
      do n=1,nrho
         yprot=1-yneu(n)
         ymu=yprot-y_el(n)
         x1(1)=dnum(n)
         call spls3(dnum,en,nrho,x1,fx1,1,sqq,sau,2,0)
         p_nuc=x1(1)**2*fx1(1)
         call spls3(dnum,etot,nrho,x1,fx1,1,sqq,sau,2,0)
         press=x1(1)**2*fx1(1)
         p_el=dnum(n)*y_el(n)*(chem(n)-EE(n))
         p_mu=dnum(n)*(1.d0-yneu(n)-y_el(n))*(chem(n)-Em(n))
      write(8,4)dnum(n),dmass(n)/dgcm
      write(6,3)dnum(n),dmass(n),p_nuc,p_el,p_mu,p_nuc+p_el+p_mu,press
      enddo
 4    format(2E13.8)
 3    FORMAT(f8.4,6e13.5,2f8.4)
 2    FORMAT(5F10.4)
 1    FORMAT(f6.4,f7.3, 2f7.4,14f8.2)
 6    FORMAT(f4.2,f7.3,7f7.2)
 102  FORMAT('     n      rho          Pnuc         Pel     ',
     :'    Pmu           Ppart        P_tot')
 101  FORMAT('  rho   dgcm14   n_p    n_e    u_n     u_e',
     :'     ENUC   EMUON   EDMUON   ELEC   EDELEC  URCAE',
     :'   URCAM    EDENS')
      STOP      
      END  


      subroutine root(x1,x2,xacc,rtbis)
      IMPLICIT REAL*8 (a-h,o-z)      
      INTEGER JMAX
c      EXTERNAL  beta_eq
      PARAMETER (JMAX=100)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=beta_eq(x2)
      f=beta_eq(x1)
      if(f*fmid.ge.0.d0) pause 'root must be bracketed in rtbis'
      if(f.lt.0.d0)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5d0
        xmid=rtbis+dx
        fmid=beta_eq(xmid)
        if(fmid.le.0.d0)rtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.d0) return
11    continue
      pause 'too many bisections in rtbis'
      END

      REAL FUNCTION beta_eq*8(x)
      IMPLICIT REAL*8 (a-h,o-z)      
      common/leptons/wmu,wel
      common/equation/b,a_s
      COMMON/MEV/HC,HC3,wn
      PARAMETER (pi=3.141592653589793d0)
      uh2=(4.d0*a_s*x)**2
      wm2=wmu**2/(hc**2)
      we2=wel**2/(hc**2)
      xkfmu=0.d0
      if(uh2.gt.wm2)xkfmu=dsqrt(uh2-wm2)
      xkfel=0.d0
      if(uh2.gt.we2)xkfel=dsqrt(uh2-we2)
      term=(xkfmu**3+xkfel**3)/(3.d0*pi**2)
      beta_eq=b*(1.d0-x)/2.d0-term
      return
      END

      subroutine neutr_frac(baryon,sym,dmuon,pfmu,delec,u,XX)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/equation/b,a_s
      common/leptons/wmu,wel
      COMMON/MEV/HC,HC3,wn
      PARAMETER (pi=3.141592653589793d0)
      pi2=pi*pi
      r=BARYON/HC3
      a_s=sym/hc
      b=r
      x1=0.d0 
      x2=1.d0
      xacc=0.0001d0
      call root(x1,x2,xacc,asym)
      u=4.d0*sym*asym
      U2=U*U  
      U3=U2*U
      DELEC=U3/(3.D0*PI2)
      IF(U2.LE.WMU**2)THEN
         PFMU=0.D0 
      ELSE
         PFMU=(U2-WMU**2)**.5D0
      ENDIF
      DMUON=PFMU*PFMU*PFMU/(3.D0*PI2)
      xx=(asym+1.d0)/2.d0
      return
      end

      REAL FUNCTION dedr *8 (XX,BARYON)
      IMPLICIT REAL*8 (A-H,O-Z)
c      calc d(nE)/dn
      COMMON/MEV/HC,HC3,wn
      common/eos/sym_e(101),e_nuc(101),dens(101),ndens
      PARAMETER (ldi=101)
      real*8 ex(ldi)
      r=BARYON/HC3   
      do i=1,ndens
         ex(i)=dens(i)*(e_nuc(i)+sym_e(i)*xx**2)
      enddo
      call derivative(dens,ex,ndens,r,dedr)
      RETURN
      END    

      subroutine derivative(x,f,n,x0,dfdx)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 df(1),xx(1),a(250),b(250),f(n),x(n)
      xx(1)=x0
      call spls3(x,f,n,xx,df,1,a,b,2,0)      
      dfdx=df(1)
      RETURN
      END    
      
      subroutine interp(x,f,n,x0,f0)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 ff(1),xx(1),a(250),b(250),f(n),x(n)
      xx(1)=x0
      call spls3(x,f,n,xx,ff,1,a,b,1,0)      
      f0=ff(1)
      RETURN
      END    



      SUBROUTINE SPLS3(X,Y,N,XI,FI,M,Q,AU,IGO,ISPL)
      IMPLICIT REAL*8(A-H,O-Z)
C     ******************************************************************
C
C     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE
C     BOUNDARIES OF THE APPROXIMATION INTERVAL.
C
C     IGO = 0      BUILD UP SPLINE ONLY.
C     IGO = 1      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.
C     IGO = 2      BUILD UP SPLINE AND COMPUTE DERIVATIVES AT M POINTS.
C
C     SINGLE PRECISION VERSION.        J.GALONSKA, 15.12.1971
c                   7.95 from muether
C     ******************************************************************
C
      DIMENSION X(2), Y(2), XI(1), Q(1), AU(1), FI(1)
C
      ZERO = 0.0
      THREE = 3.0
      SIX = 6.0
      FACT = 0.1666667                      
      IF (ISPL.NE.0)  GO TO 30
C
C
      AU(1) = ZERO
      AU(N) = ZERO
      Q(1) = ZERO
      HK = X(2) - X(1)
      YSAVE = (Y(2)-Y(1)) / HK
      AUX = ZERO
      NN = N - 1
      DO 10  K = 2,NN
      HX = X(K+1) - X(K-1)
      DIVQ = (HK*Q(K-1)+HX+HX)
      HK = X(K+1) - X(K)
      YK = (Y(K+1)-Y(K)) / HK
      Q(K) = - HK / DIVQ
      AU(K) = (SIX*(YK-YSAVE)-AUX) / DIVQ
      YSAVE = YK
      AUX = AU(K) * HK
   10 CONTINUE
C
      NN2 = NN + 2
      DO 20  KK = 2,NN
      K = NN2 - KK
   20 AU(K) = Q(K) * AU(K+1) + AU(K)
C
      IF (IGO.EQ.0)  RETURN
C
C     ******************************************************************
C
C     INTERPOLATION OR COMPUTATION OF DERIVATIVES.
C
C     IGO = 1      INTERPOLATE FOR M POINTS.
C     IGO = 2      COMPUTE DERIVATIVES AT M POINTS.
C
C     ******************************************************************
C
   30 DO 100  J = 1,M
      IF (X(1).GT.XI(J))  GO TO 110
      IF (XI(J).GT.X(N))  GO TO 120            
      M1 = 1
      M2 = N
   50 M3 = (M2+M1)/2
      IF (XI(J).GE.X(M3))  GO TO 70
      M2 = M3
      GO TO 80
   70 M1 = M3
   80 IF (M1+1-M2.NE.0)  GO TO 50
   90 DIJ = X(M2) - XI(J)
      DIM1J = X(M1) - XI(J)
      HI = X(M2) - X(M1)
      HI2 = HI * HI
      IF (IGO.GE.2)  GO TO 95
      DIJ3 = DIJ * DIJ * DIJ
      DIM1J3 = DIM1J * DIM1J * DIM1J
      FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(SIX*Y(M1)-HI2*AU(M1))
     1        *DIJ-(SIX*Y(M2)-HI2*AU(M2))*DIM1J) / HI
      GO TO 100
   95 FI(J) = FACT * (THREE*(AU(M2)*DIM1J*DIM1J-AU(M1)*DIJ*DIJ)
     1       -SIX*(Y(M1)-Y(M2))+HI2*(AU(M1)-AU(M2))) / HI
      GOTO 100
  110 M1 = 1
      M2 = 2
      GO TO 90
C
  120 M1 = N - 1
      M2 = N
      GO TO 90
  100 CONTINUE
      RETURN
C
C
      END

      subroutine interp1d(x,fx,np,npol,x1,fx1)
c----------------------------------------------------------------------------- 
c     input arrays: x(np),fx(np)   ->function f(x)  
c           value  x1 the point 
c     output value   fx1           -> f(x1)
c     mp: x(mp)< x1 <x(mp+1) 
c     mp=1 if x1 < x(1)
c     mp=np if x1 > x(np)
c     npol= # of interpolation points
c------------------------------------------------------------------------------
      implicit none
      real*8 x(np),fx(np),x1,fx1,df,a,b
      integer mp,np,l,npol,i
      a=dabs(x1-x(1))
      b=dabs(x1-x(np))
c      if(x1.eq.x(1))then
c         fx1=fx(1)
c      elseif(x1.eq.x(np))then
c         fx1=fx(np)
      if(a.le.0.0001d0)then
         fx1=fx(1)
      elseif(b.le.0.0001d0)then
         fx1=fx(np)
      else
         CALL LOCATE(x,NP,x1,MP)
         L=MIN(MAX(MP-(npol-1)/2,1),NP+1-npol)
         if(mp.eq.0.or.mp.eq.np)then
            if(mp.eq.0)then
               write(6,*)'warning x1= ',x1,' and  x(',1,')= ',x(1)
            else
               write(6,*)'warning x1= ',x1,' and  x(',np,')= ',x(np)
            endif
            write(6,*)np,npol
            do i=1,np
               write(6,*)i,x(i),fx(i)
            enddo
         endif
         CALL POLINT(x(L),fx(L),npol,x1,fx1,df)
      endif
      return
      end




      SUBROUTINE locate(xx,n,x,j)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION xx(n)
      jl=0
      ju=n+1
10    IF(ju-jl.gt.1) THEN 
         jm=(ju+jl)/2
         IF((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) THEN           
            jl=jm
         ELSE
            ju=jm
         ENDIF
      GOTO 10
      ENDIF
      j=jl
      RETURN
      END


      SUBROUTINE polint (xa,ya,n,x,y,dy)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION xa(n),ya(n),c(101),d(101)
      PARAMETER (NMAX=101)
      ns=1
      dif=dabs(x-xa(1))
      DO i=1,n
         dift = dabs(x-xa(i))
         IF(dift.lt.dif) THEN 
            ns=i
            dif=dift
         ENDIF 
         c(i)=ya(i)
         d(i)=ya(i)
      ENDDO 
      y=ya(ns)
      ns=ns-1
      DO m=1,n-1
         DO i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den = ho - hp
            IF(den.eq.0.d0) PAUSE 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         ENDDO
         IF (2*ns.lt.n-m) THEN
            dy=c(ns+1)
         ELSE 
            dy=d(ns)
            ns=ns-1
         ENDIF 
         y=y+dy
      ENDDO
      RETURN 
      END
       





