      SUBROUTINE MU_CAPT(Q,IZz,IAa,IU,IW,ISIGN,IME,
     :                   NN,LN,JN,NP,LP,JP,TERM1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/INTTI/M,XMAX,IZ,IA,IO,RBEE,RKAP
      DATA ALPHA,RMASS/7.30D-3,9.39D-1/
      DATA RMUM/1.0566D-1/              ! muon mass in GeV
      M=100 
      XMAX=15.
      iz=izz
      ia=iaa
      IO=0
      JF=IU
      CALL JBINOM
      RHBARW=45.0*(FLOAT(IA)**(-0.3333333))-25.0*(FLOAT(IA)**
     &       (-0.6666667))
      RBEE=197.0/SQRT(940.0*RHBARW)    ! b in fm
      REDM=RMUM*RMASS*DBLE(IA)*5.07D0/(RMUM+RMASS*DBLE(IA))
      RKAP=ALPHA*DBLE(IZ)*REDM*DBLE(RBEE)
c      WRITE(6,*) 'nn ln 2jn np lp 2jp : BARE  SP MATRIX ELEMENTS:'
      TERM1=SPME(IW,IU,NN,LN,JN,NP,LP,JP,Q,IME,ISIGN,JF)
c      WRITE(6,'(6I3,F12.6)')NN,LN,JN,NP,LP,JP,TERM1

      END


C
C                       TEOR:<SUHONEN>MUONSPMET.FOR
C
******************************************************************
      REAL FUNCTION SPME*8
     &         (IW,IU,NN,LN,JN,NP,LP,JP,Q,IME,ISIGN,JF)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*  Definition of the reduced single-particle matrix elements      *
*  (j_n // T_1 // j_p) :                                          *
*                                                                 *
*  IME=1 <--> <1wu>, <1wu+/-> ; IME=2 <--> <0ww>, <0ww+/-> ;      *
*  IME=3 <--> <0wwp> ; IME=4 <--> <1wup>                          * 
*                                                                 *
*  ISIGN=0 <-> j_w ; ISIGN=1 <-> j_w + A*j_(w-1)                  *
*  ISIGN=-1 <-> j_w - A*j_(w+1)                                   *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/INTTI/M,XMAX,IZ,IA,IO,RBEE,RKAP
      DATA FPII/0.0795775/
*********************************************************************
      SPME=0.0
      IF(IU.NE.JF)RETURN
      IF((IME.EQ.2).AND.(IW.NE.IU))RETURN
      IF((ISIGN.EQ.1).AND.(IW.EQ.0))RETURN
      IF(IW.LT.0)RETURN
      IF((IME.EQ.3).AND.(IW.NE.IU))RETURN
C --------------------------------------------
      CALL GEOM(IW,IU,LN,JN,LP,JP,G1WU,G0WWP,G1WUPA,G1WUPB)
C      PRINT*,'G1WU,G0WWP,G1WUPA,G1WUPB',G1WU,G0WWP,G1WUPA,G1WUPB
      IF(IME.EQ.1)SPME=SQRT(3.0)*FPII*G1WU*
     &RINTEG(NN,LN,NP,LP,0,IW,Q,ISIGN)               ! <1wu,ISIGN>
C      PRINT*,'D110,INTEGRAL',SQRT(3.0)*FPII*G1WU,
C     &        RINTEG(NN,LN,NP,LP,0,IW,Q,ISIGN)
C                    <1wu>
      IF(IME.EQ.2)SPME=FPII*G0WWP*
     &RINTEG(NN,LN,NP,LP,0,IW,Q,ISIGN)  ! <0ww,ISIGN>
      IF(IME.LT.3)GOTO 100
      APU1=0.0
      APU2=0.0
      APU1=RINTEG(NN,LN,NP,LP,1,IW,Q,0)/RBEE
     &   -2.0*SQRT(FLOAT(NP+LP)+1.5)*RINTEG(NN,LN,NP,LP+1,0,IW,Q,0)/RBEE
      APU=0.0
      IF(LP.GE.1)APU=2.0*SQRT(FLOAT(NP+LP)+0.5)*
     &               RINTEG(NN,LN,NP,LP-1,0,IW,Q,0)/RBEE
      APU2=APU-RINTEG(NN,LN,NP,LP,1,IW,Q,0)/RBEE
C      PRINT*,'b,APU1,APU,APU2',RBEE,APU1,APU,APU2
      IF(IME.EQ.3)THEN                               ! <0wwp>
      IF(JP.EQ.(2*LP+1))THEN
         SPME=-FPII*G0WWP*APU1   ! in fm**(-1)
      ELSE
         SPME=-FPII*G0WWP*APU2   ! in fm**(-1)
      END IF
      END IF                                         ! <0wwp>
      IF(IME.EQ.4)THEN                               ! <1wup>
      SQ1=FLOAT(LP+1)/FLOAT(2*LP+1)
      SQ1=SQRT(SQ1)
      SQ2=FLOAT(LP)/FLOAT(2*LP+1)
      SQ2=SQRT(SQ2)
      PHAS=1.0
      INDI=IW+IU
      IF(MOD(INDI,2).NE.0)PHAS=-PHAS
      SPME=SQRT(3.0)*FPII*PHAS*
     &       (SQ1*APU1*G1WUPB-SQ2*APU2*G1WUPA)   ! in fm**(-1)
      END IF                                         ! <1wup>
100   RETURN
      END
*****************************************************************
*****************************************************************
      SUBROUTINE GEOM(IW,IU,LN,JN,LP,JP,G1WU,G0WWP,G1WUPA,G1WUPB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CGH0,CG,RACAH
      REAL*8 DIU,DIW,DLN,DLP,DJN,DJP
* * * * * * * * * * * * * * * * * * * * * * * * * * * *
      G1WU=0.0
      G0WWP=0.0
      G1WUPA=0.0
      G1WUPB=0.0
      DIU=DBLE(IU)
      DIW=DBLE(IW)
      DLN=DBLE(LN)
      DLP=DBLE(LP)
      DJN=DBLE(JN)/2.D0
      DJP=DBLE(JP)/2.D0
      RJPHAT=SQRT(FLOAT(JP+1))
      RJNHAT=SQRT(FLOAT(JN+1))
      RLPHAT=SQRT(FLOAT(2*LP+1))
      RLNHAT=SQRT(FLOAT(2*LN+1))
      RCGU=SNGL(CGH0(DJN,DJP,DIU,1))
      RCGW=SNGL(CGH0(DJN,DJP,DIW,1))
      RCG2=0.0
      IF(LP.GE.1)RCG2=SNGL(CG(DLN,0.D0,DLP-1.D0,0.D0,DIW))*
     &SNGL(RACAH(DLP,DLN,1.D0,DIW,DIU,DLP-1.D0))*SQRT(FLOAT(2*LP-1))
      RCG3=SNGL(CG(DLN,0.D0,DLP+1.D0,0.D0,DIW))*
     &SNGL(RACAH(DLP,DLN,1.D0,DIW,DIU,DLP+1.D0))*SQRT(FLOAT(2*LP+3))
      RCG5=SNGL(CG(DIU,1.D0,1.D0,-1.D0,DIW))
      RCG6=SNGL(CG(DIU,0.D0,1.D0,0.D0,DIW))
      RAC=SNGL(RACAH(DJP,DJN,DLP,DLN,DIU,0.5D0))
      PH1=1.0
      PH2=1.0
      PH3=1.0
      PH4=1.0
      PH5=1.0
      IND1=LN+(JN+1)/2
      IF(MOD(IND1,2).NE.0)PH1=-PH1
      IND2=LP+IU+IW+1
      IF(MOD(IND2,2).NE.0)PH2=-PH2
      IND3=IU+(JP+JN)/2
      IF(MOD(IND3,2).NE.0)PH3=-PH3
      IND4=IW+(JN+1)/2
      IF(MOD(IND4,2).NE.0)PH4=-PH4
      IND5=IND1+1+LP
      IF(MOD(IND5,2).NE.0)PH5=-PH5
      APU1=PH2*RJPHAT*RJNHAT*RCGU
      APU2=0.0
      IF(IU.GT.0)APU2=(FLOAT(JN+1)+PH3*FLOAT(JP+1))
     &                *RCG5/SQRT(FLOAT(2*IU*(IU+1)))
      APU3=PH1*RCG6
      G1WU=APU1*(APU2+APU3)
      G0WWP=-PH4*RJPHAT*RJNHAT*RCGW
      APU4=PH5*RLPHAT*RLNHAT*RJPHAT*RJNHAT*SQRT(FLOAT(2*IU+1))*RAC
      G1WUPA=APU4*RCG2    ! L=lp-1
      G1WUPB=APU4*RCG3    ! L=lp+1
      RETURN
      END





C
C                       TEOR:<SUHONEN>MUONINT.FOR
C
******************************************************************
C
* SIMPSONIAN INTGRATION OF THE FUNCTION  exp(-KAPPA*X)*F(X), WITH
*
* F(X)=X**(2+K)*j_w(qX)*g(nn,ln,X)*g(np,Lp,X), or
* F(X)=X**(2+K)*(j_w(qX)+/-(KAPPA/q)j_(w-/+1))*g(nn,ln,X)*g(np,Lp,X) 
*
      REAL FUNCTION RINTEG*8(NN,LN,NP,LP,K,IW,Q,ISIGN)
********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FUNCT,DH,DX,SUM,DMAX,RKAP,QFM,QQ
      COMMON/INTTI/M,XMAX,IZ,IA,IO,RBEE,RKAP
*********************************************************************
C      READ(5,*)NN,LN,NP,LP,K,IW      ! IW=0,1,2 
C      READ(5,*)M,XMAX      ! No. of int. points, maximum range of int.
C      READ(5,*)Q,IZ,IA,IO   ! q in MeV, data of the parent, IO-param.
* * * * * * * * * * * * * * * * * * * * * * * * * * * *
      QFM=DBLE(Q)*5.07D-3            ! q in fm**(-1)
      QQ=QFM*DBLE(RBEE)
C      PRINT*,'q (in fm**(-1))=',QFM
*********************************************
      CALL RFACTOR
      CALL RFACTA
*********************************************     
      DMAX=DBLE(XMAX)
      DH=DMAX/DBLE(M)
      SUM=0.D0
      DX=0.D0
      SUM=FUNCT(NN,LN,NP,LP,K,IW,QQ,RBEE,RKAP,DX,IO,ISIGN)+
     &    FUNCT(NN,LN,NP,LP,K,IW,QQ,RBEE,RKAP,DMAX,IO,ISIGN)
      SUM=SUM/2.0D0
      DO 100 I=1,M-1
      DX=DX+DH
      SUM=SUM+FUNCT(NN,LN,NP,LP,K,IW,QQ,RBEE,RKAP,DX,IO,ISIGN)
100   CONTINUE
      RINTEG=SNGL(SUM*DH)
CCC      WRITE(6,200)NN,LN,NP,LP,K,IW,IZ,IA,Q,RINTEG
200   FORMAT(/,1X,'nn,ln,np,lp,k,w',1X,6(I3),3X,'Z,A=',2(I3),2X,
     &       'q,I(q)=',1X,F8.3,1X,1PG14.7,/)
      RETURN
      END
***************************************************************
***************************************************************
      REAL FUNCTION FUNCT*8(NN,LN,NP,LP,K,IW,QFM,
     &                                BEE,RKAP,DX,IO,ISIGN)
*
*  ISIGN=0 <-> j_w ; ISIGN=1 <-> j_w + A*j_(w-1)
*  ISIGN=-1 <-> j_w - A*j_(w+1)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RKAP,DX,QX,BESSEL(0:10),QFM,BESFACT,FACTOR,BESADD
      REAL*8 SIGN
***************************************************************
      QX=QFM*DX
      RX=SNGL(DX)
      IF(IO.EQ.1)PRINT*,'r=',RX
      DO 100 I=1,10
      BESSEL(I)=0.D0   
100   CONTINUE
      BESSEL(0)=1.D0       ! w=0
      IF(DX.GT.1.0D-10)THEN
      BESSEL(0)=DSIN(QX)/QX                                 ! w=0
      BESSEL(1)=DSIN(QX)/(QX*QX)-DCOS(QX)/QX                ! w=1
      BESSEL(2)=3.D0*DSIN(QX)/(QX*QX*QX)-3.D0*DCOS(QX)/(QX*QX)
     &       -DSIN(QX)/QX                                   ! w=2        
      BESSEL(3)=15.D0*DSIN(QX)/(QX*QX*QX)-15.D0*DCOS(QX)/(QX*QX)
     &       -6.D0*DSIN(QX)/QX+DCOS(QX)                     ! w=3
      BESSEL(3)=BESSEL(3)/QX                                ! w=3
      BESSEL(4)=105.D0*(DSIN(QX)/QX-DCOS(QX))/(QX*QX*QX)-
     &45.D0*DSIN(QX)/(QX*QX)+10.D0*DCOS(QX)/QX+DSIN(QX)     ! w=4
      BESSEL(4)=BESSEL(4)/QX                                ! w=4
C      PRINT*,'X=',QX
C      PRINT*,'BESSEL 0,1,2,3,4=',BESSEL(0),BESSEL(1),BESSEL(2),
C     &BESSEL(3),BESSEL(4)
C ********************** OTHER METHOD
      DO 200 I=5,10
      BESSEL(I)=(2.D0*DBLE(I)-1.D0)*BESSEL(I-1)/QX-BESSEL(I-2)
C      PRINT*,'REKURSIO: I,BESSEL=',I,BESSEL(I)
200   CONTINUE
      END IF
      IF(ISIGN.EQ.0)BESFACT=BESSEL(IW)
      IF(IABS(ISIGN).GT.0)THEN
         SIGN=1.D0
         IF(ISIGN.LT.0)SIGN=-1.D0
         FACTOR=RKAP/QFM
         BESADD=0.D0
         IND=IW-ISIGN
         IF(IND.GE.0)BESADD=SIGN*FACTOR*BESSEL(IND)
         BESFACT=BESSEL(IW)+BESADD
      END IF
      FUNCT=(DX**(2+K))*BESFACT*DBLE(RADIAL(NN,LN,BEE,RX))*
     &   DBLE(RADIAL(NP,LP,BEE,RX))
      FUNCT=FUNCT*DEXP(-RKAP*DX)
      IF(IO.EQ.1)PRINT*,'BESSEL(w),R(NN,LN),R(NP,LP)=',BESSEL(IW),
     &         RADIAL(NN,LN,BEE,RX),RADIAL(NP,LP,BEE,RX)
      RETURN
      END


C
C                       TEOR:<SUHONEN>RADIAL.FOR
C
C      IMPLICIT REAL*8(A-H,T-Z)

      REAL FUNCTION RADIAL*8(N,L,BEE,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RHO,VEE(0:6,0:12)
      REAL*8 FACTO,FACT
      COMMON/FACTI/FACTO(0:15)
      COMMON/FACTORI/FACT(0:10)
      DATA PII/3.141592654/
**********************************************
C      READ(5,*)N,L,BEE,X      ! n,l,b,r
      RNU=1/(BEE*BEE)
CCC   CALL FACTOR
CCC   CALL FACTORIAL
CCC      RHO=DBLE(RNU*X*X)      ! here original units
      RHO=DBLE(X*X)             ! here scaled units for integration
      CALL VEET(RHO,VEE)
      FACTO1=SNGL(FACTO(L+N))
      FACTO2=SNGL(FACTO(L))
      FACTO3=SNGL(FACT(N))
CCC      APU1=SQRT(2.0*RNU)     ! here original units
      APU1=SQRT(2.0)            ! here scaled units for integration
      APU1=APU1**(2*L+3)
      APU2=SQRT(2.0/PII)
      RNORM=FACTO1*APU1*APU2/((2.0**N)*FACTO3)
      RNORM=SQRT(RNORM)/FACTO2
C      PRINT*,'FACTO1,FACTO2,FACTO3,APU1,APU2:',
C     & FACTO1,FACTO2,FACTO3,APU1,APU2
      VNL=SNGL(VEE(N,L))      
      RADAPU=1.0
      IF(L.GT.0)RADAPU=X**L
CCC      RADIAL=RNORM*EXP(-0.5*RNU*X*X)*RADAPU*VNL     ! original
      RADIAL=RNORM*EXP(-0.5*X*X)*RADAPU*VNL            ! scaled
C      PRINT*,'RNORM,RADAPU,VNL:',RNORM,RADAPU,VNL
C      WRITE(6,100)N,L,X,RADIAL
100   FORMAT(/,2X,'n,l,r,R_nl(r):  ',I1,2X,I1,2X,F8.3,5X,1PG15.8)
      RETURN
      END
C
      SUBROUTINE VEET(X,VEE)
* CALCULATE THE V(N,L,X) FUNCTIONS OF HORIE & SASAKI (RELATED
* TO THE ASSOCIATED LAGUERRE POLYNOMIALS) BY USING THE FOLLOWING
* RECURSION RELATIONS:
*     V(N,L-1,X)=V(N-1,L-1,X)-2*X*V(N-1,L,X)/(2*L+1),
*     V(N,L,X)=((2*L+1)*V(N,L-1,X)+2*N*V(N-1,L,X))/(2*N+2*L+1).
* STORE THESE IN TO THE MATRIX VEE(N,L).
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,VEE(0:6,0:12)
*****************************************************************
      DO 100 L=0,12
      VEE(0,L)=1.D0
100   CONTINUE
      N=1
200   CONTINUE
      DO 300 L=1,12
      VEE(N,L-1)=VEE(N-1,L-1)-2.D0*X*VEE(N-1,L)/DFLOAT(2*L+1)
      VEE(N,L)=(DFLOAT(2*L+1)*VEE(N,L-1)+DFLOAT(2*N)*
     &         VEE(N-1,L))/DFLOAT(2*N+2*L+1)
300   CONTINUE
      N=N+1
      IF(N.GE.7)GOTO 400
      GOTO 200
400   RETURN
      END
*
******************************************************************
*          CALCULATE THE FACTORIALS (2L+1)!! , L=0,...,15 .      *
******************************************************************
      SUBROUTINE RFACTOR
c      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8(A-H,T-Z)
      COMMON/FACTI/FACTO(0:15)
      FACTO(0)=1.D0
      FACTO(1)=3.D0
      DO 100 I=2,15
      FACTO(I)=DBLE(2*I+1)*FACTO(I-1)
100   CONTINUE
      RETURN
      END
*
******************************************************************
*          CALCULATE THE FACTORIALS n! , n=0,...,10 .      *
******************************************************************
      SUBROUTINE RFACTA
c     IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8(A-H,T-Z)
      COMMON/FACTORI/FACT(0:10)
      FACT(0)=1.D0
      FACT(1)=1.D0
      DO 100 I=2,10
      FACT(I)=DBLE(I)*FACT(I-1)
100   CONTINUE
      RETURN
      END



C
C             DELTA.GEOM
C
      REAL FUNCTION CGH0*8(FJ1,FJ2,FK,IX)
      IMPLICIT REAL*8 (A-H,O-Z)
C     CLEBSH-GORDAN COEF.         
C     IX=1 -- (FJ1,0.5,FJ2,-0.5/FK,0)
C     IX=2 -- (FJ1,  0,FJ2,   0/FK,0)
C     IMPLICIT REAL*8(A-H,O-Z)       
      COMMON /CFACTR/FC(30),FCH(30),ICF
      CGH0=0.D0                        
      IF(DMIN1(FJ1+FJ2-FK,FK-DABS(FJ1-FJ2)).LT.-0.1D0) RETURN
      IF(ICF-12345) 1,2,1                                    
    1 CONTINUE                                               
      F=1.D0                                                 
      FH=0.5D0                                               
      DO 10 I=1,30                                           
      FI=I                                                   
      FC(I)=F                                                
      FCH(I)=FH                                              
      F=F*FI                                                 
   10 FH=FH*(FI+0.5D0)                                       
      ICF=12345                                              
    2 CONTINUE                                               
      IG2=FJ1+FJ2+FK+0.1D0                                   
      IG=IG2/2                                               
      G=IG                                                   
      K=FK+0.1D0                                             
      IG1=IG+1                                               
      IGK=IG-K+1                                             
      FGK=FCH(IGK)/FC(IGK)                                   
      SN=(-1)**(IG-K)                                        
      GO TO (20,30),IX                                       
   20 J1=FJ1+0.51D0                                          
      J2=FJ2+0.51D0                                          
      FG=FC(IG1+1)/FCH(IG1)*FGK                              
      IF(IG2.NE.2*IG) GO TO 25                               
      IGJ1=IG-J1+1                                           
      IGJ2=IG-J2+1                                           
      F1=FCH(IGJ1)/FC(IGJ1)*FCH(IGJ2)/FC(IGJ2)*FG/(G+1)      
      CGH0=SN*SQRT((FK+0.5D0)/(J1*J2*(G-FK+0.5D0))*F1)      
      RETURN                                                 
   25 IGJ1=IG-J1+2                                           
      IGJ2=IG-J2+2                                           
      F1=FCH(IGJ1)/FC(IGJ1)*FCH(IGJ2)/FC(IGJ2)*FG            
      CGH0=SN*SQRT((FK+0.5D0)/(J1*J2*(G-FJ1+1)*(G-FJ2+1))*F1)
      RETURN                                                  
   30 IF(IG2.NE.2*IG) RETURN                                  
      IGJ1=IG-FJ1+1.1D0                                       
      IGJ2=IG-FJ2+1.1D0                                       
      F1=FC(IG1)/FCH(IG1)*FCH(IGJ1)/FC(IGJ1)*FCH(IGJ2)/FC(IGJ2) 
      CGH0=SN*SQRT((FK+0.5D0)/(G-FJ1+0.5D0)/(G-FJ2+0.5D0)/(G-FK+0.5D0)
     1  *FGK*F1)                                                       
      RETURN                                                           
      END                                                              
      REAL FUNCTION RACAH*8(A,B,C,D,E,F) 
c      RACAH COEFFICIENT                                                
      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON /ZBINOM/BB(50,50)                                         
      RACAH=0.                                                         
      IF(DMIN1(E-ABS(A-B),E-ABS(C-D),F-ABS(A-C),F-ABS(B-D)).LT.-0.1)   
     1  RETURN                                                         
      A1=A+1.01D0                                                      
      D1=D+1.01D0                                                      
      F1=F+1.01D0                                                      
      I1=A1+B-E                                                        
      I2=C+D1-E                                                        
      I3=A1+C-F                                                        
      I4=B+D1-F                                                        
      IF(MIN0(I1,I2,I3,I4).LT.1) RETURN                                
      IZN=DMAX1(0.0D0,A+D-E-F,B+C-E-F)+1.01D0                          
      IZX=MIN0(I1,I2,I3,I4)                                            
      IF(IZN.GT.IZX) RETURN                                            
      I0=A1+B+C+D1                                                     
      J1=A1+E-B                                                        
      J3=C+F1-A                                                        
      J4=D1+F-B                                                        
      K1=A1+B+E                                                        
      K4=B+D1+F                                                        
      L1=J1+J3-1                                                       
      IA=A1+A                                                          
      ID=D1+D                                                          
      SUM=0.                                                           
      SN=(-1)**IZN                                                     
      DO 10 IZ=IZN,IZX                                                 
      SN=-SN                                                           
      IZM=IZ-1                                                         
      M1=I2-IZM                                                        
      M2=I3-IZM                                                        
      SUM=SUM+SN*BB(I1,IZ)*BB(J3,M1)*BB(J1,M2)*BB(I4,IZ)/BB(I0,IZ)     
   10 CONTINUE                                                         
      FAC=BB(I0,I1)*BB(I0,I4)*BB(L1,J1)/(K1*K4*BB(K1,IA)*BB(IA,I3)     
     1  *BB(K4,ID)*BB(ID,I2)*BB(L1,J4))                                
      RACAH=SQRT(FAC)*SUM                                             
      RETURN                                                           
      END                                                              
      SUBROUTINE JBINOM                                                
      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON /ZBINOM/BB(50,50)                                         
      DO 5 I=1,50                                                      
      DO 5 J=1,50                                                      
    5 BB(I,J)=0.D0                                                     
      BB(1,1)=1.D0                                                     
      DO 20 L=2,50                                                     
      BB(L,1)=1.D0                                                     
      BB(L,L)=1.D0                                                     
      IF(L.EQ.2) GO TO 20                                              
      LM=(L+1)/2                                                       
      S=1.0                                                            
      DO 10 M=2,LM                                                     
      LL=L-M+1                                                         
      FLL=LL                                                           
      FM1=M-1                                                          
      S=S*FLL/FM1                                                      
      BB(L,M)=S                                                        
   10 BB(L,LL)=S                                                       
   20 CONTINUE                                                         
      RETURN                                                           
      END                             
                                 
      REAL FUNCTION COEX*8(A,B,C,D,E,F,G,H,O)
C     9-J COEFFICIENT 
      IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZBINOM/BB(50,50)
      COEX=0.D0               
      X1=A+B-C                
      X2=A+D-G                
      X3=D+E-F                
      X4=B+H-E                
      X5=H+O-G                
      X6=F+O-C                
      IF(DMIN1(X1,X2,X3,X4,X5,X6,C-DABS(A-B),G-DABS(A-D),F-DABS(D-E)
     1  ,E-DABS(B-H),G-DABS(H-O),C-DABS(F-O)).LT.-0.1D0) RETURN     
      I1=X1+1.01D0                                                  
      I2=X2+1.01D0                                                  
      I3=X3+1.01D0                                                  
      I4=X4+1.01D0                                                  
      I5=X5+1.01D0                                                  
      I6=X6+1.01D0                                                  
      L1=B+C-A+1.01D0                                               
      L2=D+G-A+1.01D0                                               
      L3=D+F-E+1.01D0                                               
      L4=E+H-B+1.01D0                                               
      L5=G+H-O+1.01D0                                               
      L6=C+F-O+1.01D0                                               
      J1=A+D+H+O+2.01D0                                             
      J2=B+D+F+H+2.01D0                                             
      J3=A+B+F+O+2.01D0                                             
      K3=D+E+F+1.01D0                                               
      K5=G+H+O+1.01D0                                               
      K6=C+F+O+1.01D0                                               
      JA=A+A+1.01D0                                                 
      JB=B+B+1.01D0                                                 
      JD=D+D+1.01D0                                                 
      JF=F+F+1.01D0                                                 
      JH=H+H+1.01D0                                                 
      FKN=DMAX1(DABS(A-O),DABS(D-H),DABS(B-F))                      
      FKX=DMIN1(A+O,D+H,B+F)                                        
      KN=FKX-FKN+1.01D0                                             
      SUM=0.                                                        
      DO 100 KK=1,KN                                                
      FK=FKN+KK-1                                                   
      N1=A+O-FK+1.01D0                                              
      N2=B+F-FK+1.01D0                                              
      N3=D+H-FK+1.01D0                                              
      M1=D+G-O+FK+1.01D0                                            
      M2=E+H-F+FK+1.01D0                                            
      M3=B+C-O+FK+1.01D0                                            
      II1=A-O+FK+1.01D0                                             
      II2=B-F+FK+1.01D0                                             
      JJ1=A+O+FK+1.01D0                                             
      JJ2=B+F+FK+1.01D0                                             
      FF=(2*FK+1)/(JJ1*BB(JJ1,JA))*SQRT(BB(J1,N1)*BB(M1,L2)*       
     1  BB(J2,N2)*BB(M2,L4)*BB(J3,N1)*BB(M3,L1)/(JJ2*BB(JH,N3)*      
     2  BB(M1,L5)*BB(JJ2,JB)*BB(JD,N3)*BB(M2,I3)*BB(JF,N2)*BB(M3,L6)))
      IXN=DMAX1(0.D0,A+H-FK-G,D+O-FK-G)+1.01D0                        
      IXX=MIN0(N1,N3,I2,I5)                                           
      SX=0.D0                                                         
      ZX=(-1)**IXN                                                    
      DO 10 IX=IXN,IXX                                                
      ZX=-ZX                                                          
      MM1=N3-IX+1                                                     
      MM2=I2-IX+1                                                     
   10 SX=SX+ZX*BB(N1,IX)*BB(L2,MM1)*BB(II1,MM2)*BB(I5,IX)/BB(J1,IX)   
      IYN=DMAX1(0.D0,B+D-FK-E,F+H-FK-E)+1.01D0                        
      IYX=MIN0(N2,N3,I4,L3)                                           
      SY=0.D0                                                         
      ZY=(-1)**IYN                                                    
      DO 20 IY=IYN,IYX                                                
      ZY=-ZY                                                          
      MM1=N3-IY+1                                                     
      MM2=I4-IY+1                                                     
   20 SY=SY+ZY*BB(N2,IY)*BB(L4,MM1)*BB(II2,MM2)*BB(L3,IY)/BB(J2,IY)   
      IZN=DMAX1(0.D0,A+F-FK-C,B+O-FK-C)+1.01D0                        
      IZX=MIN0(N1,N2,I1,I6)                                           
      SZ=0.D0                                                         
      ZZ=(-1)**IZN                                                    
      DO 30 IZ=IZN,IZX                                                
      ZZ=-ZZ                                                          
      MM1=N2-IZ+1                                                     
      MM2=I1-IZ+1                                                     
   30 SZ=SZ+ZZ*BB(N1,IZ)*BB(L1,MM1)*BB(II1,MM2)*BB(I6,IZ)/BB(J3,IZ)   
      SUM=SUM+FF*SX*SY*SZ                                             
  100 CONTINUE                                                        
      FF=BB(J1,I5)*BB(J2,L3)*BB(J3,I6)/(K3*K5*K6*BB(JA,I2)*BB(K5,JH)  
     1  *BB(JB,I4)*BB(K3,JD)*BB(JA,I1)*BB(K6,JF))                     
      COEX=SUM*SQRT(FF)                                              
      RETURN                                                          
      END                                                             
      REAL FUNCTION CG*8(FJ1,FM1,FJ2,FM2,FJ)     
C    CLEBSCH-GORDAN COEFFICIENT                                      
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION FC(32)                                                
      FC(1)=1.D0                                                      
      DO 1 I=1,31                                                     
      I1=I+1                                                          
    1 FC(I1)=FC(I)*I                                                  
      CG=0.D0                                                         
      JM1=FJ1-FM1+1.01D0                                              
      JP1=FJ1+FM1+1.01D0                                              
      JM2=FJ2-FM2+1.01D0                                              
      JP2=FJ2+FM2+1.01D0                                              
      J12=FJ1+FJ2-FJ+1.01D0                                           
      J13=FJ1+FJ-FJ2+1.01D0                                           
      J23=FJ2+FJ-FJ1+1.01D0                                           
      JM=FJ-FM1-FM2+1.01D0                                            
      JP=FJ+FM1+FM2+1.01D0                                            
      IF(MIN0(JM1,JP1,JM2,JP2,J12,J13,J23,JM,JP).LT.1) RETURN         
      J123=FJ1+FJ2+FJ+2.01D0                                          
      JJM1=(FJ-FJ2+FM1)*1.01D0                                        
      JJM2=(FJ-FJ1-FM2)*1.01D0                                        
      IZX=MIN0(J12,JM1,JP2)                                           
      IZN=MAX0(0,-JJM1,-JJM2)+1                                       
      SUM=0.D0                                                        
      SN=-(-1)**IZN                                                   
      DO 10 IZ1=IZN,IZX                                               
      IZ=IZ1-1                                                        
      SUM=SUM+SN/(FC(IZ1)*FC(J12-IZ)*FC(JM1-IZ)*FC(JP2-IZ)*           
     1  FC(JJM1+IZ1)*FC(JJM2+IZ1))                                    
   10 SN=-SN                                                          
      FF=(2*FJ+1)*FC(JP1)*FC(JM1)*FC(JP2)*FC(JM2)*FC(JP)*FC(JM)       
     1  *FC(J12)*FC(J13)*FC(J23)/FC(J123)                             
      CG=SUM*SQRT(FF)                                                
      RETURN                                                          
      END                                                             
C
C          DELTA.ACOEFF
C
      REAL FUNCTION AHM*8(N,L,ND,LD,MM)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ZBINOM/BB(50,50)
      COMMON /AHMOF/CFA(50)
      DIMENSION C(50),CD(50)
C
      NP=N+1
      NDP=ND+1
      F=1.0
      G=2*(L+N)+1
      DO 10 M=1,NP
      IM=NP-M+1
      C(IM)=BB(NP,IM)*F
      F=F*G
      G=G-2.0
   10 CONTINUE
      F=1.0
      G=2*(LD+ND)+1
      DO 11 M=1,NDP
      IM=NDP-M+1
      CD(IM)=BB(NDP,IM)*F
      F=F*G
      G=G-2.0
   11 CONTINUE
      MXS=N+ND+1
      DO 20 IS=1,MXS
      SUM=0.0
      DO 21 M=1,IS
      MD=IS+1-M
      IF(M.GT.NP.OR.MD.GT.NDP) GO TO 21
      SUM=SUM+C(M)*CD(MD)
   21 CONTINUE
      CFA(IS)=(-1.0)**(IS-1)*SUM
   20 CONTINUE
      JJ=(MM-L-LD)/2+1
      AHM=CFA(JJ)
      RETURN
      END

