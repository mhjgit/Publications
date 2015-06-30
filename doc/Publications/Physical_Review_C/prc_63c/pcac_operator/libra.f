C      FUNCTION PAR(A)
C      FUNCTION IPAR(I)
C      FUNCTION DELTA(X)
C      FUNCTION DELTAI(I)
C      FUNCTION IDELTA(I)
C      FUNCTION TR(A1,A2,A3)
C      FUNCTION ITR(IA1,IA2,IA3)
C      FUNCTION BINOM(A,B)
C      FUNCTION FACS(X)
C      FUNCTION FAC(X)
C      FUNCTION FACI(IX)
C      FUNCTION FACIS(IX)
C      FUNCTION DFAC(A)

C      FUNCTION CG(X1,X2,X3,X4,X5,X6)
C      FUNCTION TJ(A1,A2,A3,B1,B2,B3)
C      FUNCTION SJ(A1,A2,B2,B1,A3,B3)
C      FUNCTION W(X1,X2,X3,X4,X5,X6)
C      FUNCTION QJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
C      FUNCTION XLSJJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)

C      FUNCTION CLEBR(A,B,C,D,E,F)
C      FUNCTION RACAH(JAD,JBD,JCD,JDD,JED,JFD)

C      FUNCTION WUE(IA,L)
C      FUNCTION WUM(IA,L)

C      FUNCTION FINDEX(XN,XL,XJ)
C      FUNCTION FINDEXZ(XN,XL,XJ)
C      SUBROUTINE INDEX
C                         COMMON/INDEX/XL(3,200)
C      SUBROUTINE INDEXZ
C                         COMMON/INDEXZ/XL(3,200)

      FUNCTION PAR(A)
      I = A
      PAR = -1.
      IF(I.EQ.2*(I/2)) PAR = 1.
      RETURN
      END

      FUNCTION IPAR(I)
      IPAR = -1
      IF(I.EQ.2*(I/2)) IPAR = 1
      RETURN
      END


      FUNCTION DELTA(X)
      IF(X) 10,20,10
20    DELTA = 1.0
      RETURN
10    DELTA = 0.0
      RETURN
      END

      FUNCTION DELTAI(I)
      DELTAI = 0.
      IF(I.EQ.0) DELTAI = 1.
      RETURN
      END

      FUNCTION IDELTA(I)
      IF(I) 10,20,10
20    IDELTA = 1
      RETURN
10    IDELTA = 0
      RETURN
      END

      FUNCTION TR(A1,A2,A3)
      IF(A1.LT.0.) GO TO 2
      IF(A2.LT.0.) GO TO 2
      IF(A3.LT.0.) GO TO 2
      AMIN = ABS(A1-A2)
      AMAX = A1+A2
      XA = AMIN
10    CONTINUE
      IF(A3.EQ.XA) GO TO 20
      XA = XA+1.
      IF(XA.GT.AMAX) GO TO 2
      GO TO 10
20    TR=1.
      RETURN
2     TR=-1.
      RETURN
      END

      FUNCTION ITR(IA1,IA2,IA3)
      IF(IA1.LT.0.) GO TO 2
      IF(IA2.LT.0.) GO TO 2
      IF(IA3.LT.0.) GO TO 2
      IAMIN = IABS(IA1-IA2)
      IAMAX = IA1+IA2
      IXA = IAMIN
10    CONTINUE
      IF(IA3.EQ.IXA) GO TO 20
      IXA = IXA+1
      IF(IXA.GT.IAMAX) GO TO 2
      GO TO 10
20    ITR=1
      RETURN
2     ITR=-1
      RETURN
      END

        FUNCTION BINOM(A,B)
        IF(A.EQ.B)GO TO 12
        IF(B.EQ.0.)GO TO 12
        IF(A-B.GT.B)GO TO 9
        M=A-B
        GO TO 11
9       M=B
11      BINOM=1.
        X=A-FLOAT(M)
        DO 10 I=1,M
        AI=I
10      BINOM=BINOM*(X+AI)/AI
        RETURN
12      BINOM=1.
        RETURN
        END

      FUNCTION DFACN(A)
      DFACN = DFAC(A)
      RETURN
      END

      FUNCTION FACS(X)
      DIMENSION FACS_SAVE(100)
      IS = X+1
      IF(X) 3,4,4
3     FACS = 0.0
      RETURN
4     FACS = FACS_SAVE(IS)
      IF(FACS.NE.0.) RETURN
      S = 1.
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*SQRT(Y)
2     FACS = S
      FACS_SAVE(IS) = FACS
      RETURN
      END

      FUNCTION FAC(X)
      DIMENSION FAC_SAVE(100)
      IS = X+1
      IF(X) 3,4,4
3     FAC = 0.0
      RETURN
4     FAC = FAC_SAVE(IS)
      IF(FAC.NE.0.) RETURN
      S = 1.
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*Y
2     FAC = S
      FAC_SAVE(IS) = FAC
      RETURN
      END

      FUNCTION FACI(IX)
      X = IX
      FACI = FAC(X)
      RETURN
      END

      FUNCTION FACIS(IX)
      X = IX
      FACIS = FACS(X)
      RETURN
      END

      FUNCTION DFAC(A)
      DIMENSION DFAC_SAVE(100)
      IS = A
      IF(A)6,7,8
6     DFAC=0.
      RETURN
7     DFAC=1.
      RETURN
8     N=A
      IF(PAR(A).EQ.1.) TYPE 77
77    FORMAT(1X,'ERROR IN DFAC: A MUST BE ODD')
      IF(PAR(A).EQ.1.) STOP
      DFAC = DFAC_SAVE(IS)
      IF(DFAC.NE.0.) RETURN
      DFAC=1.
      DO 9 I=1,N,2
      AI=I
9     DFAC=AI*DFAC
      DFAC_SAVE(IS) = DFAC
      RETURN
      END

      FUNCTION CG(X1,X2,X3,X4,X5,X6)
      CG = CLEBR(X1,X2,X3,X4,X5,X6)
      RETURN
      END

      FUNCTION TJ(A1,A2,A3,B1,B2,B3)
      TJ = PAR(A1-A2-B3)*CLEBR(A1,B1,A2,B2,A3,-B3)/SQRT(2.*A3+1.)
      RETURN
      END

      FUNCTION SJ(A1,A2,B2,B1,A3,B3)
      SJ = PAR(A1+A2+A3+B1)*RACAHR(A1,A2,A3,B1,B2,B3)
      RETURN
      END

      FUNCTION W(X1,X2,X3,X4,X5,X6)
      W = RACAHR(X1,X2,X3,X4,X5,X6)
      RETURN
      END

      FUNCTION QJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      QJ = 0.
      J1 = 2.*X1
      J2 = 2.*X2
      J3 = 2.*X3
      J4 = 2.*X4
      J5 = 2.*X5
      J6 = 2.*X6
      J7 = 2.*X7
      J8 = 2.*X8
      J9 = 2.*X9
      KMIN=MAX0(IABS(J1-J9),IABS(J2-J6),IABS(J4-J8))
      KMAX=MIN0(J1+J9,J2+J6,J4+J8)
      IF (KMIN.GT.KMAX) RETURN
      DO 100 K=KMIN,KMAX,2
      XK = K/2.
      A=SJ(X1,X4,X7,X8,X9,XK)
      B=SJ(X2,X5,X8,X4,XK,X6)
      C=SJ(X3,X6,X9,XK,X1,X2)
      QJ=QJ+PAR(2.*XK)*(2.*XK+1.)*A*B*C
100   CONTINUE
      RETURN
      END

      FUNCTION XLSJJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      X = (2.*X3+1.)*(2.*X6+1.)*(2.*X7+1.)*(2.*X8+1.)
      XLSJJ = SQRT(X)*QJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      RETURN
      END

C
C ***********************  CLEB, CLEBI, CLEB2I ************************
C
      FUNCTION CLEBR(A,B,C,D,E,F)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
      PARAMETER (LFACTC=200)
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
      LOGICAL FIRST
C      REAL*8 FACLOG
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=NINT(2.*A)
      IB=NINT(2.*C)
      IC=NINT(2.*E)
      ID=NINT(2.*B)
      IE=NINT(2.*D)
      IF=NINT(2.*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBI(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEB(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL TH_FACINIT
      RAC=0.0
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1.0
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1.0)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1.)**((IAM-IDM)/2)*SQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1.0)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5*(ALOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0.0
      S1=(-1.0)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*EXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEB=RAC
      RETURN
      END
      FUNCTION RACAH(JAD,JBD,JCD,JDD,JED,JFD)
C
C        CALCULATES RACAH COEFFICIENTS
C
C        ENTRIES : RACAH , RACAHI , RACAHR
C            RACAH  : INTEGER ARGUMENTS = 2*J
C            RACAHI : ARGUMENTS = TRUE INTEGER VALUE
C            RACAHR : ARGUMENTS = TRUE REAL VALUE
C        EXTERNAL : TH_FACINIT , GENERATES FACTORIAL TABLE
C
      DIMENSION I(16)
      LOGICAL FIRST
C      REAL*8 G,S
      COMMON / LOGFAC / FIRST,LFACT,G(1)
      EQUIVALENCE(I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),
     1 (I(6),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11),
     2 (I(12),I12),(I(13),I13),(I(14),I14),(I(15),I15),(I(16),I16)
C        MAKE USEFULL COMBINATIONS
      K=JAD+JBD-JED+2
      I1=K/2
      IF((2*I1).NE.K) GOTO 300
      K=JCD+JDD-JED+2
      I4=K/2
      IF((2*I4).NE.K) GOTO 300
      K=JAD+JCD-JFD+2
      I7=K/2
      IF((2*I7).NE.K) GOTO 300
      K=JBD+JDD-JFD+2
      I10=K/2
      IF((2*I10).NE.K) GOTO 300
      I13=I1+JED
      I14=I4+JED
      I15=I7+JFD
      I16=I10+JFD
      I2=I13-JAD
      I3=I13-JBD
      I5=I14-JCD
      I6=I14-JDD
      I8=I15-JAD
      I9=I15-JCD
      I11=I16-JBD
      I12=I16-JDD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,2,2
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    2 IL=MAX(I13,I14,I15,I16)
      IF(MIN(JAD,JBD,JCD,JDD,JED,JFD)) 300,20,1
C    ..............
      ENTRY RACAHI(JA1,JB1,JC1,JD1,JE1,JF1)
C        MAKE USEFULL COMBINATIONS
      I13=JA1+JB1+JE1+1
      I14=JC1+JD1+JE1+1
      I15=JA1+JC1+JF1+1
      I16=JB1+JD1+JF1+1
      I1=I13-JE1*2
      I2=I13-JA1*2
      I3=I13-JB1*2
      I4=I14-JE1*2
      I5=I14-JC1*2
      I6=I14-JD1*2
      I7=I15-JF1*2
      I8=I15-JA1*2
      I9=I15-JC1*2
      I10=I16-JF1*2
      I11=I16-JB1*2
      I12=I16-JD1*2
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,4,4
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    4 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA1,JB1,JC1,JD1,JE1,JF1)
      IF(LMIN)300,20,1
C     ............
      ENTRY RACAHR(A,B,C,D,E,F)
C     CONVERT ARGUMENTS TO INTEGER 
      JA=NINT(2.*A)
      JB=NINT(2.*B)
      JC=NINT(2.*C)
      JD=NINT(2.*D)
      JE=NINT(2.*E)
      JF=NINT(2.*F)
C        MAKE USEFULL COMBINATIONS
      K=JA+JB-JE+2
      I1=K/2
      IF((2*I1-K).NE.0) GOTO 300
      K=JC+JD-JE+2
      I4=K/2
      IF((2*I4-K).NE.0) GOTO 300
      K=JA+JC-JF+2
      I7=K/2
      IF((2*I7-K).NE.0) GOTO 300
      K=JB+JD-JF+2
      I10=K/2
      IF((2*I10-K).NE.0) GOTO 300
      I13=I1+JE
      I14=I4+JE
      I15=I7+JF
      I16=I10+JF
      I2=I13-JA
      I3=I13-JB
      I5=I14-JC
      I6=I14-JD
      I8=I15-JA
      I9=I15-JC
      I11=I16-JB
      I12=I16-JD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,3,3
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    3 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA,JB,JC,JD,JE,JF)
      IF(LMIN)300,20,1
C      ------------
    1 IF(FIRST) CALL TH_FACINIT
      IF(IL.GE.LFACT) STOP 'RACAH: LENGTH FACTORIAL TABLE INSUFFICIENT'
      J1=IL-I13+1 
      J2=IL-I14+1 
      J3=IL-I15+1 
      J4=IL-I16+1
      J5=I13+I4-IL 
      J6=I15+I5-IL 
      J7=I16+I6-IL
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      H=PH*EXP ((G(I1)+G(I2)+G(I3)-G(I13+1)+G(I4)+G(I5)+G(I6)-
     1G(I14+1)+G(I7)+G(I8)+G(I9)-G(I15+1)+G(I10)+G(I11)+G(I12)-G(I16+1))
     2*.5+G(IL+1)-G(J1)-G(J2)-G(J3)-G(J4)-G(J5)-G(J6)-G(J7))
      IF(N)300,110,120
C
  110 RACAH=H 
      RETURN
C
  120 S=1.
      K=N-1
      KL=IL+1
      J5=J5-1
      J6=J6-1
      J7=J7-1
      DO 130 J=1,N   
C K=N-J
      S=1.-((KL+K)*(J5-K)*(J6-K)*(J7-K))*S/((J1+K)*(J2+K)*(J3+K)*(J4+K))
      K=K-1
  130 CONTINUE  
      RACAH=H*S
      RETURN
C
C      ONE OF THE ARGUMENTS =0
   20 IAD=IL
      IBD=IL
      DO 21 J=13,16
      IF(IAD.LT.I(J)) GOTO 22
      IF(IAD.LT.IBD) IBD=IAD
      IAD=I(J)
      GOTO 21
   22 IF(IBD.GT.I(J)) IBD=I(J)
   21 CONTINUE
      J5=I13+I4-IL 
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      RACAH=PH/SQRT(FLOAT(IAD*IBD))
      RETURN
C
C      IMPOSSIBLE COMBINATION OF ARGUMENTS
  300 RACAH=0.
      RETURN
      END

      SUBROUTINE TH_FACINIT
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
      DATA FIRST/.TRUE./,LFACT/LFACTC/
C      REAL*8 FACLOG,FN
      FIRST=.FALSE.
      FACLOG(1)=0.0
      FACLOG(2)=0.0
      FN=1.0
      DO 10 I=3,LFACTC
      FN=FN+1.0
      FACLOG(I)=FACLOG(I-1)+LOG(FN)
   10 CONTINUE
      RETURN
      END


      FUNCTION WUE(IA,L)
      PI = 3.14159
      XL = L
      A = IA
      AEL = (1.2*(A**(1./3.)))**(2*XL)
      Y = 3./(3.+XL)
      WUE = (1./(4.*PI))*Y*Y*AEL
      RETURN
      END

      FUNCTION WUM(IA,L)
      PI = 3.14159
      XL = L
      A = IA
      AML = (1.2*(A**(1./3.)))**(2*XL-2.)
      Y = 3./(3.+XL)
      WUM = (10./(PI))*Y*Y*AML
      RETURN
      END

      FUNCTION XXM(K)
      M = 0.99+(SQRT(8.*K+1.)-1.)/2.
      XXM = M
      RETURN
      END

      FUNCTION XXL(K)
      M = XXM(K)
      L = -K+M*(M+1)/2
      IF((L+M).EQ.2*((L+M)/2)) L = L + 1
      XXL = L
      RETURN
      END

      FUNCTION XXJ(K)
      M = XXM(K)
      XXJ=0.5-K+M*(M+1)/2.
      RETURN
      END

      FUNCTION XXN(K)
CBAB N = 1,2,3,4...
      M = XXM(K)
      L = XXL(K)
      XXN = 1+(M-L-1)/2
      RETURN
      END

      FUNCTION FINDEX(XN,XL,XJ)
CBAB XN STARTS AT 1
      FINDEX = 0
      IF(XN.EQ.0.) RETURN
      IF(TR(XL,XJ,0.5).EQ.-1.) RETURN
      N = XN-1.
      L = XL
      J2 = 2.*XJ
      FINDEX = ((2*N+L)*(2*N+L+3)-(J2-1)+2)/2
      RETURN
      END

      FUNCTION IFINDEXI(NN,L,J2)
CBAB NN STARTS AT 1
      IFINDEXI = 0
      IF(NN.EQ.0.) RETURN
      IF(J2.GT.2*L+1) RETURN
      IF(J2.LT.2*L-1) RETURN
      N = NN-1
      IFINDEXI = ((2*N+L)*(2*N+L+3)-(J2-1)+2)/2
      RETURN
      END

      FUNCTION FINDEXZ(XN,XL,XJ)
CBAB XN STARTS AT 0
      FINDEX = 0
      IF(XN.LT.0.) RETURN
      IF(TR(XL,XJ,0.5).EQ.-1.) RETURN
      N = XN
      L = XL
      J2 = 2.*XJ
      FINDEX = ((2*N+L)*(2*N+L+3)-(J2-1)+2)/2
      RETURN
      END

      SUBROUTINE INDEXZ
      COMMON/INDEXZSK/XL(3,200)
      I = 0
      DO 100 I1 = 1,1000
      N = I1-1
      LN = 0
      L = N-2*LN
200   IF(L.LT.0) GO TO 100
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L+1
      XL(1,I) = LN
      XL(2,I) = L
      XL(3,I) = XJ2/2.
      IF(L.EQ.0) GO TO 201
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L-1
      XL(1,I) = LN
      XL(2,I) = L
      XL(3,I) = XJ2/2.
201   L = L-2
      LN = LN + 1
      GO TO 200
100   CONTINUE
101   CONTINUE
      CALL LABELN
      RETURN
      END

      SUBROUTINE INDEX
      COMMON/INDEXSK/XL(3,200)
      I = 0
      DO 100 I1 = 1,1000
      N = I1-1
      LN = 0
      L = N-2*LN
200   IF(L.LT.0) GO TO 100
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L+1
      XL(1,I) = LN+1
      XL(2,I) = L
      XL(3,I) = XJ2/2.
      IF(L.EQ.0) GO TO 201
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L-1
      XL(1,I) = LN+1
      XL(2,I) = L
      XL(3,I) = XJ2/2.
201   L = L-2
      LN = LN + 1
      GO TO 200
100   CONTINUE
101   CONTINUE
      CALL LABELN
      RETURN
      END


      SUBROUTINE LABELN
      CHARACTER*2 ZLABEL 
      CHARACTER*1 LABEL
      COMMON/CLABEL/LABEL(20)
      COMMON/ZLABEL/ZLABEL(200)

      LABEL(1) = 's'
      LABEL(2) = 'p'
      LABEL(3) = 'd'
      LABEL(4) = 'f'
      LABEL(5) = 'g'
      LABEL(6) = 'h'
      LABEL(7) = 'i'
      LABEL(8) = 'j'
      LABEL(9) = 'k'
      LABEL(10) = 'l'
      LABEL(11) = 'm'
      LABEL(12) = 'n'
      LABEL(13) = 'o'
      LABEL(14) = 'p'
      LABEL(15) = 'q'
      LABEL(16) = 'r'
      LABEL(17) = 's'
      LABEL(18) = 't'
      LABEL(19) = 'u'
      LABEL(20) = 'v'
      ZLABEL(1) = 'H'
      ZLABEL(2) = 'He'
      ZLABEL(3) = 'Li'
      ZLABEL(4) = 'Be'
      ZLABEL(5) = 'B'
      ZLABEL(6) = 'C'
      ZLABEL(7) = 'N'
      ZLABEL(8) = 'O'
      ZLABEL(9) = 'F'
      ZLABEL(10) = 'Ne'
      ZLABEL(11) = 'Na'
      ZLABEL(12) = 'Mg'
      ZLABEL(13) = 'Al'
      ZLABEL(14) = 'Si'
      ZLABEL(15) = 'P'
      ZLABEL(16) = 'S'
      ZLABEL(17) = 'Cl'
      ZLABEL(18) = 'Ar'
      ZLABEL(19) = 'K'
      ZLABEL(20) = 'Ca'
      ZLABEL(21) = 'Sc'
      ZLABEL(22) = 'Ti'
      ZLABEL(23) = 'V'
      ZLABEL(24) = 'Cr'
      ZLABEL(25) = 'Mn'
      ZLABEL(26) = 'Fe'
      ZLABEL(27) = 'Co'
      ZLABEL(28) = 'Ni'
      ZLABEL(29) = 'Cu'
      ZLABEL(30) = 'Zn'
      ZLABEL(31) = 'Ga'
      ZLABEL(32) = 'Ge'
      ZLABEL(33) = 'As'
      ZLABEL(34) = 'Se'
      ZLABEL(35) = 'Br'
      ZLABEL(36) = 'Kr'
      ZLABEL(37) = 'Rb'
      ZLABEL(38) = 'Sr'
      ZLABEL(39) = 'Y'
      ZLABEL(40) = 'Zr'
      ZLABEL(41) = 'Nb'
      ZLABEL(42) = 'Mo'
      ZLABEL(43) = 'Tc'
      ZLABEL(44) = 'Ru'
      ZLABEL(45) = 'Rh'
      ZLABEL(46) = 'Pd'
      ZLABEL(47) = 'Ag'
      ZLABEL(48) = 'Cd'
      ZLABEL(49) = 'In'
      ZLABEL(50) = 'Sn'
      ZLABEL(51) = 'Sb'
      ZLABEL(52) = 'Te'
      ZLABEL(53) = 'I'
      ZLABEL(54) = 'Xe'
      ZLABEL(55) = 'Cs'
      ZLABEL(56) = 'Ba'
      ZLABEL(57) = 'La'
      ZLABEL(58) = 'Ce'
      ZLABEL(59) = 'Pr'
      ZLABEL(60) = 'Nd'
      ZLABEL(61) = 'Pm'
      ZLABEL(62) = 'Sm'
      ZLABEL(63) = 'Eu'
      ZLABEL(64) = 'Gd'
      ZLABEL(65) = 'Tb'
      ZLABEL(66) = 'Dy'
      ZLABEL(67) = 'Ho'
      ZLABEL(68) = 'Er'
      ZLABEL(69) = 'Tm'
      ZLABEL(70) = 'Yb'
      ZLABEL(71) = 'Lu'
      ZLABEL(72) = 'Hf'
      ZLABEL(73) = 'Ta'
      ZLABEL(74) = 'W'
      ZLABEL(75) = 'Re'
      ZLABEL(76) = 'Os'
      ZLABEL(77) = 'Ir'
      ZLABEL(78) = 'Pt'
      ZLABEL(79) = 'Au'
      ZLABEL(80) = 'Hg'
      ZLABEL(81) = 'Tl'
      ZLABEL(82) = 'Pb'
      ZLABEL(83) = 'Bi'
      ZLABEL(84) = 'Po'
      ZLABEL(85) = 'At'
      ZLABEL(86) = 'Rn'
      ZLABEL(87) = 'Fr'
      ZLABEL(88) = 'Ra'
      ZLABEL(89) = 'Ac'
      ZLABEL(90) = 'Th'
      ZLABEL(91) = 'Pa'
      ZLABEL(92) = 'U'
      ZLABEL(93) = 'Np'
      ZLABEL(94) = 'Pu'
      ZLABEL(95) = 'Am'
      ZLABEL(96) = 'Cm'
      ZLABEL(97) = 'Bk'
      ZLABEL(98) = 'Cf'
      ZLABEL(99) = 'Es'
      ZLABEL(100) = 'Fm'
      ZLABEL(101) = 'Md'
      ZLABEL(102) = 'No'
      ZLABEL(103) = 'Lr'
      ZLABEL(104) = 'X1'
      ZLABEL(105) = 'X2'
      ZLABEL(106) = 'X3'
      ZLABEL(107) = 'X4'
      ZLABEL(108) = 'X5'
      ZLABEL(109) = 'X6'
      END

