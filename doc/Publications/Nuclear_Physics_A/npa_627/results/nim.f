From m.h.jensen@fys.uio.no Wed Apr  9 20:26:54 1997
Received: from pat.uio.no by norsci0.nordita.dk (4.1/SMI-4.1)
	id AA23252; Wed, 9 Apr 97 20:26:47 +0200
Resent-Message-Id: <9704091826.AA23252@norsci0.nordita.dk>
Received: from ulrik.uio.no by pat.uio.no id <11529-0@pat.uio.no>;
          Wed, 9 Apr 1997 20:26:32 +0200
Return-Path: <machleid@phys.uidaho.edu>
Old-Received: from 129.101.65.111 by pat.uio.no with SMTP (PP) id 
              <11498-0@pat.uio.no>; Wed, 9 Apr 1997 20:26:06 +0200
Old-Received: by tamaluit.phys.uidaho.edu (5.61/UW-NDC Revision: 2.21 ) id 
              AA24840; Wed, 9 Apr 97 11:26:00 -0700
Date: Wed, 9 Apr 97 11:26:00 -0700
From: Ruprecht Machleidt <machleid@phys.uidaho.edu>
Message-Id: <9704091826.AA24840@tamaluit.phys.uidaho.edu>
To: m.h.jensen@fys.uio.no
Subject: updated nijmegen pot.
Cc: machleid@phys.uidaho.edu
Resent-Date: Wed, 9 Apr 1997 20:26:32 +0200
Resent-From: m.h.jensen@fys.uio.no
Resent-To: mjensen@nordita.dk
Status: R

Dear Morten,

find below a revides version of my nijmegen code,
in which the Nijm-II 1P1 is updated according to Stoks
(and because the earlier version of the Nijm-II 1P1
produces an (unphysical) bound state).

Otherwise, this code is identical to the earlier one
and, thus, also calculates the Nijm93 and Nijm-I pots
in the old fashion.

I replotted the half-off shel 1P1 momnetum-space potentials:
the new Nijm-II 1P1 is now almost identical to the Reid93
1P1, i.e. not crazy anymore.

All the best,
Ruprecht.





      subroutine nijm2
c
c**** 4/8/97
c**** UPDATED VERSION (see NOTE below) of
c**** The NIJMEGEN POTENTIALS in momentum space
c
c**** interface routine to adjust the NIJMEGEN POTENTIAL CODE
c**** (see attached) to the Bonn application routines.
c
c**** this package is selfcontained.
c
c     R. Machleidt;
c     April 8, 1997
c
c     NOTE: this package contains the UPDATED Nijm-II potential
c     in which the 1P1 is changed as compared to the version of 1994;
c     this updated Nijm-II 1P1 was received from Dirk Hueber on 4/8/97;
c     he got it from V. Stoks (who talks to Hueber but not to me).
c     otherwise this routine is identical to the earlier routine nijm.f;
c     that is, this routine also contains the Nijm93 and Nijm-I
c     potentials in their original version of 1994.
c
c
c
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      common/cpot/ v(6),xmev,ymev
      common/cstate/ j,heform,sing,trip,coup,endep,label
      COMMON/EMANHP/PHNAME
      COMMON/CHOICE/IDPAR
      COMMON/RELKIN/NONREL
      logical heform,sing,trip,coup,endep
      CHARACTER PHNAME*3, TYPE*2
      LOGICAL NONREL
      REAL*8 VPOT(2,2)
      dimension vv(6)
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)                         
      character*4 name(3),nname(15)
      data pi/3.141592653589793d0/
      logical index
      data index/.false./
      data jj/-1/
c
c
      if (index) go to 50
      index=.true.
c
c
c        read in parameters for Nijmegen potential
c
      write (kwrite,10000)
10000 format (///' Updated Nijmegen Potential',
     1' (Nijm-II 1P1 updated, 4/8/97)'/
     2           ' --------------------------')
      read (kread,10001) nname
10001 format (15a4)
      write (kwrite,10002) nname
10002 format (' ',15a4)
c        IDPAR = 0, 1, 2 for Nijm 93, I, II(local), respectively.
      read (kread,10003) name,IDPAR
10003 format (2a4,a2,i1)
      write (kwrite,10004) name,IDPAR
10004 format (' ',2a4,a2,i1)
c        TYPE = 'PP', 'NN', 'NP', or 'PN'.
      read (kread,10005) name,TYPE
10005 format (2a4,a2,a2)
      write (kwrite,10006) name,TYPE
10006 format (' ',2a4,a2,a2)
      read (kread,10007) name,NONREL
10007 format (2a4,a2,l1)
      write (kwrite,10008) name,NONREL
10008 format (' ',2a4,a2,l1)
      read (kread,10009) name,label 
10009 format (2a4,a2,a4)
      write (kwrite,10010) name,label 
10010 format (' ',2a4,a2,a4///)
c
c
      endep=.false.
      fa=1./(2.d0*pi*pi)
c
c
   50 if (j.eq.0.or.j.eq.jj) go to 90
      jj=j
c
c
      if (j.gt.9) write (kwrite,19001)
19001 format (///' warning. nijmegen potential for j greater 9'/
     1' not defined.'
     2' the potential is set to zero.'/
     3' execution continued.'///)
c
c
c
      aj=dfloat(j)                                               
      aj1=dfloat(j+1)                                            
      a2j1=dfloat(2*j+1)                                                
      aaj6=dsqrt(aj*aj1)                                                  
c                                                                       
c        coefficient matrix for the translations into lsj formalism
c                                                                       
      adminv(1,1)=aj1                                                   
      adminv(1,2)=aj                                                    
      adminv(1,3)=-aaj6                                                 
      adminv(1,4)=-aaj6                                                 
      adminv(2,1)=aj                                                    
      adminv(2,2)=aj1                                                   
      adminv(2,3)=aaj6                                                 
      adminv(2,4)=aaj6                                                  
      adminv(3,1)=aaj6                                                 
      adminv(3,2)=-aaj6                                                  
      adminv(3,3)=aj1                                                 
      adminv(3,4)=-aj                                                   
      adminv(4,1)=aaj6                                                 
      adminv(4,2)=-aaj6                                                  
      adminv(4,3)=-aj                                                   
      adminv(4,4)=aj1                                                 
c                                                                       
c       inversion                                                     
c                                                                       
      call dminv (adminv,4,deter,ldminv,mdminv)                         
c                                                                       
c
c
c
   90 do 95 iv=1,6
   95 vv(iv)=0.d0
c
c
      if (j.gt.9) go to 2000
c
c
      j1=j+1
      do 295 i=1,3
      if (i.eq.1.and..not.sing) go to 295
      if (i.eq.2.and..not.trip) go to 295
      if (i.eq.3.and..not.coup) go to 295
c
c
      go to (100,110,120,130,140,150,160,170,180,190),j1
c
c
c        j = 0
c
  100 go to (101,102,295),i
  101 phname='1S0'
      go to 200
  102 phname='3P0'
      go to 200
c
c
c        j = 1
c
  110 go to (111,112,113),i
  111 phname='1P1'
      go to 200
  112 phname='3P1'
      go to 200
  113 phname='3C1'
      go to 200
c
c
c        j = 2
c
  120 go to (121,122,123),i
  121 phname='1D2'
      go to 200
  122 phname='3D2'
      go to 200
  123 phname='3C2'
      go to 200
c
c
c        j = 3
c
  130 go to (131,132,133),i
  131 phname='1F3'
      go to 200
  132 phname='3F3'
      go to 200
  133 phname='3C3'
      go to 200
c
c
c        j = 4
c
  140 go to (141,142,143),i
  141 phname='1G4'
      go to 200
  142 phname='3G4'
      go to 200
  143 phname='3C4'
      go to 200
c
c
c        j = 5
c
  150 go to (151,152,153),i
  151 phname='1H5'
      go to 200
  152 phname='3H5'
      go to 200
  153 phname='3C5'
      go to 200
c
c
c        j = 6
c
  160 go to (161,162,163),i
  161 phname='1I6'
      go to 200
  162 phname='3I6'
      go to 200
  163 phname='3C6'
      go to 200
c
c
c        j = 7
c
  170 go to (171,172,173),i
  171 phname='1J7'
      go to 200
  172 phname='3J7'
      go to 200
  173 phname='3C7'
      go to 200
c
c
c        j = 8
c
  180 go to (181,182,183),i
  181 phname='1K8'
      go to 200
  182 phname='3K8'
      go to 200
  183 phname='3C8'
      go to 200
c
c
c        j = 9
c
  190 go to (191,192,193),i
  191 phname='1L9'
      go to 200
  192 phname='3L9'
      go to 200
  193 phname='3C9'
      go to 200
c
c
c
c
  200 call pnymlsj (ymev,xmev,type,vpot)
c
c
c
c
      go to (201,202,203),i
  201 vv(1)=vpot(1,1)
      go to 295
  202 vv(2)=vpot(1,1)
      go to 295
  203 vv(3)=vpot(2,2)
      vv(4)=vpot(1,1)
      vv(5)=vpot(2,1)
      vv(6)=vpot(1,2)
  295 continue
c
c
      if (j.ne.0) go to 1000
      vv(3)=vv(2)
      vv(2)=0.d0
      go to 2000
c
c
c
 1000 if (.not.heform) go to 2000                                      
c                                                                       
c                                                                       
c         translation into (combination of) helicity states
c
c
      do 1005 i=1,4
 1005 vl(i)=vv(i+2)
c
      do 1020 ii=1,4
      iii=ii+2
      vv(iii)=0.d0
c
      do 1015 i=1,4                                                      
 1015 vv(iii)=vv(iii)+adminv(ii,i)*vl(i)                                  
 1020 vv(iii)=vv(iii)*a2j1                                      
c
c                                                                       
c
c
c
c         over-all factors
c
 2000 do 2005 iv=1,6
 2005 v(iv)=vv(iv)*fa
c
c
      return
      end
      SUBROUTINE PNYMLSJ(QI,QF,TYPE,VPOT)
c
c       this is the original code as obtained by anonymous ftp
c       from nijmegen on 10/20/94
c       the only change made: 
c       the write statement that writes out the parameters
c       is commented out.
c
************************************************************************
**    NIJMEGEN NUCLEON-NUCLEON POTENTIAL PROGRAM                      **
**    ------------------------------------------                      **
**    Version: June 1994                                              **
**    E-mail: thefalg@sci.kun.nl                                      **
**    Reference: Stoks et al. Phys.Rev. C49 (1994) June               **
**                                                                    **
**    Refined and extended version of the 1978 Nijmegen NN potential  **
**    in momentum space on LSJ basis.                                 **
**    Basic references for formulas and set up of this program can be **
**    found in Phys.Rev.D17 (1978) 768 and Phys.Rev.C40 (1989) 2226.  **
**    For momentum-space projection see THEF-NYM-91.05 preprint       **
**    (available via anonymous ftp from thef-nym.sci.kun.nl).         **
**                                                                    **
**                                                                    **
**    INPUT :  QI    center of mass momentum initial state in MeV     **
**    -----    QF    center of mass momentum final   state in MeV     **
**             TYPE  'PP', 'NN', 'NP', or 'PN' (character*2)          **
**             Name partial wave via COMMON/EMANHP/PHNAME (see below) **
**             Maximum total angular momentum J=9 !!                  **
**                                                                    **
**    OUTPUT:  central VC, spin-spin VS, tensor VT, spin-orbit VLS,   **
**    ------   asymmetric spin-orbit VLA, quadratic spin-orbit VQ,    **
**             and extra quadratic spin-orbit part VQP.               **
**             All potentials are calculated via a partial-wave       **
**             decomposition from Jmin up to Jmax, communicated via   **
**                COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2), **
**                              VLA(-2:2),VQ(-3:3),VQP(-2:2)          **
**             The subroutine returns a 2x2 potential matrix VPOT     **
**             in MeV**-2 which is the partial-wave momentum-space    **
**             potential for the partial wave PHNAME (see below)      **
**--------------------------------------------------------------------**
**    Defining the K-matrix as   2i*mu*q*K = (1-S)(1+S)^-1            **
**    (so for singlet channel   tan(delta)=-2*mu*q*K )                **
**    the partial-wave Lippmann-Schwinger equation reads              **
**                                                                    **
**       K(q'q) = V(q'q) + 2/pi int dk k^2 V(q'k) G(q,k) K(kq)        **
**    with                                                            **
**       G(q,k) = P / (E(q) - k^2/2/mu)                               **
**       V(q'k) = 1 / (4*pi) * VPOT(QI=k,QF=q')                       **
**--------------------------------------------------------------------**
**    Potential decomposition in momentum space plane-wave basis:     **
**    V(QF,QI) = VC                                                   **
**             + VS   (SIG1.SIG2)                                     **
**             + VT   [(SIG1.K)(SIG2.K)-K2/3(SIG1.SIG2)]              **
**             + VLS  (i/2)(SIG1+SIG2).N                              **
**             + VLA  (i/2)(SIG1-SIG2).N            (NOT USED !!!)    **
**             + VQ   (SIG1.N)(SIG2.N)                                **
**             + VQP  [(SIG1.Q)(SIG2.Q)-Q2(SIG1.SIG2)                 **
**                    -(1/4)(SIG1.K)(SIG2.K)+(1/4)K2(SIG1.SIG2)]      **
**                                                                    **
**          K = QF - QI ,   Q = (QF+QI)/2 ,   N = QI X QF = Q X K     **
**                                                                    **
**    NOTE: In the partial wave decomposition we used the             **
**          SYM-convention.                                           **
**          If you use another convention in your Lippmann-Schwinger  **
**          programm, you may need an extra minus sign for the        **
**          the off-diagonal tensor potential VPOT(1,2) and VPOT(2,1) **
**                                                                    **
**    NQ12  integer which opts the full Fourier transform for the     **
**          quadratic spin-orbit Q12 operator                         **
**       0: Fourier from (SIG1.N)(SIG2.N) --> Q12 ==> approximation   **
**       1: exact Fourier Q12 --> (SIG1.N)(SIG2.N)+extra              **
**          so same potential in coordinate space as in momentum space**
**    ( Presently NQ12=1 is included as a DATA statement )            **
**                                                                    **
**    COMMON-blocks which have to be filled beforehand:               **
**    +  COMMON/CHOICE/IDPAR                                          **
**              IDPAR is an integer and denotes the various different **
**                    models that can be chosen.                      **
**              IDPAR=0: nijm93: potential for pp and np together.    **
**                               including a phenomenological         **
**                               parameter to give the 1S0 pp and np  **
**                               phase shift/scattering length        **
**                               difference                           **
**                    2: nijmI : Reidlike model, each partial wave has**
**                               its own parameterset                 **
**                    3: nijmII: like nijmI, but fully local          **
**    +  COMMON/EMANHP/PHNAME                                         **
**              PHNAME is character*3 and contains the name of the    **
**              partial wave in the spectral notation.                **
**              - singlets:           1S0  1P1  1D2  1F3  1G4 ...     **
**              - triplets uncoupled: 3P0  3P1  3D2  3F3  3G4 ...     **
**              - triplets coupled:        3C1  3C2  3C3  3C4 ...     **
**                where 3C1 denotes  3S1 -- EPS1 -- 3D1 channel       **
**                      3C2 denotes  3P2 -- EPS2 -- 3F2 channel ...   **
**    +  COMMON/RELKIN/NONREL                                         **
**              NONREL is a logical used in the IDPAR=1 and 2 options.**
**              NONREL=.TRUE.  gives the deuteron binding energy of   **
**                             B=2.224575 MeV using non-relativistic  **
**                             kinematics.                            **
**              NONREL=.FALSE. gives the deuteron binding energy of   **
**                             B=2.224575 MeV using relativistic      **
**                             kinematics.                            **
**              Model IDPAR=0 gives the deuteron only using           **
**              relativistic kinematics ==> NONREL=.FALSE.            **
**                                                                    **
**    NOTE: ALL potential models use a fixed fpi**2=0.075 for the     **
**    ----- pion-nucleon coupling constant at the pion pole, which    **
**          is represented by the DATA FPPPI0/0.075D0/ statement.     **
**                                                                    **
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VPOT(2,2)
      INTEGER SPIN
      CHARACTER PHNAME*3, PHNAM0*3, TYPE*2
      COMMON/EMANHP/PHNAME
      COMMON/CHOICE/IDPAR
      COMMON/JFAC/  JMMM,JMM,JM,J,JP,JPP,JPPP
      COMMON/QFAC/  QI2,QF2,QIQF,QI2F2,QI2F22,QIPF2
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      DATA HBARC/197.327053D0/
      DATA NQ12/1/, PHNAM0/'***'/,  PI/3.14159265358979323846D0/
      SAVE NCHAN,SPIN,L,ISO,  TJMM,TJM,TJ,TJP,TJPP,TJJ
      SAVE E00,F00,G00,E11,F11,G11,EMM,EPP,GMM,GPP,FMM,FPP,FPM

      IF(DABS(QI-QF).GT.30*HBARC) THEN
C*      Necessary to avoid underflows
        VPOT(1,1) = 0D0
        VPOT(1,2) = 0D0
        VPOT(2,1) = 0D0
        VPOT(2,2) = 0D0
        RETURN
      ENDIF

      NLOC=1
C*    Only IDPAR=2 is fully local model
      IF(IDPAR.EQ.2) NLOC=0

      IF(PHNAME.NE.PHNAM0) THEN
        PHNAM0=PHNAME
        NCHAN=1
        IF(PHNAME(2:2).EQ.'C') NCHAN=2
        IF(PHNAME(1:1).EQ.'1') SPIN=0
        IF(PHNAME(1:1).EQ.'3') SPIN=1
        READ(PHNAME,'(2X,I1)') J
        IF(J.GT.9) WRITE(*,*)
     .             '**** Partial wave exceeds allowable maximum J=9'
        IF(J.GT.9) STOP
        L=J
        IF(PHNAME.EQ.'3P0') L=1
        IF(NCHAN.EQ.2) L=J-1
        ISO=MOD(SPIN+L+1,2)
        JMMM=J-3
        JMM =J-2
        JM  =J-1
        JP  =J+1
        JPP =J+2
        JPPP=J+3
        TJMM=2D0*J-3D0
        TJM =2D0*J-1D0
        TJ  =2D0*J+1D0
        TJP =2D0*J+3D0
        TJPP=2D0*J+5D0
        TJJ =DSQRT(J*(J+1D0))
C**     J-dependent coefficients for quadratic spin-orbit
        E00 = J*JM/TJM/TJ
        F00 =-2D0*(J*J+JM)/TJM/TJP
        G00 = JP*JPP/TJ/TJP
        E11 = JM*JPP/TJM/TJ
        F11 =-2D0*JM*JPP/TJM/TJP
        G11 = JM*JPP/TJ/TJP
        EMM =-JM*JMM/TJM/TJMM
        EPP =-J*(2D0*J*J+7D0*JP)/TJ/TJ/TJP
        GMM =-(2D0*J*J-3D0*J+2D0)*JP/TJM/TJ/TJ
        GPP =-JPP*JPPP/TJP/TJPP
        FMM = 2D0*(2D0*J*J*J-3D0*J*J-2D0*JM)/TJ/TJ/TJMM
        FPM = 2D0*TJJ/TJ/TJ
        FPP = 2D0*(2D0*J*J*J+9D0*J*J+10D0*J+1D0)/TJ/TJ/TJPP
      ENDIF

      QI2=QI*QI
      QF2=QF*QF
      QIQF=QI*QF
      QI2F2=QI2+QF2
      QI2F22=QI2F2*QI2F2
      QIPF2=(QI+QF)*(QI+QF)
      S2PSI=2D0*QIQF/QI2F2
      SPSI2=QF2/QI2F2
      CPSI2=QI2/QI2F2

      CALL VMOM(TYPE,NLOC,NQ12,ISO)

      IF(NCHAN.EQ.1) THEN
        IF(SPIN.EQ.0) THEN
          VPOT(1,1) = VC(0) - 3D0*VS(0) +
     .                QI2*QF2*(E00*VQ(-2)+F00*VQ(0)+G00*VQ(2))
        ELSEIF(L.EQ.J) THEN
          VPOT(1,1) = VC(0) + VS(0) - QIQF*(VLS(-1)-VLS(1))/TJ +
     .                2D0/3D0*QI2F2*
     .                (VT(0)-0.5D0*S2PSI*(TJP*VT(-1)+TJM*VT(1))/TJ) +
     .                QI2*QF2*(E11*VQ(-2)+F11*VQ(0)+G11*VQ(2))
        ELSEIF(PHNAME.EQ.'3P0') THEN
          VPOT(1,1) = VC(1) + VS(1) - QIQF*JPP/TJP*(VLS(0)-VLS(2)) +
     .                2D0/3D0*QI2F2*JPP/TJ*
     .                (-VT(1)+0.5D0*S2PSI*(TJPP*VT(0)+TJ*VT(2))/TJP) +
     .                QI2*QF2*(EPP*VQ(-1) + FPP*VQ(1) + GPP*VQ(3))
        ENDIF
      ELSE
        VPOT(1,1) = VC(-1) + VS(-1) + QIQF*JM/TJM*(VLS(-2)-VLS(0)) +
     .              2D0/3D0*QI2F2*JM/TJ*
     .              (-VT(-1)+0.5D0*S2PSI*(TJMM*VT(0)+TJ*VT(-2))/TJM) +
     .              QI2*QF2*(EMM*VQ(-3) + FMM*VQ(-1) + GMM*VQ(1))
        VPOT(1,2) =-2D0*QI2F2*TJJ/TJ*
     .              (-S2PSI*VT(0)+CPSI2*VT(-1)+SPSI2*VT(1)) -
     .              QI2*QF2*FPM*(VQ(1)-VQ(-1))
        VPOT(2,1) =-2D0*QI2F2*TJJ/TJ*
     .              (-S2PSI*VT(0)+SPSI2*VT(-1)+CPSI2*VT(1)) -
     .              QI2*QF2*FPM*(VQ(1)-VQ(-1))
        VPOT(2,2) = VC(1) + VS(1) - QIQF*JPP/TJP*(VLS(0)-VLS(2)) +
     .              2D0/3D0*QI2F2*JPP/TJ*
     .              (-VT(1)+0.5D0*S2PSI*(TJPP*VT(0)+TJ*VT(2))/TJP) +
     .              QI2*QF2*(EPP*VQ(-1) + FPP*VQ(1) + GPP*VQ(3))
      ENDIF

      IF(NQ12.EQ.1) THEN
C**     Extra contribution from inverse Fourier transform of Q12
C**     so momentum space exactly equivalent to configuration space
        IF(NCHAN.EQ.1) THEN
          IF(SPIN.EQ.0) THEN
            VPOT(1,1) = VPOT(1,1) - 2D0*QIQF/TJ*(J*VQP(-1)+JP*VQP(1))
          ELSEIF(L.EQ.J) THEN
            VPOT(1,1) = VPOT(1,1) + QIQF*(-VQP(-1)+VQP(1))/TJ
          ELSEIF(PHNAME.EQ.'3P0') THEN
            VPOT(1,1) = VPOT(1,1) + QIQF*
     .                  ((2D0*J*J+5D0*J+4D0)/TJ*VQP(0)+JPP*VQP(2))/TJP
          ENDIF
        ELSE
          VPOT(1,1) = VPOT(1,1) + QIQF*
     .                (JM*VQP(-2)+(2D0*J*J-JM)/TJ*VQP(0))/TJM
          VPOT(1,2) = VPOT(1,2) + 2D0*QIQF*TJJ/TJ*VQP(0)
          VPOT(2,1) = VPOT(2,1) + 2D0*QIQF*TJJ/TJ*VQP(0)
          VPOT(2,2) = VPOT(2,2) + QIQF*
     .                ((2D0*J*J+5D0*J+4D0)/TJ*VQP(0)+JPP*VQP(2))/TJP
        ENDIF
      ENDIF

      RETURN
      END
************************************************************************
      SUBROUTINE VMOM(TYPE,NLOC,NQ12,ISO)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NEUTM
      CHARACTER PHNAME*3,PHNAM0*3, TYPE*2,TYP0*2
      LOGICAL FIRST
      COMMON/EMANHP/PHNAME
      COMMON/CHOICE/IDPAR
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      COMMON/MESONM/AMPI,AMETA,AMETP, AMRO,AMOM,AMFI, AMA0,AMEP,AMF0,
     .              AMPIC,AMROC,AMA0C,AWPIC,AWROC,AWA0C, AVSC
      COMMON/PARAMS/PAR(6,5)
      COMMON/COPLNG/ALPV,THPV,PV1, FPI,FETA,FETP,FPI2,FETA2,FETP2,
     .              ALVE,THV ,GV1, GRO,GOM,GFI,  GRO2,GOM2,GFI2,
     .              ALVM,     FV1, FRO,FOM,FFI,  FRO2,FOM2,FFI2,
     .              ALGS,THGS,GS1, GA0,GEP,GF0,  GA02,GEP2,GF02,
     .              GFPRO,GFPOM,GFPFI, GFMRO,GFMOM,GFMFI,
     .              FPIC2,GROC2,FROC2,GFPROC,GFMROC,GA0C2
      COMMON/YUKEFF/ARO,AMR1,AROC,AMRC1,AWRC1,BRO,AMR2,BROC,AMRC2,AWRC2,
     .              AVSC1,AVSC2, AEPS,AME1,BEPS,AME2
      COMMON/BROADM/GAMRO,THRRO,GAMRC,THRRC,GAMEP,THRE0,THREC,NRO,NEP
      COMMON/SCALIN/AMT,AMPV
      COMMON/CUTOFF/ALAM,ALAMP,ALAMV,ALAMS
      COMMON/POMRON/GPOM,FPOM2,AMPOM,AMPOM2,AMPOM4,
     .              GA2D,FA2D2,AMA2D,AMA2D2,AMA2D4
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2
      DATA FPPPI0/0.075D0/, FIRST/.FALSE./
      DATA TYP0/'XX'/, PHNAM0/'***'/
      DATA PI/3.14159265358979323846D0/
      SAVE ISIGN,CONV,PROTM,NEUTM,HBARC

      IF(FIRST) GOTO 10
      FIRST = .TRUE.
      CONV  = PI/180D0
**    Nucleon and meson masses (Particle Data Group 1990)
      HBARC = 197.327053D0
      PROTM = 938.27231D0
      NEUTM = 939.56563D0
      AMT   = 938.27231D0
      AMPV  = 139.5675D0
      AMPIC = 139.5675D0
      AMPI  = 134.9739D0
      AMETA = 548.8D0
      AMETP = 957.5D0
      AMROC = 768.3D0
      AMRO  = 768.7D0
      AMOM  = 781.95D0
      AMFI  =1019.412D0
      AMA0C = 983.3D0
      AMA0  = 983.3D0
      AMF0  = 975.6D0
      AVSC  = ((NEUTM-PROTM)/AMROC)**2
      FAC   = 2D0*DSQRT(PROTM*NEUTM)/(PROTM+NEUTM)
      AWPIC = FAC*DSQRT(AMPIC**2 - (NEUTM-PROTM)**2)
      AWROC = FAC*DSQRT(AMROC**2 - (NEUTM-PROTM)**2)
      AWA0C = FAC*DSQRT(AMA0C**2 - (NEUTM-PROTM)**2)

**    Broad rho-meson: spectral density to two effective Yukawa's
**    Yukawa's fitted to PHI0C 0.001 - 2 fm (steps 0.002, ALAMV=825 MeV)
      GAMRO = 152.4D0
      THRRO = AMPIC+AMPIC
      GAMRC = 149.1D0
      THRRC = AMPIC+AMPI
      NRO  = 1
      ARO  = 0.2655205D0
      BRO  = 0.5607493D0
      AMR1 = 645.3772D0
      AMR2 = 878.3667D0
      AROC = 0.3875515D0
      BROC = 0.4508341D0
      AMRC1= 674.1521D0
      AMRC2= 929.9742D0
      AWRC1 = FAC*DSQRT(AMRC1**2 - (NEUTM-PROTM)**2)
      AWRC2 = FAC*DSQRT(AMRC2**2 - (NEUTM-PROTM)**2)
      AVSC1 = ((NEUTM-PROTM)/AMRC1)**2
      AVSC2 = ((NEUTM-PROTM)/AMRC2)**2

   10 CONTINUE
      IF(TYPE.EQ.TYP0) THEN
        IF((IDPAR.EQ.1 .OR. IDPAR.EQ.2) .AND. PHNAME.NE.PHNAM0) GOTO 15
        GOTO 20
      ENDIF
      TYP0=TYPE
      IF(TYPE.EQ.'PP') THEN
        AMY = PROTM
        AMN = PROTM
        I3Y = 1
        I3N = 1
      ELSEIF(TYPE.EQ.'NP') THEN
        AMY = NEUTM
        AMN = PROTM
        I3Y =-1
        I3N = 1
      ELSEIF(TYPE.EQ.'PN') THEN
        AMY = PROTM
        AMN = NEUTM
        I3Y = 1
        I3N =-1
      ELSEIF(TYPE.EQ.'NN') THEN
        AMY = NEUTM
        AMN = NEUTM
        I3Y =-1
        I3N =-1
      ENDIF
      ISIGN = I3Y*I3N
      CALL AMFACS

   15 PHNAM0=PHNAME
      CALL NYMPAR(ISIGN)
C              pseudovec  vector    tensor   scalar    pomeron
C              FPI        GRO       FRO      GA0       GPOM
C              PV1        GV1       FOM      GEP       GA2D
C              ALPV       ALVE      FFI      GF0       AMPOM
C              THPV       THV       RHOC%    THGS      AMA2D
C     Cut-off  ALAMP      ALAMV     ALAM     ALAMS     AMEP
C     2 Yukawa AEPS       AME1      BEPS     AME2      GAMEP

**    (joined) cutoffs for pseudoscalar, vector, scalar
      ALAMP = PAR(5,1)
      ALAMV = PAR(5,2)
      ALAMS = PAR(5,4)
      ALAM  = PAR(5,3)
      IF(ALAMP.EQ.0D0) ALAMP=ALAM
      IF(ALAMV.EQ.0D0) ALAMV=ALAM
      IF(ALAMS.EQ.0D0) ALAMS=ALAM
**    pseudovector couplings
C     FPI   = PAR(1,1)
      FPI   = DSQRT(FPPPI0*FDEXP(-(AMPI/ALAMP)**2))
      PV1   = PAR(2,1)
      ALPV  = PAR(3,1)
      THPV  = PAR(4,1)*CONV
**    vector couplings
      GRO   = PAR(1,2)
      GV1   = PAR(2,2)
      ALVE  = PAR(3,2)
      THV   = PAR(4,2)*CONV
**    tensor couplings
      FRO   = PAR(1,3)
      FOM   = PAR(2,3)
      FFI   = PAR(3,3)
**    scalar couplings
      GA0   = PAR(1,4)
      GEP   = PAR(2,4)
      GF0   = PAR(3,4)
      THGS  = PAR(4,4)*CONV
**    diffractive contribution
      GPOM  = PAR(1,5)
      GA2D  = PAR(2,5)
      AMPOM = PAR(3,5)
      AMA2D = PAR(4,5)
      IF(AMA2D.EQ.0D0) AMA2D=AMPOM
      AMPOM2= AMPOM*AMPOM
      AMPOM4= AMPOM2*AMPOM2
      AMA2D2= AMA2D*AMA2D
      AMA2D4= AMA2D2*AMA2D2

      CALL NYMCOP(ISIGN)

**    Broad epsilon-meson: spectral density to two effective Yukawa's
**    Yukawa's fitted to PHI0C from 0.001 to 2 fm (steps 0.005)
      AMEP  = PAR(5,5)
      GAMEP = PAR(6,5)
      THRE0 = AMPI+AMPI
      THREC = AMPIC+AMPIC
      NEP  = 0
      AEPS = PAR(6,1)
      AME1 = PAR(6,2)
      BEPS = PAR(6,3)
      AME2 = PAR(6,4)

**    isospin dependence
   20 FACISO= 4D0*ISO-2D0
      GA0C2 = FACISO*GA0*GA0
      FPIC2 = FACISO*FPI*FPI/(AMPV*AMPV)
      GROC2 = FACISO*GRO*GRO
      FROC2 = FACISO*FRO*FRO/(4D0*AMT*AMT)
      GFPROC= FACISO*(GRO*FRO+FRO*GRO)/(2D0*AMT)
      GFMROC= FACISO*(GRO*FRO-FRO*GRO)/(2D0*AMT)
      IF(PHNAME.EQ.'1S0' .AND. (TYPE.EQ.'NP'.OR.TYPE.EQ.'PN')) THEN
**      Phenomenological phase-shift difference in PP and NP 1S0
        GROC2=GROC2*(1D0+PAR(4,3))
        FROC2=FROC2*(1D0+PAR(4,3))
      ENDIF

      FPOM2 = (GPOM/AMT)**2
      FA2D2 = (GA2D/AMT)**2

      DO 50 KK=-1,1
        VC(KK) = 0D0
   50   VS(KK) = 0D0
      DO 51 KK=-2,2
        VT(KK)  = 0D0
        VLS(KK) = 0D0
        VLA(KK) = 0D0
        VQ(KK)  = 0D0
   51   VQP(KK) = 0D0
      VQ(-3) = 0D0
      VQ( 3) = 0D0

      CALL PSVECT(TYPE,NLOC)
      CALL VECTOR(TYPE,NLOC,NQ12)
      CALL SCALAR(TYPE,NLOC,NQ12)
      CALL DIFRAC(TYPE,NLOC,NQ12,ISIGN,FACISO)

      RETURN
      END
************************************************************************
      SUBROUTINE AMFACS
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2

      AMY2  = AMY*AMY
      AMYI  = 1D0/AMY
      AMYI2 = AMYI*AMYI
      AMN2  = AMN*AMN
      AMNI  = 1D0/AMN
      AMNI2 = AMNI*AMNI
      AMYPN = AMY+AMN
      AMYMN = AMY-AMN
      AMYPN2= AMYPN*AMYPN
      AMYMN2= AMYMN*AMYMN
      AY2PN2= AMY2+AMN2
      AY2MN2= AMY2-AMN2
      AMYN  = AMY*AMN
      AMYNI = AMYI*AMNI
      AMYNI2= AMYNI*AMYNI
      AYPNI2= AMYI2+AMNI2
      AYMNI2= AMYI2-AMNI2
      REDM = AMYN/AMYPN
      ALF  = 4D0*REDM/AMYPN
      IF(AMY.EQ.AMN) ALF = 1D0
      RETURN
      END
************************************************************************
      SUBROUTINE NYMPAR(ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAMS/PAR(6,5)
      COMMON/CHOICE/IDPAR
      DIMENSION PARTR(5,6)
C              pseudovec  vector    tensor   scalar    pomeron
C              FPI        GRO       FRO      GA0       GPOM
C              PV1        GV1       FOM      GEP       GA2D
C              ALPV       ALVE      FFI      GF0       AMPOM
C              THPV       THV       RHOC%    THGS      AMA2D
C     Cut-off  ALAMP      ALAMV     ALAM     ALAMS     AMEP
C     2 Yukawa AEPS       AME1      BEPS     AME2      GAMEP
C----------------------------------------------------------------------
 
**    Parameter-set with separate cutoffs, fitted to all data    Nijm93
      DATA PARTR/
     1 .2720668D+00,.9209319D+00,.3770582D+01,.1384689D+01,.5228672D+01
     2,.1595311D+00,.2594356D+01,.5816365D+00,.5310001D+01,.2204600D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.3484505D+01,.2081618D+03
     4,-.230000D+02,.3750000D+02,.4371144D-01,.3790000D+02,.0000000D+00
     5,.1177107D+04,.9045040D+03,.0000000D+00,.5544013D+03,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
 
      IF(IDPAR.EQ.0) THEN
        DO 1 I=1,6
          DO 1 J=1,5
  1         PAR(I,J)=PARTR(J,I)
      ELSEIF(IDPAR.EQ.1) THEN
        CALL PHSRDL(ISIGN)
      ELSEIF(IDPAR.EQ.2) THEN
        CALL PHSLOC(ISIGN)
      ENDIF
 
c**** WRITE(*,*) ' NIJMEGEN POTENTIAL PARAMETERS ARE:'
c**** DO 4 I=1,6
c*4     WRITE(*,5) (PAR(I,J),J=1,5)
  5   FORMAT(5(1X,D15.7))
      RETURN
      END
************************************************************************
C     This subroutine reads the parameters of the 0-350 MeV fitted
C     Reidlike potential, Jan 1993, chi**2/data=1.03             Nijm I
C-----------------------------------------------------------------------
      SUBROUTINE PHSRDL(ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER PHNAME*3
      LOGICAL NONREL
      COMMON/PARAMS/PAR(6,5)
      COMMON/EMANHP/PHNAME
      COMMON/RELKIN/NONREL
      DIMENSION APP1S0(5,6),ANP1S0(5,6),PAR1P1(5,6),PAR1D2(5,6),
     .          PAR1F3(5,6),PAR1G4(5,6),PAR3P1(5,6),PAR3D2(5,6),
     .          PAR3F3(5,6),PAR3G4(5,6),PAR3P0(5,6),PAR3C1(5,6),
     .          PAR3C2(5,6),PAR3C3(5,6),PAR3C4(5,6),PARRST(5,6),
     .          PNR3C1(5,6)
      DATA APP1S0/
     1 .2702427D+00,.6723103D+00,.3697581D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4846532D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4744299D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA ANP1S0/
     1 .2702427D+00,.6723103D+00,.3423514D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4942984D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1P1/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.3290519D+01,.5473693D-01,.4784616D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1D2/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.0000000D+00,.5473693D-01,.2348432D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.3639503D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1F3/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.4812761D+01,.5473693D-01,.6316315D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1G4/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.4142482D+01,.5473693D-01,.4859055D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3P1/
     1 .2702427D+00,.6723103D+00,.0000000D+00,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4980991D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4923498D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3D2/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5880449D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.1693534D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3F3/
     1 .2702427D+00,.6723103D+00,.8961465D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5458122D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3G4/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.3201615D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3P0/
     1 .2702427D+00,.6723103D+00,.3160965D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.3726932D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C1/
     1 .2706937D+00,.2598577D+01,.2125742D+01,.8334762D+00,.2751792D+01
     2,.0000000D+00,.2536422D+01,.5473693D-01,.4906761D+01,.44371888D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8848614D+03,.6465221D+03,.0000000D+00,.6990612D+03,.7600000D+03
     6,.1690008D+00,.5831763D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C2/
     1 .2702427D+00,.6723103D+00,.5188648D+01,.8334762D+00,.3668014D+01
     2,.1341669D+00,.0000000D+00,.5473693D-01,.3995761D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8275346D+03,.0000000D+00,.5831699D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C3/
     1 .2702427D+00,.6723103D+00,.3877361D+01,.8334762D+00,.4723766D+01
     2,.2811080D+00,.0000000D+00,.5473693D-01,.4818122D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C4/
     1 .2702427D+00,.6723103D+00,.7551377D+01,.8334762D+00,.2360920D+01
     2,.2871299D+00,.0000000D+00,.5473693D-01,.4855317D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PARRST/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4859055D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
C     Parameters for a nonlocal potential with non-relativistic deuteron
      DATA PNR3C1/
     1 .2706937D+00,.2588931D+01,.2120162D+01,.8334762D+00,.2751792D+01
     2,.0000000D+00,.2521404D+01,.5473693D-01,.4889314D+01,.44371806D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8848614D+03,.6465221D+03,.0000000D+00,.6990612D+03,.7600000D+03
     6,.1690008D+00,.5831292D+03,.6130152D+00,.1021139D+04,.6400000D+03/
 
      DO 1 I=1,6
        DO 1 J=1,5
          IF(PHNAME.EQ.'1S0') THEN
            IF(ISIGN.EQ. 1) PAR(I,J)=APP1S0(J,I)
            IF(ISIGN.EQ.-1) PAR(I,J)=ANP1S0(J,I)
          ELSEIF(PHNAME.EQ.'1P1') THEN
            PAR(I,J)=PAR1P1(J,I)
          ELSEIF(PHNAME.EQ.'1D2') THEN
            PAR(I,J)=PAR1D2(J,I)
          ELSEIF(PHNAME.EQ.'1F3') THEN
            PAR(I,J)=PAR1F3(J,I)
          ELSEIF(PHNAME.EQ.'1G4') THEN
            PAR(I,J)=PAR1G4(J,I)
          ELSEIF(PHNAME.EQ.'3P0') THEN
            PAR(I,J)=PAR3P0(J,I)
          ELSEIF(PHNAME.EQ.'3P1') THEN
            PAR(I,J)=PAR3P1(J,I)
          ELSEIF(PHNAME.EQ.'3D2') THEN
            PAR(I,J)=PAR3D2(J,I)
          ELSEIF(PHNAME.EQ.'3F3') THEN
            PAR(I,J)=PAR3F3(J,I)
          ELSEIF(PHNAME.EQ.'3G4') THEN
            PAR(I,J)=PAR3G4(J,I)
          ELSEIF(PHNAME.EQ.'3C1') THEN
            IF(.NOT.NONREL) PAR(I,J)=PAR3C1(J,I)
            IF(NONREL) PAR(I,J)=PNR3C1(J,I)
          ELSEIF(PHNAME.EQ.'3C2') THEN
            PAR(I,J)=PAR3C2(J,I)
          ELSEIF(PHNAME.EQ.'3C3') THEN
            PAR(I,J)=PAR3C3(J,I)
          ELSEIF(PHNAME.EQ.'3C4') THEN
            PAR(I,J)=PAR3C4(J,I)
          ELSE
            PAR(I,J)=PARRST(J,I)
          ENDIF
    1 CONTINUE
      RETURN
      END
************************************************************************
C     This subroutine reads the parameters of the 0-350 MeV fitted local
C     potential (without Q**2), Jan 1993, chi**2/data=1.03       Nijm II
C-----------------------------------------------------------------------
      SUBROUTINE PHSLOC(ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER PHNAME*3
      LOGICAL NONREL
      COMMON/PARAMS/PAR(6,5)
      COMMON/EMANHP/PHNAME
      COMMON/RELKIN/NONREL
      DIMENSION APP1S0(5,6),ANP1S0(5,6),PAR1P1(5,6),PAR1D2(5,6),
     .          PAR1F3(5,6),PAR1G4(5,6),PAR3P1(5,6),PAR3D2(5,6),
     .          PAR3F3(5,6),PAR3G4(5,6),PAR3P0(5,6),PAR3C1(5,6),
     .          PAR3C2(5,6),PAR3C3(5,6),PAR3C4(5,6),PARRST(5,6),
     .          PNR3C1(5,6)
      DATA APP1S0/
     1 .2702427D+00,.6723103D+00,.4006116D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4669822D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.3909783D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4692497D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA ANP1S0/
     1 .2702427D+00,.6723103D+00,.7407998D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4626459D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2196159D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4722381D+03,.6130152D+00,.1021139D+04,.6400000D+03/
c
c******** below is the revised 1p1, as received from D. Hueber on 4/8/97
c******** who got it from V. Stoks.
      DATA PAR1P1/
     1 .2702427D+00,.6723103D+00,.0000000D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.3430827D+01,.5473693D-01,.2825877D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+03,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
cold  DATA PAR1P1/
cold 1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
cold 2,.2871299D+00,.2849628D+01,.5473693D-01,.5411714D+01,.4437220D+00
cold 3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
cold 4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
cold 5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
cold 6,.1690008D+00,.5217623D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1D2/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4138573D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2243917D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4367031D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1F3/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.3983000D+01,.5473693D-01,.5627977D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR1G4/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4859055D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2037620D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3P1/
     1 .2702427D+00,.6723103D+00,.0000000D+00,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4171550D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.3368384D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4530824D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3D2/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5469270D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.1847244D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3F3/
     1 .2702427D+00,.6723103D+00,.6012926D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5530460D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3G4/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.3663270D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3P0/
     1 .2702427D+00,.6723103D+00,.2761025D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.3041218D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8275346D+03,.0000000D+00,.1134832D+04,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C1/
     1 .2702427D+00,.1607944D+01,.1841778D+01,.5469244D+00,.3469472D+01
     2,.2871299D+00,.2240543D+01,.5473693D-01,.4035077D+01,.4437213D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.5151821D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8275346D+03,.0000000D+00,.8044237D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C2/
     1 .2702427D+00,.6723103D+00,.5816373D+01,.8334762D+00,.3957678D+01
     2,.2353573D+00,.0000000D+00,.5473693D-01,.4143714D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2781205D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8275346D+03,.0000000D+00,.6121468D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C3/
     1 .2702427D+00,.6723103D+00,.4050335D+01,.8334762D+00,.4316501D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5048592D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PAR3C4/
     1 .2702427D+00,.6723103D+00,.7347855D+01,.8334762D+00,.2579081D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.5157279D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
      DATA PARRST/
     1 .2702427D+00,.6723103D+00,.4728635D+01,.8334762D+00,.2751792D+01
     2,.2871299D+00,.2849628D+01,.5473693D-01,.4859055D+01,.4437220D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.2579522D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.0000000D+00,.0000000D+00,.8275346D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
C     Parameters for a local potential with non-relativistic deuteron
      DATA PNR3C1/
     1 .2702427D+00,.1710842D+01,.1781765D+01,.5362515D+00,.3461562D+01
     2,.2871299D+00,.2247159D+01,.5473693D-01,.4028595D+01,.44387996D+00
     3,.3550000D+00,.1000000D+01,.0000000D+00,.8389363D+00,.5184792D+03
     4,-.230000D+02,.3750000D+02,.0000000D+00,.3790000D+02,.0000000D+00
     5,.8275346D+03,.0000000D+00,.8044237D+03,.0000000D+00,.7600000D+03
     6,.1690008D+00,.4878179D+03,.6130152D+00,.1021139D+04,.6400000D+03/
 
      DO 1 I=1,6
        DO 1 J=1,5
          IF(PHNAME.EQ.'1S0') THEN
            IF(ISIGN.EQ. 1) PAR(I,J)=APP1S0(J,I)
            IF(ISIGN.EQ.-1) PAR(I,J)=ANP1S0(J,I)
          ELSEIF(PHNAME.EQ.'1P1') THEN
            PAR(I,J)=PAR1P1(J,I)
          ELSEIF(PHNAME.EQ.'1D2') THEN
            PAR(I,J)=PAR1D2(J,I)
          ELSEIF(PHNAME.EQ.'1F3') THEN
            PAR(I,J)=PAR1F3(J,I)
          ELSEIF(PHNAME.EQ.'1G4') THEN
            PAR(I,J)=PAR1G4(J,I)
          ELSEIF(PHNAME.EQ.'3P0') THEN
            PAR(I,J)=PAR3P0(J,I)
          ELSEIF(PHNAME.EQ.'3P1') THEN
            PAR(I,J)=PAR3P1(J,I)
          ELSEIF(PHNAME.EQ.'3D2') THEN
            PAR(I,J)=PAR3D2(J,I)
          ELSEIF(PHNAME.EQ.'3F3') THEN
            PAR(I,J)=PAR3F3(J,I)
          ELSEIF(PHNAME.EQ.'3G4') THEN
            PAR(I,J)=PAR3G4(J,I)
          ELSEIF(PHNAME.EQ.'3C1') THEN
            IF(.NOT.NONREL) PAR(I,J)=PAR3C1(J,I)
            IF(NONREL) PAR(I,J)=PNR3C1(J,I)
          ELSEIF(PHNAME.EQ.'3C2') THEN
            PAR(I,J)=PAR3C2(J,I)
          ELSEIF(PHNAME.EQ.'3C3') THEN
            PAR(I,J)=PAR3C3(J,I)
          ELSEIF(PHNAME.EQ.'3C4') THEN
            PAR(I,J)=PAR3C4(J,I)
          ELSE
            PAR(I,J)=PARRST(J,I)
          ENDIF
    1 CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE NYMCOP(ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/COPLNG/ALPV,THPV,PV1, FPI,FETA,FETP,FPI2,FETA2,FETP2,
     .              ALVE,THV ,GV1, GRO,GOM,GFI,  GRO2,GOM2,GFI2,
     .              ALVM,     FV1, FRO,FOM,FFI,  FRO2,FOM2,FFI2,
     .              ALGS,THGS,GS1, GA0,GEP,GF0,  GA02,GEP2,GF02,
     .              GFPRO,GFPOM,GFPFI, GFMRO,GFMOM,GFMFI,
     .              FPIC2,GROC2,FROC2,GFPROC,GFMROC,GA0C2
      COMMON/SCALIN/AMT,AMPV
      DATA SR3/1.7320508075688772935D0/

**    pseudovector coupling constants
      PV8  = FPI * (4D0*ALPV-1D0)/SR3
      COST = DCOS(THPV)
      SINT = DSIN(THPV)
      FETA = COST*PV8 - SINT*PV1
      FETP = SINT*PV8 + COST*PV1
      FPI2 = FPI*FPI/(AMPV*AMPV) * ISIGN
      FETA2= FETA*FETA/(AMPV*AMPV)
      FETP2= FETP*FETP/(AMPV*AMPV)

**    vector coupling constants
      GV8  = GRO * (4D0*ALVE-1D0)/SR3
      COST = DCOS(THV)
      SINT = DSIN(THV)
      GFI  = COST*GV8 - SINT*GV1
      GOM  = SINT*GV8 + COST*GV1
      GRO2 = GRO*GRO * ISIGN
      GOM2 = GOM*GOM
      GFI2 = GFI*GFI

**    tensor coupling constants
      COST = DCOS(THV)
      SINT = DSIN(THV)
      FV8  = COST*(GFI+FFI) + SINT*(GOM+FOM) - GV8
      FV1  =-SINT*(GFI+FFI) + COST*(GOM+FOM) - GV1
      ALVM = (SR3*(GV8+FV8)+(GRO+FRO))/(4D0*(GRO+FRO))
      FRO2 = FRO*FRO/(4D0*AMT*AMT) * ISIGN
      FOM2 = FOM*FOM/(4D0*AMT*AMT)
      FFI2 = FFI*FFI/(4D0*AMT*AMT)

      GFPRO = (GRO*FRO+FRO*GRO)/(2D0*AMT) * ISIGN
      GFPOM = (GOM*FOM+FOM*GOM)/(2D0*AMT)
      GFPFI = (GFI*FFI+FFI*GFI)/(2D0*AMT)
      GFMRO = (GRO*FRO-FRO*GRO)/(2D0*AMT) * ISIGN
      GFMOM = (GOM*FOM-FOM*GOM)/(2D0*AMT)
      GFMFI = (GFI*FFI-FFI*GFI)/(2D0*AMT)

**    scalar coupling constants
      COST = DCOS(THGS)
      SINT = DSIN(THGS)
      GS8  = COST*GEP + SINT*GF0
      GS1  =-SINT*GEP + COST*GF0
      ALGS = (SR3*GS8+GA0)/(4D0*GA0)
      GA02 = GA0*GA0 * ISIGN
      GEP2 = GEP*GEP
      GF02 = GF0*GF0

      RETURN
      END
************************************************************************
      SUBROUTINE PSVECT(TYPE,NLOC)
C*    Partial-wave momentum-space potentials for pseudoscalar exchange
C*       NLOC=0,1: no nonlocal (q^2+k^2/4) contributions
C*              2: nonlocal contributions in spin-spin and tensor
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER TYPE*2
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      COMMON/JFAC/  JMMM,JMM,JM,J,JP,JPP,JPPP
      COMMON/QFAC/  QI2,QF2,QIQF,QI2F2,QI2F22,QIPF2
      COMMON/MESONM/AMPI,AMETA,AMETP, AMRO,AMOM,AMFI, AMA0,AMEP,AMF0,
     .              AMPIC,AMROC,AMA0C,AWPIC,AWROC,AWA0C, AVSC
      COMMON/COPLNG/ALPV,THPV,PV1, FPI,FETA,FETP,FPI2,FETA2,FETP2,
     .              ALVE,THV ,GV1, GRO,GOM,GFI,  GRO2,GOM2,GFI2,
     .              ALVM,     FV1, FRO,FOM,FFI,  FRO2,FOM2,FFI2,
     .              ALGS,THGS,GS1, GA0,GEP,GF0,  GA02,GEP2,GF02,
     .              GFPRO,GFPOM,GFPFI, GFMRO,GFMOM,GFMFI,
     .              FPIC2,GROC2,FROC2,GFPROC,GFMROC,GA0C2
      COMMON/CUTOFF/ALAM,ALAMP,ALAMV,ALAMS
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2
      DIMENSION U(-3:12), R(-3:12), ELN(15), ERN(15)
      DATA U/16*0D0/, R/16*0D0/
      DATA PI/3.14159265358979323846D0/

      XCOM=0.5D0*QI2F2/QIQF
      ALAM2=ALAMP*ALAMP
      Y=2D0*QIQF/ALAM2
      JMAX=J+2
      KOM=1
      FAC=2D0*PI/QIQF

C**     Pseudoscalar mesons: pi, eta, eta'
      DO 1000 IN=1,3
        IF(IN.EQ.1) THEN
          AMES=AMPI
          FP2 =FPI2
        ELSEIF(IN.EQ.2) THEN
          AMES=AMETA
          FP2 =FETA2
        ELSEIF(IN.EQ.3) THEN
          AMES=AMETP
          FP2 =FETP2
        ENDIF
        AMES2=AMES*AMES
        X=XCOM+0.5D0*AMES2/QIQF
        CALL SEX(X,Y,AMES/ALAMP,JMAX,KOM,ELN,ERN)
        XPS =-FAC*FP2*QI2F2/3D0
        YPS = FAC*FP2*2D0*QIQF/3D0
        XPT =-FAC*FP2
        XNLS= FAC*FP2*AMYNI*QI2F22/12D0
        YNLS=-FAC*FP2*AMYNI*QI2F2*QIQF/6D0
        XNLT= FAC*FP2*AMYNI*QI2F2/4D0
        DO 5 K=0,JMAX
          U(K)=ELN(K+1)
          R(K)=ERN(K+1)
    5   CONTINUE
        VS(-1) = VS(-1) + (XPS+X*YPS)*U(JM) - YPS*R(JM)
        VS( 0) = VS( 0) + (XPS+X*YPS)*U(J)  - YPS*R(J)
        VS( 1) = VS( 1) + (XPS+X*YPS)*U(JP) - YPS*R(JP)
        VT(-2) = VT(-2) +  XPT*U(JMM)
        VT(-1) = VT(-1) +  XPT*U(JM)
        VT( 0) = VT( 0) +  XPT*U(J)
        VT( 1) = VT( 1) +  XPT*U(JP)
        VT( 2) = VT( 2) +  XPT*U(JPP)
        IF(NLOC.EQ.2) THEN
          VS(-1) = VS(-1) + (XNLS+X*YNLS)*U(JM) - YNLS*R(JM)
          VS( 0) = VS( 0) + (XNLS+X*YNLS)*U(J)  - YNLS*R(J)
          VS( 1) = VS( 1) + (XNLS+X*YNLS)*U(JP) - YNLS*R(JP)
          VT(-2) = VT(-2) + XNLT*U(JMM)
          VT(-1) = VT(-1) + XNLT*U(JM)
          VT( 0) = VT( 0) + XNLT*U(J)
          VT( 1) = VT( 1) + XNLT*U(JP)
          VT( 2) = VT( 2) + XNLT*U(JPP)
        ENDIF
 1000 CONTINUE

      IF(TYPE.EQ.'PP' .OR. TYPE.EQ.'NN') RETURN

C**     Charged pseudoscalar: pi+
      AMES=AWPIC
      FP2 =FPIC2
      AMES2=AMES*AMES
      X=XCOM+0.5D0*AMES2/QIQF
      CALL SEX(X,Y,AMES/ALAMP,JMAX,KOM,ELN,ERN)
      XPS =-FAC*FP2/ALF*QI2F2/3D0
      YPS = FAC*FP2/ALF*2D0*QIQF/3D0
      XPT =-FAC*FP2/ALF
      XNLS= FAC*FP2/ALF*AYPNI2*QI2F22/24D0
      YNLS=-FAC*FP2/ALF*AYPNI2*QI2F2*QIQF/12D0
      XNLT= FAC*FP2/ALF*AYPNI2*QI2F2/8D0
      DO 6 K=0,JMAX
        U(K)=ELN(K+1)
        R(K)=ERN(K+1)
    6 CONTINUE
      VS(-1) = VS(-1) + (XPS+X*YPS)*U(JM) - YPS*R(JM)
      VS( 0) = VS( 0) + (XPS+X*YPS)*U(J)  - YPS*R(J)
      VS( 1) = VS( 1) + (XPS+X*YPS)*U(JP) - YPS*R(JP)
      VT(-2) = VT(-2) +  XPT*U(JMM)
      VT(-1) = VT(-1) +  XPT*U(JM)
      VT( 0) = VT( 0) +  XPT*U(J)
      VT( 1) = VT( 1) +  XPT*U(JP)
      VT( 2) = VT( 2) +  XPT*U(JPP)
      IF(NLOC.EQ.2) THEN
        VS(-1) = VS(-1) + (XNLS+X*YNLS)*U(JM) - YNLS*R(JM)
        VS( 0) = VS( 0) + (XNLS+X*YNLS)*U(J)  - YNLS*R(J)
        VS( 1) = VS( 1) + (XNLS+X*YNLS)*U(JP) - YNLS*R(JP)
        VT(-2) = VT(-2) + XNLT*U(JMM)
        VT(-1) = VT(-1) + XNLT*U(JM)
        VT( 0) = VT( 0) + XNLT*U(J)
        VT( 1) = VT( 1) + XNLT*U(JP)
        VT( 2) = VT( 2) + XNLT*U(JPP)
      ENDIF

      RETURN
      END
************************************************************************
      SUBROUTINE VECTOR(TYPE,NLOC,NQ12)
C*    Partial-wave momentum-space potentials for vector exchange
C*       NQ12=0: inexact Fourier of quadratic spin-orbit Q12 operator
C*            1: exact Fourier Q12 from coordinate to momentum
C*       NLOC=0: no nonlocal (q^2+k^2/4) contributions
C*            1: nonlocal contributions in central potential only
C*            2: nonlocal contributions also in spin-spin (only for NP)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER TYPE*2
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      COMMON/JFAC/  JMMM,JMM,JM,J,JP,JPP,JPPP
      COMMON/QFAC/  QI2,QF2,QIQF,QI2F2,QI2F22,QIPF2
      COMMON/MESONM/AMPI,AMETA,AMETP, AMRO,AMOM,AMFI, AMA0,AMEP,AMF0,
     .              AMPIC,AMROC,AMA0C,AWPIC,AWROC,AWA0C, AVSC
      COMMON/COPLNG/ALPV,THPV,PV1, FPI,FETA,FETP,FPI2,FETA2,FETP2,
     .              ALVE,THV ,GV1, GRO,GOM,GFI,  GRO2,GOM2,GFI2,
     .              ALVM,     FV1, FRO,FOM,FFI,  FRO2,FOM2,FFI2,
     .              ALGS,THGS,GS1, GA0,GEP,GF0,  GA02,GEP2,GF02,
     .              GFPRO,GFPOM,GFPFI, GFMRO,GFMOM,GFMFI,
     .              FPIC2,GROC2,FROC2,GFPROC,GFMROC,GA0C2
      COMMON/YUKEFF/ARO,AMR1,AROC,AMRC1,AWRC1,BRO,AMR2,BROC,AMRC2,AWRC2,
     .              AVSC1,AVSC2, AEPS,AME1,BEPS,AME2
      COMMON/CUTOFF/ALAM,ALAMP,ALAMV,ALAMS
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2
      DIMENSION U(-3:12), R(-3:12), S(-3:12), G(-3:12), ELN(15), ERN(15)
      DATA U/16*0D0/, R/16*0D0/, S/16*0D0/, G/16*0D0/
      DATA PI/3.14159265358979323846D0/

      XCOM=0.5D0*QI2F2/QIQF
      ALAM2=ALAMV*ALAMV
      Y=2D0*QIQF/ALAM2
      JMAX=J+3
      KOM=1
      FAC=2D0*PI/QIQF

C**     Vector mesons: broad rho (2 Yukawa's), omega, fi
      DO 1000 IN=1,4
        IF(IN.EQ.1) THEN
          AMES=AMR1
          GV2 =GRO2 *ARO
          FV2 =FRO2 *ARO
          GFV =GFPRO*ARO
          GFM =GFMRO*ARO
        ELSEIF(IN.EQ.2) THEN
          AMES=AMR2
          GV2 =GRO2 *BRO
          FV2 =FRO2 *BRO
          GFV =GFPRO*BRO
          GFM =GFMRO*BRO
        ELSEIF(IN.EQ.3) THEN
          AMES=AMOM
          GV2 =GOM2
          FV2 =FOM2
          GFV =GFPOM
          GFM =GFMOM
        ELSEIF(IN.EQ.4) THEN
          AMES=AMFI
          GV2 =GFI2
          FV2 =FFI2
          GFV =GFPFI
          GFM =GFMFI
        ENDIF
        AMES2=AMES*AMES
        X=XCOM+0.5D0*AMES2/QIQF
        CALL SEX(X,Y,AMES/ALAMV,JMAX+1,KOM,ELN,ERN)
        XVC = FAC * (GV2*(1D0-AMYNI/2D0/ALF*QI2F2) - GFV/4D0/REDM*QI2F2
     .             + FV2*AMYNI/4D0*QI2F22)
        YVC = FAC * (GV2*AMYNI/ALF*QIQF+GFV/2D0/REDM*QIQF
     .             - FV2*AMYNI*QIQF*QI2F2)
        ZVC = FAC * FV2*AMYNI*QI2*QF2
        XVS =-FAC * ((GV2+AMYPN*GFV+4D0*AMYN*FV2)*AMYNI/6D0*QI2F2
     .             - FV2*AMYNI/12D0*QI2F22)
        YVS = FAC * ((GV2+AMYPN*GFV+4D0*AMYN*FV2)*AMYNI/3D0*QIQF
     .             - FV2*AMYNI/3D0*QIQF*QI2F2)
        ZVS = FAC * FV2*AMYNI/3D0*QI2*QF2
        XVT = FAC * ((GV2+AMYPN*GFV+4D0*AMYN*FV2)*AMYNI/4D0
     .             - FV2*AMYNI/8D0*QI2F2)
        YVT = FAC * FV2*AMYNI/4D0*QIQF
        XVLS=-FAC * ((GV2*(0.5D0+1D0/ALF)+AMYPN*GFV)*AMYNI
     .             - FV2*(0.5D0+1D0/ALF)*AMYNI*QI2F2)
        YVLS=-FAC * FV2*(1D0+2D0/ALF)*AMYNI*QIQF
        XVQ =-FAC * (GV2+4D0*AMYPN*GFV+8D0*AMYPN2*FV2)*AMYNI2/16D0
        XVQ2= XVQ * QIQF
        XVNL= FAC * GV2*(AYPNI2+AMYNI)/4D0*QI2F2
        DO 5 K=1,JMAX
          U(K)=ELN(K+1)
          R(K)=ERN(K+1)
          S(K)=(K*ERN(K)+(K+1)*ERN(K+2))/(2*K+1)
          G(K)=(ELN(K+2)-ELN(K))/(2*K+1)
    5   CONTINUE
        U(0)=ELN(1)
        R(0)=ERN(1)
        S(0)=ERN(2)
        ARG=(QIPF2+AMES2)/ALAM2
        G(0)=(ELN(2)-ELN(1))-FDEXP(-QIPF2/ALAM2)*CE1(ARG)
        VC(-1) = VC(-1) + (XVC+X*YVC+X*X*ZVC)*U(JM)
     .                  - (YVC+X*ZVC)*R(JM) - ZVC*S(JM)
        VC( 0) = VC( 0) + (XVC+X*YVC+X*X*ZVC)*U(J)
     .                  - (YVC+X*ZVC)*R(J)  - ZVC*S(J)
        VC( 1) = VC( 1) + (XVC+X*YVC+X*X*ZVC)*U(JP)
     .                  - (YVC+X*ZVC)*R(JP) - ZVC*S(JP)
        VS(-1) = VS(-1) + (XVS+X*YVS+X*X*ZVS)*U(JM)
     .                  - (YVS+X*ZVS)*R(JM) - ZVS*S(JM)
        VS( 0) = VS( 0) + (XVS+X*YVS+X*X*ZVS)*U(J)
     .                  - (YVS+X*ZVS)*R(J)  - ZVS*S(J)
        VS( 1) = VS( 1) + (XVS+X*YVS+X*X*ZVS)*U(JP)
     .                  - (YVS+X*ZVS)*R(JP) - ZVS*S(JP)
        VT(-2) = VT(-2) + (XVT+X*YVT)*U(JMM)- YVT*R(JMM)
        VT(-1) = VT(-1) + (XVT+X*YVT)*U(JM) - YVT*R(JM)
        VT( 0) = VT( 0) + (XVT+X*YVT)*U(J)  - YVT*R(J)
        VT( 1) = VT( 1) + (XVT+X*YVT)*U(JP) - YVT*R(JP)
        VT( 2) = VT( 2) + (XVT+X*YVT)*U(JPP)- YVT*R(JPP)
        VLS(-2)= VLS(-2)+ (XVLS+X*YVLS)*U(JMM)- YVLS*R(JMM)
        VLS(-1)= VLS(-1)+ (XVLS+X*YVLS)*U(JM) - YVLS*R(JM)
        VLS( 0)= VLS( 0)+ (XVLS+X*YVLS)*U(J)  - YVLS*R(J)
        VLS( 1)= VLS( 1)+ (XVLS+X*YVLS)*U(JP) - YVLS*R(JP)
        VLS( 2)= VLS( 2)+ (XVLS+X*YVLS)*U(JPP)- YVLS*R(JPP)
        VQ(-3) = VQ(-3) +  XVQ*U(JMMM)
        VQ(-2) = VQ(-2) +  XVQ*U(JMM)
        VQ(-1) = VQ(-1) +  XVQ*U(JM)
        VQ( 0) = VQ( 0) +  XVQ*U(J)
        VQ( 1) = VQ( 1) +  XVQ*U(JP)
        VQ( 2) = VQ( 2) +  XVQ*U(JPP)
        VQ( 3) = VQ( 3) +  XVQ*U(JPPP)
        IF(NLOC.NE.0) THEN
          VC(-1) = VC(-1) + XVNL*U(JM)
          VC( 0) = VC( 0) + XVNL*U(J)
          VC( 1) = VC( 1) + XVNL*U(JP)
        ENDIF

        IF(NQ12.EQ.1) THEN
C**       Extra contribution from inverse Fourier transform of Q12
          VQP(-2) = VQP(-2) + XVQ2*G(JMM)
          VQP(-1) = VQP(-1) + XVQ2*G(JM)
          VQP( 0) = VQP( 0) + XVQ2*G(J)
          VQP( 1) = VQP( 1) + XVQ2*G(JP)
          VQP( 2) = VQP( 2) + XVQ2*G(JPP)
        ENDIF

 1000 CONTINUE

      IF(TYPE.EQ.'PP' .OR. TYPE.EQ.'NN') RETURN

C**     Charged vector: broad rho+ (2 Yukawa's)
      DO 2000 IN=1,2
        IF(IN.EQ.1) THEN
          AMES=AWRC1
          GV2 =GROC2 *AROC
          FV2 =FROC2 *AROC
          GFV =GFPROC*AROC
          GFM =GFMROC*AROC
          FFAC=FAC*AVSC1
        ELSEIF(IN.EQ.2) THEN
          AMES=AWRC2
          GV2 =GROC2 *BROC
          FV2 =FROC2 *BROC
          GFV =GFPROC*BROC
          GFM =GFMROC*BROC
          FFAC=FAC*AVSC2
        ENDIF
        AMES2=AMES*AMES
        X=XCOM+0.5D0*AMES2/QIQF
        CALL SEX(X,Y,AMES/ALAMV,JMAX+1,KOM,ELN,ERN)
        XVC = FAC * (GV2*(ALF-AMYNI/2D0*QI2F2) - ALF*
     .        (GFV*AMYPN+FV2*AMYMN2)*AYPNI2/8D0*QI2F2 + FV2*AMYNI/4D0*
     .        (1D0+(AY2MN2**2+2D0*AMYN*AMYMN2)*AMYNI2/16D0)*QI2F22)
        YVC = FAC * (GV2*AMYNI*QIQF + ALF*(GFV*AMYPN+FV2*AMYMN2)*
     .        AYPNI2/4D0*QIQF - FV2*AMYNI*
     .        (1D0+(AY2MN2**2+2D0*AMYN*AMYMN2)*AMYNI2/16D0)*QI2F2*QIQF)
        ZVC = FAC * FV2*AMYNI*
     .        (1D0+(AY2MN2**2+2D0*AMYN*AMYMN2)*AMYNI2/16D0)*QI2*QF2
        XVS = FAC * (-(GV2+AMYPN*GFV+AMYPN2*FV2)*
     .        (AMYNI/6D0+ALF*AMYMN2*AMYNI2/16D0)*QI2F2 +
     .        FV2/ALF*AYPNI2/24D0*QI2F22)
        YVS = FAC * ((GV2+AMYPN*GFV+AMYPN2*FV2)*
     .        (AMYNI/3D0+ALF*AMYMN2*AMYNI2/8D0)*QIQF -
     .        FV2/ALF*AYPNI2/6D0*QIQF*QI2F2)
        ZVS = FAC * FV2/ALF*AYPNI2/6D0*QI2*QF2
        XVT = FAC * ((GV2+AMYPN*GFV+AMYPN2*FV2)*AMYNI/4D0 -
     .        FV2/ALF*AYPNI2/16D0*QI2F2)
        YVT = FAC* FV2/ALF*AYPNI2/8D0*QIQF
        XVLS= FAC * (-(GV2*(2D0-ALF/2D0)+2D0*AY2PN2/AMYPN*GFV +
     .               2D0*AMYMN2*FV2)*AMYNI +
     .        FV2*((2D0/ALF-0.5D0)/ALF-AMYMN2*AMYNI/8D0)*AMYNI*QI2F2)
        YVLS=-FAC* FV2*((4D0/ALF-1D0)/ALF-AMYMN2*AMYNI/4D0)*AMYNI*QIQF
        XVQ =-FAC * ALF*(GV2+4D0*AMYPN*GFV+8D0*AMYPN2*FV2)*AMYNI2/16D0
        XVNLC= FAC * (ALF*GV2*(AYPNI2+AMYNI)/2D0+GFV*AMYMN2/AMYPN*AMYNI+
     .                FV2*(4D0/ALF-2D0*ALF*AY2PN2*AMYNI))*QI2F2/2D0
        XVNLS= FAC * ALF*AMYMN2*AMYNI2/8D0*(GV2+AMYPN*GFV+AMYPN2*FV2)*
     .               QI2F2
C**       Second part of vector-vector potential (same as scalar)
        XVC  = XVC  + FFAC * ALF*GV2*(1D0+AMYNI/4D0*QI2F2)
        YVC  = YVC  - FFAC * ALF*GV2*AMYNI/2D0*QIQF
        XVLS = XVLS + FFAC * ALF*GV2*AMYNI/2D0
        XVQ  = XVQ  - FFAC * ALF*GV2*AMYNI2/16D0
        XVNLC= XVNLC+ FFAC * ALF*GV2*(AYPNI2/4D0-AMYNI)*QI2F2/2D0

        XVQ2= XVQ * QIQF

        DO 6 K=1,JMAX
          U(K)=ELN(K+1)
          R(K)=ERN(K+1)
          S(K)=(K*ERN(K)+(K+1)*ERN(K+2))/(2*K+1)
          G(K)=(ELN(K+2)-ELN(K))/(2*K+1)
    6   CONTINUE
        U(0)=ELN(1)
        R(0)=ERN(1)
        S(0)=ERN(2)
        ARG=(QIPF2+AMES2)/ALAM2
        G(0)=(ELN(2)-ELN(1))-FDEXP(-QIPF2/ALAM2)*CE1(ARG)
        VC(-1) = VC(-1) + (XVC+X*YVC+X*X*ZVC)*U(JM)
     .                  - (YVC+X*ZVC)*R(JM) - ZVC*S(JM)
        VC( 0) = VC( 0) + (XVC+X*YVC+X*X*ZVC)*U(J)
     .                  - (YVC+X*ZVC)*R(J)  - ZVC*S(J)
        VC( 1) = VC( 1) + (XVC+X*YVC+X*X*ZVC)*U(JP)
     .                  - (YVC+X*ZVC)*R(JP) - ZVC*S(JP)
        VS(-1) = VS(-1) + (XVS+X*YVS+X*X*ZVS)*U(JM)
     .                  - (YVS+X*ZVS)*R(JM) - ZVS*S(JM)
        VS( 0) = VS( 0) + (XVS+X*YVS+X*X*ZVS)*U(J)
     .                  - (YVS+X*ZVS)*R(J)  - ZVS*S(J)
        VS( 1) = VS( 1) + (XVS+X*YVS+X*X*ZVS)*U(JP)
     .                  - (YVS+X*ZVS)*R(JP) - ZVS*S(JP)
        VT(-2) = VT(-2) + (XVT+X*YVT)*U(JMM)- YVT*R(JMM)
        VT(-1) = VT(-1) + (XVT+X*YVT)*U(JM) - YVT*R(JM)
        VT( 0) = VT( 0) + (XVT+X*YVT)*U(J)  - YVT*R(J)
        VT( 1) = VT( 1) + (XVT+X*YVT)*U(JP) - YVT*R(JP)
        VT( 2) = VT( 2) + (XVT+X*YVT)*U(JPP)- YVT*R(JPP)
        VLS(-2)= VLS(-2)+ (XVLS+X*YVLS)*U(JMM)- YVLS*R(JMM)
        VLS(-1)= VLS(-1)+ (XVLS+X*YVLS)*U(JM) - YVLS*R(JM)
        VLS( 0)= VLS( 0)+ (XVLS+X*YVLS)*U(J)  - YVLS*R(J)
        VLS( 1)= VLS( 1)+ (XVLS+X*YVLS)*U(JP) - YVLS*R(JP)
        VLS( 2)= VLS( 2)+ (XVLS+X*YVLS)*U(JPP)- YVLS*R(JPP)
        VQ(-3) = VQ(-3) +  XVQ*U(JMMM)
        VQ(-2) = VQ(-2) +  XVQ*U(JMM)
        VQ(-1) = VQ(-1) +  XVQ*U(JM)
        VQ( 0) = VQ( 0) +  XVQ*U(J)
        VQ( 1) = VQ( 1) +  XVQ*U(JP)
        VQ( 2) = VQ( 2) +  XVQ*U(JPP)
        VQ( 3) = VQ( 3) +  XVQ*U(JPPP)
        IF(NLOC.NE.0) THEN
          VC(-1) = VC(-1) + XVNLC*U(JM)
          VC( 0) = VC( 0) + XVNLC*U(J)
          VC( 1) = VC( 1) + XVNLC*U(JP)
        ENDIF
        IF(NLOC.EQ.2) THEN
          VS(-1) = VS(-1) + XVNLS*U(JM)
          VS( 0) = VS( 0) + XVNLS*U(J)
          VS( 1) = VS( 1) + XVNLS*U(JP)
        ENDIF

        IF(NQ12.EQ.1) THEN
C**       Extra contribution from inverse Fourier transform of Q12
          VQP(-2) = VQP(-2) + XVQ2*G(JMM)
          VQP(-1) = VQP(-1) + XVQ2*G(JM)
          VQP( 0) = VQP( 0) + XVQ2*G(J)
          VQP( 1) = VQP( 1) + XVQ2*G(JP)
          VQP( 2) = VQP( 2) + XVQ2*G(JPP)
        ENDIF

 2000 CONTINUE

      RETURN
      END
************************************************************************
      SUBROUTINE SCALAR(TYPE,NLOC,NQ12)
C*    Partial-wave momentum-space potentials for scalar exchange
C*       NQ12=0: inexact Fourier of quadratic spin-orbit Q12 operator
C*            1: exact Fourier Q12 from coordinate to momentum
C*       NLOC=0: no nonlocal (q^2+k^2/4) contributions
C*            1: nonlocal contributions in central potential
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER TYPE*2
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      COMMON/JFAC/  JMMM,JMM,JM,J,JP,JPP,JPPP
      COMMON/QFAC/  QI2,QF2,QIQF,QI2F2,QI2F22,QIPF2
      COMMON/MESONM/AMPI,AMETA,AMETP, AMRO,AMOM,AMFI, AMA0,AMEP,AMF0,
     .              AMPIC,AMROC,AMA0C,AWPIC,AWROC,AWA0C, AVSC
      COMMON/COPLNG/ALPV,THPV,PV1, FPI,FETA,FETP,FPI2,FETA2,FETP2,
     .              ALVE,THV ,GV1, GRO,GOM,GFI,  GRO2,GOM2,GFI2,
     .              ALVM,     FV1, FRO,FOM,FFI,  FRO2,FOM2,FFI2,
     .              ALGS,THGS,GS1, GA0,GEP,GF0,  GA02,GEP2,GF02,
     .              GFPRO,GFPOM,GFPFI, GFMRO,GFMOM,GFMFI,
     .              FPIC2,GROC2,FROC2,GFPROC,GFMROC,GA0C2
      COMMON/YUKEFF/ARO,AMR1,AROC,AMRC1,AWRC1,BRO,AMR2,BROC,AMRC2,AWRC2,
     .              AVSC1,AVSC2, AEPS,AME1,BEPS,AME2
      COMMON/CUTOFF/ALAM,ALAMP,ALAMV,ALAMS
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2
      DIMENSION U(-3:12), R(-3:12), G(-3:12), ELN(15), ERN(15)
      DATA U/16*0D0/, R/16*0D0/, G/16*0D0/
      DATA PI/3.14159265358979323846D0/

      XCOM=0.5D0*QI2F2/QIQF
      ALAM2=ALAMS*ALAMS
      Y=2D0*QIQF/ALAM2
      JMAX=J+3
      KOM=1
      FAC=2D0*PI/QIQF

C**     Scalar mesons: a0, broad epsilon (2 Yukawa's), f0
      DO 1000 IN=1,4
        IF(IN.EQ.1) THEN
          AMES=AMA0
          GS2 =GA02
        ELSEIF(IN.EQ.2) THEN
          AMES=AME1
          GS2 =GEP2*AEPS
        ELSEIF(IN.EQ.3) THEN
          AMES=AME2
          GS2 =GEP2*BEPS
        ELSEIF(IN.EQ.4) THEN
          AMES=AMF0
          GS2 =GF02
        ENDIF
        AMES2=AMES*AMES
        X=XCOM+0.5D0*AMES2/QIQF
        CALL SEX(X,Y,AMES/ALAMS,JMAX+1,KOM,ELN,ERN)
        XSC =-FAC * GS2*(1D0+AYPNI2/8D0*QI2F2)
        YSC = FAC * GS2*AYPNI2/4D0*QIQF
        XSLS=-FAC * GS2*AYPNI2/4D0
        XSQ = FAC * GS2*AMYNI2/16D0
        XSNL= FAC * GS2*AMYNI/4D0*QI2F2
        XSQ2= XSQ * QIQF
        DO 5 K=1,JMAX
          U(K)=ELN(K+1)
          R(K)=ERN(K+1)
          G(K)=(ELN(K+2)-ELN(K))/(2*K+1)
    5   CONTINUE
        U(0)=ELN(1)
        R(0)=ERN(1)
        ARG=(QIPF2+AMES2)/ALAM2
        G(0)=(ELN(2)-ELN(1))-FDEXP(-QIPF2/ALAM2)*CE1(ARG)
        VC(-1) = VC(-1) + (XSC+X*YSC)*U(JM) - YSC*R(JM)
        VC( 0) = VC( 0) + (XSC+X*YSC)*U(J)  - YSC*R(J)
        VC( 1) = VC( 1) + (XSC+X*YSC)*U(JP) - YSC*R(JP)
        VLS(-2)= VLS(-2)+ XSLS*U(JMM)
        VLS(-1)= VLS(-1)+ XSLS*U(JM)
        VLS( 0)= VLS( 0)+ XSLS*U(J)
        VLS( 1)= VLS( 1)+ XSLS*U(JP)
        VLS( 2)= VLS( 2)+ XSLS*U(JPP)
        VQ(-3) = VQ(-3) + XSQ*U(JMMM)
        VQ(-2) = VQ(-2) + XSQ*U(JMM)
        VQ(-1) = VQ(-1) + XSQ*U(JM)
        VQ( 0) = VQ( 0) + XSQ*U(J)
        VQ( 1) = VQ( 1) + XSQ*U(JP)
        VQ( 2) = VQ( 2) + XSQ*U(JPP)
        VQ( 3) = VQ( 3) + XSQ*U(JPPP)
        IF(NLOC.NE.0) THEN
          VC(-1) = VC(-1) + XSNL*U(JM)
          VC( 0) = VC( 0) + XSNL*U(J)
          VC( 1) = VC( 1) + XSNL*U(JP)
        ENDIF

        IF(NQ12.EQ.1) THEN
C**       Extra contribution from inverse Fourier transform of Q12
          VQP(-2) = VQP(-2) + XSQ2*G(JMM)
          VQP(-1) = VQP(-1) + XSQ2*G(JM)
          VQP( 0) = VQP( 0) + XSQ2*G(J)
          VQP( 1) = VQP( 1) + XSQ2*G(JP)
          VQP( 2) = VQP( 2) + XSQ2*G(JPP)
        ENDIF

 1000 CONTINUE

      IF(TYPE.EQ.'PP' .OR. TYPE.EQ.'NN') RETURN

C**     Charged scalar: a0+
      AMES=AWA0C
      GS2 =GA0C2
      AMES2=AMES*AMES
      X=XCOM+0.5D0*AMES2/QIQF
      CALL SEX(X,Y,AMES/ALAMS,JMAX+1,KOM,ELN,ERN)
      XSC =-FAC * ALF*GS2*(1D0+AMYNI/4D0*QI2F2)
      YSC = FAC * ALF*GS2*AMYNI/2D0*QIQF
      XSLS=-FAC * ALF*GS2*AMYNI/2D0
      XSQ = FAC * ALF*GS2*AMYNI2/16D0
      XSNL=-FAC * ALF*GS2*(AYPNI2/4D0-AMYNI)*QI2F2/2D0
      XSQ2= XSQ * QIQF
      DO 6 K=1,JMAX
        U(K)=ELN(K+1)
        R(K)=ERN(K+1)
        G(K)=(ELN(K+2)-ELN(K))/(2*K+1)
    6 CONTINUE
      U(0)=ELN(1)
      R(0)=ERN(1)
      ARG=(QIPF2+AMES2)/ALAM2
      G(0)=(ELN(2)-ELN(1))-FDEXP(-QIPF2/ALAM2)*CE1(ARG)
      VC(-1) = VC(-1) + (XSC+X*YSC)*U(JM) - YSC*R(JM)
      VC( 0) = VC( 0) + (XSC+X*YSC)*U(J)  - YSC*R(J)
      VC( 1) = VC( 1) + (XSC+X*YSC)*U(JP) - YSC*R(JP)
      VLS(-2)= VLS(-2)+ XSLS*U(JMM)
      VLS(-1)= VLS(-1)+ XSLS*U(JM)
      VLS( 0)= VLS( 0)+ XSLS*U(J)
      VLS( 1)= VLS( 1)+ XSLS*U(JP)
      VLS( 2)= VLS( 2)+ XSLS*U(JPP)
      VQ(-3) = VQ(-3) + XSQ*U(JMMM)
      VQ(-2) = VQ(-2) + XSQ*U(JMM)
      VQ(-1) = VQ(-1) + XSQ*U(JM)
      VQ( 0) = VQ( 0) + XSQ*U(J)
      VQ( 1) = VQ( 1) + XSQ*U(JP)
      VQ( 2) = VQ( 2) + XSQ*U(JPP)
      VQ( 3) = VQ( 3) + XSQ*U(JPPP)
      IF(NLOC.NE.0) THEN
        VC(-1) = VC(-1) + XSNL*U(JM)
        VC( 0) = VC( 0) + XSNL*U(J)
        VC( 1) = VC( 1) + XSNL*U(JP)
      ENDIF

      IF(NQ12.EQ.1) THEN
C**     Extra contribution from inverse Fourier transform of Q12
        VQP(-2) = VQP(-2) + XSQ2*G(JMM)
        VQP(-1) = VQP(-1) + XSQ2*G(JM)
        VQP( 0) = VQP( 0) + XSQ2*G(J)
        VQP( 1) = VQP( 1) + XSQ2*G(JP)
        VQP( 2) = VQP( 2) + XSQ2*G(JPP)
      ENDIF

      RETURN
      END
************************************************************************
      SUBROUTINE DIFRAC(TYPE,NLOC,NQ12,ISIGN,FACISO)
C*    Partial-wave momentum-space potentials for diffractive part
C*       NQ12=0: inexact Fourier of quadratic spin-orbit Q12 operator
C*            1: exact Fourier Q12 from coordinate to momentum
C*       NLOC=0: no nonlocal (q^2+k^2/4) contributions
C*            1: nonlocal contributions in central potential
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER TYPE*2
      COMMON/POTMOM/VC(-1:1),VS(-1:1),VT(-2:2),VLS(-2:2),VLA(-2:2),
     .              VQ(-3:3),VQP(-2:2)
      COMMON/JFAC/  JMMM,JMM,JM,J,JP,JPP,JPPP
      COMMON/QFAC/  QI2,QF2,QIQF,QI2F2,QI2F22,QIPF2
      COMMON/POMRON/GPOM,FPOM2,AMPOM,AMPOM2,AMPOM4,
     .              GA2D,FA2D2,AMA2D,AMA2D2,AMA2D4
      COMMON/AMCOEF/ALF,REDM, AMY,AMY2,AMYI,AMYI2, AMN,AMN2,AMNI,AMNI2,
     .              AMYPN,AMYMN, AMYPN2,AMYMN2, AY2PN2,AY2MN2,
     .              AMYN,AMYNI,AMYNI2, AYPNI2,AYMNI2
      DIMENSION R(-3:12), S(-3:12), G(-3:12), ELN(15), ERN(15)
      DATA R/16*0D0/, S/16*0D0/, G/16*0D0/
      DATA PI/3.14159265358979323846D0/

      XCOM=0.5D0*QI2F2/QIQF
      JMAX=J+3
      KOM=2
      FAC=4D0*PI

C**     Diffractive contribution: pomeron (f, f', A2), pomeron'
      DO 1000 IN=1,2
        IF(IN.EQ.1) THEN
          AMES2=AMPOM2
          GD2=FPOM2
        ELSEIF(IN.EQ.2) THEN
          AMES2=AMA2D2
          GD2=FA2D2*ISIGN
        ENDIF
        Y=0.5D0*QIQF/AMES2
        X=XCOM
        CALL SEX(X,Y,0D0,JMAX+1,KOM,ELN,ERN)
        XDC = FAC * GD2*(1D0+AYPNI2/8D0*QI2F2)
        YDC =-FAC * GD2*AYPNI2/4D0*QIQF
        XDLS= FAC * GD2*AYPNI2/4D0
        XDQ =-FAC * GD2*AMYNI2/16D0
        XDNL=-FAC * GD2*AMYNI/4D0*QI2F2
        XDQ2= XDQ * QIQF
        DO 5 K=1,JMAX
          R(K)=ERN(K+1)
          S(K)=(K*ERN(K)+(K+1)*ERN(K+2))/(2*K+1)
          G(K)=(ERN(K+2)-ERN(K))/(2*K+1)
    5   CONTINUE
        R(0)=ERN(1)
        S(0)=ERN(2)
        G(0)=(ERN(2)-ERN(1))-2D0*FDEXP(-QIPF2/4D0/AMES2)*AMES2/QIQF
        VC(-1) = VC(-1) + XDC*R(JM) + YDC*S(JM)
        VC( 0) = VC( 0) + XDC*R(J)  + YDC*S(J)
        VC( 1) = VC( 1) + XDC*R(JP) + YDC*S(JP)
        VLS(-2)= VLS(-2)+ XDLS*R(JMM)
        VLS(-1)= VLS(-1)+ XDLS*R(JM)
        VLS( 0)= VLS( 0)+ XDLS*R(J)
        VLS( 1)= VLS( 1)+ XDLS*R(JP)
        VLS( 2)= VLS( 2)+ XDLS*R(JPP)
        VQ(-3) = VQ(-3) + XDQ*R(JMMM)
        VQ(-2) = VQ(-2) + XDQ*R(JMM)
        VQ(-1) = VQ(-1) + XDQ*R(JM)
        VQ( 0) = VQ( 0) + XDQ*R(J)
        VQ( 1) = VQ( 1) + XDQ*R(JP)
        VQ( 2) = VQ( 2) + XDQ*R(JPP)
        VQ( 3) = VQ( 3) + XDQ*R(JPPP)
        IF(NLOC.NE.0) THEN
          VC(-1) = VC(-1) + XDNL*R(JM)
          VC( 0) = VC( 0) + XDNL*R(J)
          VC( 1) = VC( 1) + XDNL*R(JP)
        ENDIF

        IF(NQ12.EQ.1) THEN
C**       Extra contribution from inverse Fourier transform of Q12
          VQP(-2) = VQP(-2) + XDQ2*G(JMM)
          VQP(-1) = VQP(-1) + XDQ2*G(JM)
          VQP( 0) = VQP( 0) + XDQ2*G(J)
          VQP( 1) = VQP( 1) + XDQ2*G(JP)
          VQP( 2) = VQP( 2) + XDQ2*G(JPP)
        ENDIF

 1000 CONTINUE

      IF(TYPE.EQ.'PP' .OR. TYPE.EQ.'NN') RETURN

C**     Charged diffractive contribution: pomeron'
      AMES2=AMA2D2
      GD2=FA2D2*FACISO
      Y=0.5D0*QIQF/AMES2
      X=XCOM
      CALL SEX(X,Y,0D0,JMAX+1,KOM,ELN,ERN)
      XDC = FAC * ALF*GD2*(1D0+AMYNI/4D0*QI2F2)
      YDC =-FAC * ALF*GD2*AMYNI/2D0*QIQF
      XDLS= FAC * ALF*GD2*AMYNI/2D0
      XDQ =-FAC * ALF*GD2*AMYNI2/16D0
      XDNL= FAC * ALF*GD2*(AYPNI2/4D0-AMYNI)*QI2F2/2D0
      XDQ2= XDQ * QIQF
      DO 6 K=1,JMAX
        R(K)=ERN(K+1)
        S(K)=(K*ERN(K)+(K+1)*ERN(K+2))/(2*K+1)
        G(K)=(ERN(K+2)-ERN(K))/(2*K+1)
    6 CONTINUE
      R(0)=ERN(1)
      S(0)=ERN(2)
      G(0)=(ERN(2)-ERN(1))-2D0*FDEXP(-QIPF2/4D0/AMES2)*AMES2/QIQF
      VC(-1) = VC(-1) + XDC*R(JM) + YDC*S(JM)
      VC( 0) = VC( 0) + XDC*R(J)  + YDC*S(J)
      VC( 1) = VC( 1) + XDC*R(JP) + YDC*S(JP)
      VLS(-2)= VLS(-2)+ XDLS*R(JMM)
      VLS(-1)= VLS(-1)+ XDLS*R(JM)
      VLS( 0)= VLS( 0)+ XDLS*R(J)
      VLS( 1)= VLS( 1)+ XDLS*R(JP)
      VLS( 2)= VLS( 2)+ XDLS*R(JPP)
      VQ(-3) = VQ(-3) + XDQ*R(JMMM)
      VQ(-2) = VQ(-2) + XDQ*R(JMM)
      VQ(-1) = VQ(-1) + XDQ*R(JM)
      VQ( 0) = VQ( 0) + XDQ*R(J)
      VQ( 1) = VQ( 1) + XDQ*R(JP)
      VQ( 2) = VQ( 2) + XDQ*R(JPP)
      VQ( 3) = VQ( 3) + XDQ*R(JPPP)
      IF(NLOC.NE.0) THEN
        VC(-1) = VC(-1) + XDNL*R(JM)
        VC( 0) = VC( 0) + XDNL*R(J)
        VC( 1) = VC( 1) + XDNL*R(JP)
      ENDIF

      IF(NQ12.EQ.1) THEN
C**     Extra contribution from inverse Fourier transform of Q12
        VQP(-2) = VQP(-2) + XDQ2*G(JMM)
        VQP(-1) = VQP(-1) + XDQ2*G(JM)
        VQP( 0) = VQP( 0) + XDQ2*G(J)
        VQP( 1) = VQP( 1) + XDQ2*G(JP)
        VQP( 2) = VQP( 2) + XDQ2*G(JPP)
      ENDIF

      RETURN
      END
************************************************************************
      FUNCTION FDEXP(X)
      IMPLICIT REAL*8(A-Z)
      IF(X.LE.-100D0) THEN
        FDEXP=0D0
      ELSE
        FDEXP=DEXP(X)
      ENDIF
      RETURN
      END
************************************************************************
      SUBROUTINE SEX(X,Y,AMES,J,KOM,ELN,ERN)
C-----------------------------------------------------------------------
C     SEX calculates :   0.5*DEXP(AMES**2)*
C                    INTEGRAL(-1,+1) (PN(Z)/(X-Z)*F(X,Y)) DZ for LN
C                                    (PN(Z)*F(X,Y))       DZ for RN
C                    PN Legendre function
C                    F(X,Y)=DEXP(Y*(Z-X))
C     argument 1<|X|, via Pade approximant for X*LOG(X) [function GRENS]  
C
C     METHOD:
C            Expand PN(Z) ( N=0,12 max) in powers of Z and calculate :
C
C            INTEGRAL(-1,+1)(DEXP(Y(Z-X)*Z**N) DZ       for TN
C            INTEGRAL(-1,+1)(DEXP(Y(Z-X)*Z**N/(X-Z)) DZ for UN
C
C            INPUT:    X : (P**2+Q**2+M**2)/(2*P*Q)
C                      Y : 2*P*Q/(CUTOFF**2)
C                      J : (total angular momentum +3)+1 for V,S,D meson
C                    KOM : 1 calculate ELN and ERN
C                          2 calculate only ERN (for "DIFRAC" routine)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ELN(15),ERN(15),ALN(49),BLN(7)
      DIMENSION UN(14),TN(35)
      DATA XYHH,XH,XMIN,XMED,EPS,EPS2,YH/
     .     1D20,1D20,-160D0,5D0,1D-13,1D-6,1D20/
      DATA ALN/1D0, 1D0, 3D0,-1D0, 5D0,-3D0,
     .         35D0,-30D0,3D0, 63D0,-70D0,15D0,
     .         231D0,-315D0,105D0,-5D0, 429D0,-693D0,315D0,-35D0,
     .         6435D0,-12012D0,6930D0,-1260D0,35D0,
     .         12155D0,-25740D0,18018D0,-4620D0,315D0,
     .         46189D0,-109395D0,90090D0,-30030D0,3465D0,-63D0,
     .         88179D0,-230945D0,218790D0,-90090D0,15015D0,-693D0,
     .         676039D0,-1939938D0,2078505D0,-1021020D0,225225D0,
     .                  -18018D0,231D0/
      DATA BLN/1D0,2D0,8D0,16D0,128D0,256D0,1024D0/

C-NOTE *** dimension UN=J+1, TN=J+2+20
      IF(J.GT.13) WRITE(*,*)
     .   ' *** SEX: J exceeds allowable maximum of JMAX=9 (+3+1)'
      IF(J.GT.13) STOP

C     CALL ERRSET(208,256,-1,1)
      DO 500 I=1,15
        ELN(I)=0D0
        ERN(I)=0D0
  500 CONTINUE
C**   Calculate basic quantities
      YM=Y-Y*X
      YP=Y+Y*X
      YX=Y*X
      IF((YM-AMES*AMES).LT.XMIN) RETURN
      ENU=FDEXP(AMES*AMES)
      EP=0D0
      IF(YM.GT.XMIN) EP=FDEXP(YM)*ENU/2D0
      EM=0D0
      IF(-YP.GT.XMIN) EM=FDEXP(-YP)*ENU/2D0

C**   If same momenta but different mesonmass: calculate only UN again
      XYH=YX-AMES*AMES
      IF(DABS(XYH-XYHH).LT.1D-10*XYH .AND. Y.EQ.YH) GOTO 50
    9 XYHH=XYH
      YH=Y
      SH=(EP-EM)/Y
      CH=(EP+EM)/Y
      MAX=J+2
      IF(X.GT.XMED) MAX=MAX+20
      AMAX=MAX
      CHECK=0.434*(-AMAX+(Y+AMAX+0.5)*DLOG(1+(AMAX+1)/Y))
      IF(CHECK.GT.7D0) GOTO 10

C**   If CHECK > 1D7 use recurrence relation backward

C**   Calculate TN with recurrence relation forward
      TN(1)=SH
      DO 1 I=2,MAX
        TN(I)=CH-(I-1)*TN(I-1)/Y
        DUM=CH
        CH=SH
        SH=DUM
    1 CONTINUE
      GOTO 51

C**   Calculate TN with recurrence relation backward
C*    Calculation starting point
   10 IF(Y.LT.1D-04) MMAX=AMAX+4
      IF(Y.LT.1D-04) GOTO 13
      GRENS1=9.197/Y
      GRENS2=AMAX*DLOG((AMAX+1)/2.718/Y)/2.718/Y+8.829/Y
      GRENS1=GRENS(GRENS1)*2.718*Y
      GRENS2=GRENS(GRENS2)*2.718*Y
      K1=GRENS1+1
      K2=GRENS2+1
      MMAX=MAX0(K1,K2)
C**   CHECK: on loss of more than 7 significant digits
      CHECK=0.434*(Y*DLOG(Y)+Y+1-(Y+0.5)*DLOG(Y+1))
C-V   IF(CHECK.GT.7D0) PRINT 5003,CHECK

C*    Start backward recursion
   13 T=0D0
      IF(MOD(MMAX,2).EQ.0) MMAX=MMAX+1
      SH=ENU*DSINH(Y)*FDEXP(-YX)
      CH=EP+EM
      DO 11 II=2,MMAX
        I=MMAX-II+2
        T=(SH-Y*T)/(I-1)
        IF(I.LE.(MAX+1)) TN(I-1)=T
        DUM=SH
        SH=CH
        CH=DUM
   11 CONTINUE
      GOTO 51

C**   Calculation UN functions
   50 IF(X.GT.XMED .AND. XH.LE.XMED) GOTO 9

C**   If mesonmass differs, but X has passed critical value
C**   more TN's have to be calculated

   51 XH=X
C**   For diffractive contribution calculation of RN is enough
      IF(KOM.EQ.2) GOTO 800

      IF(X.GT.XMED) GOTO 60

C**   Calculation for X < XMED using forward recurrence relation
      CE1P=EP*CE1(-YM)
      CE1M=EM*CE1(YP)
C**   CHECK: on loss of more than 7 significant digits
      IF(CE1P.NE.0D0 .OR. CE1M.NE.0D0) THEN
        IF(DABS((CE1P-CE1M)/(CE1P+CE1M)).LT.EPS2) PRINT 5001
      ENDIF
      UN(1)=CE1P-CE1M
      MAX=J+1
      DO 52 I=2,MAX
        UN(I)=X*UN(I-1)-TN(I-1)
   52 CONTINUE
      GOTO 600

C**   Calculation for X > XMED using backward recurrence relation
   60 MAX=J+1

C**   First calculate UN(MAX) with expansion in TN
      SUM=TN(MAX)/X
      IF(SUM.EQ.0D0) GOTO 62
      XI=1D0/X
      DSUM=0D0
      DO 61 I=1,20
        SUM=SUM+DSUM
        XI=XI/X
        DSUM=TN(MAX+I)*XI
        IF(DABS(DSUM/SUM).LT.EPS) GOTO 62
   61 CONTINUE
      PRINT 5002
   62 UN(MAX)=SUM
      DO 63 II=2,MAX
        I=MAX+2-II
        UN(I-1)=(UN(I)+TN(I-1))/X
   63 CONTINUE
      GOTO 600

C**   Form linear combinations to obtain LN and RN
  600 MAX=J+1
      NB=1
      DO 31 L1=1,MAX
        II=L1+MOD(L1,2)
        DO 41 N1=1,L1,2
          ELN(L1)=ELN(L1)+ALN(NB+(N1-1)/2)*UN(L1+1-N1)/BLN(II/2)
          ERN(L1)=ERN(L1)+ALN(NB+(N1-1)/2)*TN(L1+1-N1)/BLN(II/2)
   41   CONTINUE
        NB=NB+II/2
   31 CONTINUE
C     CALL ERRSET(208,256,1,1)
      RETURN
  800 MAX=J+1
      NB=1
      DO 32 L1=1,MAX
        II=L1+MOD(L1,2)
        DO 42 N1=1,L1,2
          ERN(L1)=ERN(L1)+ALN(NB+(N1-1)/2)*TN(L1+1-N1)/BLN(II/2)
   42   CONTINUE
        NB=NB+II/2
   32 CONTINUE
C     CALL ERRSET(208,256,1,1)
      RETURN
 5001 FORMAT(1X,'**** SEX **** Loss of significant digits UN(1) ')
 5002 FORMAT(1X,'**** SEX **** UN MAX not accurate  ')
 5003 FORMAT(1X,'**** SEX **** TN not accurate ; loss of :',D10.3,
     .          ' decimal digits ')
      END
************************************************************************
      DOUBLE PRECISION FUNCTION GRENS(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA A0,A1,A2,A3/-0.589654,-0.0595734,0.649227,0.1809910/

      IF(X.GT.1.2D6) THEN
        GRENS=X**0.847
      ELSEIF(X.GT.9D4) THEN
        GRENS=X**0.825
      ELSEIF(X.GT.6900D0) THEN
        GRENS=X**0.806
      ELSEIF(X.GT.460D0) THEN
        GRENS=X**0.781
      ELSEIF(X.GT.23D0) THEN
        GRENS=X**0.751
      ELSE
        GRENS=(X*A3-A1+DSQRT((A1-X*A3)*(A1-X*A3)-4D0*A2*(A0-X)))/2D0/A2
      ENDIF

      RETURN
      END
************************************************************************
      DOUBLE PRECISION FUNCTION CE1(X)
C**     CE1 calculates the function : DEXP(X) * E1(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA EPS/10D-14/
      DATA GAM/0.577215664901532861D0/
      DATA AP0,AP1,AP2,AP3,AP4,AP5,AP6 /
     .     0.463996004278035D+01, 0.127788778637147D+03,
     .     0.735910238555843D+03, 0.139583023127254D+04,
     .     0.101614779141469D+04, 0.286647946600883D+03,
     .     0.256489038620717D+02/
      DATA AQ1,AQ2,AQ3,AQ4,AQ5,AQ6,AQ7 /
     .     0.512251050448444D+02, 0.503800829553457D+03,
     .     0.165169408854742D+04, 0.220150290642078D+04,
     .     0.127719811988873D+04, 0.312294439564262D+03,
     .     0.256489625816454D+02/
      DATA BP0,BP1,BP2,BP3,BP4,BP5     /
     .     0.335956527252693D+01, 0.204955591333077D+02,
     .     0.267757325223533D+02, 0.112883678215773D+02,
     .     0.164680678114210D+01, 0.655193572680895D-01/
      DATA BQ1,BQ2,BQ3,BQ4,BQ5,BQ6     /
     .     0.143836492361913D+02, 0.400563387674630D+02,
     .     0.366148021121537D+02, 0.128696120312766D+02,
     .     0.171232738644327D+01, 0.655193403549186D-01/
      DATA CP0,CP1,CP2,CP3,CP4         /
     .     0.298014627030798D+01, 0.113803314436134D+02,
     .     0.947288802836929D+01, 0.247747160891423D+01,
     .     0.188516317695352D+00/
      DATA CQ1,CQ2,CQ3,CQ4,CQ5         /
     .     0.988019055335016D+01, 0.189408176576544D+02,
     .     0.117618585876339D+02, 0.266598761793551D+01,
     .     0.188516320637495D+00/
      DATA DP0,DP1,DP2,DP3   /
     .     0.242331331460798D+01, 0.432777141801875D+01,
     .     0.160959648287707D+01, 0.148720388893508D+00/
      DATA DQ1,DQ2,DQ3,DQ4   /
     .     0.558734308280980D+01, 0.578865453197840D+01,
     .     0.175831677540018D+01, 0.148720389489176D+00/
      DATA EP0,EP1,EP2,EP3   /
     .     0.226526458912153D+01, 0.332000741007556D+01,
     .     0.104761178441346D+01, 0.837423061701825D-01/
      DATA EQ1,EQ2,EQ3,EQ4   /
     .     0.478887726713541D+01, 0.428387700117901D+01,
     .     0.113135408983342D+01, 0.837423061723804D-01/
      DATA FP0,FP1,FP2,FP3   /
     .     0.190053654321203D+01, 0.151285969203750D+01,
     .     0.205314346964057D+00, 0.264152351883344D-03/
      DATA FQ1,FQ2,FQ3,FQ4   /
     .     0.320887608816311D+01, 0.171790987670629D+01,
     .     0.205578499347658D+00, 0.264152351839874D-03/
 
      PQX=1D0
      IF(X.LT.1D0) THEN
        E1=-GAM-DLOG(X)+X
        N=1
        IX=1
        TERM=X
    1   N=N+1
        IX=-IX
        TERM=TERM*X/N
        IF(TERM.LT.EPS) GOTO 10
        E1=E1+IX*TERM/N
        GOTO 1
   10   PPX=FDEXP(X)*E1
      ELSEIF(X.LT.3D0) THEN
        PPX=AP0+X*(AP1+X*(AP2+X*(AP3+X*(AP4+X*(AP5+X*AP6)))))
        PQX=1D0+X*(AQ1+X*(AQ2+X*(AQ3+X*(AQ4+X*(AQ5+X*(AQ6+X*AQ7))))))
      ELSEIF(X.LT.6D0) THEN
        PPX=BP0+X*(BP1+X*(BP2+X*(BP3+X*(BP4+X*BP5))))
        PQX=1D0+X*(BQ1+X*(BQ2+X*(BQ3+X*(BQ4+X*(BQ5+X*BQ6)))))
      ELSEIF(X.LT.14D0) THEN
        PPX=CP0+X*(CP1+X*(CP2+X*(CP3+X*CP4)))
        PQX=1D0+X*(CQ1+X*(CQ2+X*(CQ3+X*(CQ4+X*CQ5))))
      ELSEIF(X.LT.25D0) THEN
        PPX=DP0+X*(DP1+X*(DP2+X*DP3))
        PQX=1D0+X*(DQ1+X*(DQ2+X*(DQ3+X*DQ4)))
      ELSEIF(X.LT.70D0) THEN
        PPX=EP0+X*(EP1+X*(EP2+X*EP3))
        PQX=1D0+X*(EQ1+X*(EQ2+X*(EQ3+X*EQ4)))
      ELSEIF(X.LT.165D0) THEN
        PPX=FP0+X*(FP1+X*(FP2+X*FP3))
        PQX=1D0+X*(FQ1+X*(FQ2+X*(FQ3+X*FQ4)))
      ELSEIF(X.GE.165D0) THEN
        Y=1D0/X
        N=0
        IN=1
        K=1
        CE1=Y
    2   N=N+1
        IN=-IN
        Y=Y/X
        K=K*N
        TERM=Y*K
        IF(TERM/CE1.LT.EPS) GOTO 20
        CE1=CE1+IN*TERM
        GOTO 2
   20   PPX=CE1
      ENDIF
      CE1=PPX/PQX
 
      RETURN
      END
c***************** this program consistently in double precision ****
c name:    dminv
c        programmbibliothek rhrz bonn        28/11/78       dminv
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c invert a matrix
c
c usage:   call dminv (a,n,d,l,m)
c
c parameters:
c
c a:       input matrix, destroyed in computation and replaced by
c          resultant inverse.
c          double precision required.
c
c n:       order of matrix a
c
c d:       resultant determinant
c          double precision required.
c
c l:       work vector of length n
c
c m:       work vector of length n
c
c remarks: matrix a must be a general matrix
c
c method:
c
c the standard gauss-jordan method is used. the determinant
c is also calculated. a determinant of zero indicates that
c the matrix is singular.
c
c programs required:
c          none
c
c author:  ibm, ssp iii
c
c**********************************************************************
      subroutine dminv (a,n,d,l,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),l(1),m(1)
c
c
c        search for largest element
c
      d=1.d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (dabs(biga)-dabs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 if(biga) 48,46,48
   46 d=0.d0
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end

