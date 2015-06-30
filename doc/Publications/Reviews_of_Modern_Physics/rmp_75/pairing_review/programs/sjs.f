

      PROGRAM  ph_screen
      IMPLICIT REAL*8 (A-H,O-Z)

      read(5,*) j1, j2
      call binom
      jc = 1
      ja = 1
      jb = 1
      jd = 1
      write (6,*) j1, j2, sjs(jc,jd,j1,jb,ja,j2) 
      END



C
C             *********************************************
C             *    FUNCTION TO CHECK # OSC. EXCITATIONS   *
C             *********************************************
C

      LOGICAL FUNCTION dencheck(i)
      IMPLICIT REAL*8 (A-H,O-Z)    
      COMMON/energy/a(5),wstart,wcn(11),e(1000),ihomega
      dencheck = ((IABS(i).EQ.0).OR.(IABS(i).GT.ihomega))

      END


C
C             *********************************************
C             *    FUNCTION TO CHECK SP LIMITS            *
C             *********************************************
C

      LOGICAL FUNCTION ptest(jp,lim)
      IMPLICIT REAL*8 (A-H,O-Z)    
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      ptest=.TRUE.
      IF(lim.LE.nn1) THEN
        IF((MOD(jp,2).EQ.0).AND.(jp.LE.nn1)) ptest=.FALSE.
        IF((MOD(jp,2).NE.0).AND.(jp.LE.np1)) ptest=.FALSE.
      ELSEIF(lim.LE.nn2) THEN
        IF((MOD(jp,2).EQ.0).AND.(jp.GT.nn1)) ptest=.False.
        IF((MOD(jp,2).NE.0).AND.((jp.GT.np1).AND.
     :      (jp.LE.np2))) ptest=.FALSE.
      ENDIF

      END

C
C             *********************************************
C             *    FUNCTION TO CHECK TRIANGULAR RELS.     *
C             *********************************************
C


      LOGICAL FUNCTION triag(i,j,k)
      IMPLICIT REAL*8 (A-H,O-Z)    
      triag = ((i-j-k)*(i-IABS(j-k)).GT.0)

      END

C
C             *********************************************
C             *    FUNCTION TO CALCULATE PHASE FACTORS    *
C             *********************************************
C


      INTEGER FUNCTION iph(n)
      IMPLICIT REAL*8 (A-H,O-Z)    
      iph=1-2*MOD(n,2)

      END



C
C             *********************************************
C             *    FUNCTION TO CALCULATE NORM OF G-MAT    *
C             *********************************************
C


      REAL FUNCTION dij*8(ja,jb)
      IMPLICIT REAL*8 (A-H,O-Z)    
      IF(ja.EQ.jb) THEN
        dij=1.4142135623731d0
      ELSE
        dij=1.D0
      ENDIF

      END




C          ********************************************
c          *         calculates ninej-symbols         *
c          ********************************************

      REAL FUNCTION snj*8(ia,ib,ie,ic,id,IF,ig,ih,it)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/binoco/ q(50,50),kh(200)
      ja=ia+1
      jb=ib+1
      jc=ic+1
      jd=id+1
      je=ie+1
      jf=IF+1
      jg=ig+1
      jh=ih+1
      jt=it+1
      i=kh(ja+jb-je+99)+kh(jb+je-ja+99)+kh(je+ja-jb+99)+kh(jc+jd-jf+99)
     : +kh(jd+jf-jc+99)+kh(jf+jc-jd+99)+kh(jg+jh-jt+99)+kh(jh+jt-jg+99)
     : +kh(jt+jg-jh+99)+kh(ja+jc-jg+99)+kh(jc+jg-ja+99)+kh(jg+ja-jc+99)
     : +kh(jb+jd-jh+99)+kh(jd+jh-jb+99)+kh(jh+jb-jd+99)+kh(je+jf-jt+99)
     : +kh(jf+jt-je+99)+kh(jt+je-jf+99)
      IF(i.GT.0) THEN
	snj=0.D0
      ELSE
	la=(ie+IF+it)/2+2
	ld=(ig+ih+it)/2+2
	ma=(ia+ic+ig)/2+2
	mc=(IF+ic+id)/2+2
	na=(ib+id+ih)/2+2
	nb=(ib+ie+ia)/2+2
	le=(ie+IF-it)/2+1
	lf=(IF+it-ie)/2+1
	lg=(it+ie-IF)/2+1
	me=(ia+ic-ig)/2+1
	mf=(ic+ig-ia)/2+1
	mg=(ig+ia-ic)/2+1
	ne=(ib+id-ih)/2+1
	nf=(id+ih-ib)/2+1
	ng=(ih+ib-id)/2+1
	lx=(it+ig-ih)/2+1
	mx=(ic+id-IF)/2+1
	nx=(ib+ie-ia)/2+1
	fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
	fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
	jsi=MAX0(IABS(je-jh),IABS(jg-jf),IABS(ja-jd))+1
	jsf=MIN0(je+jh,jg+jf,ja+jd)-1
	ps=-1-2*(jsi/2*2-jsi)
	fs=ps*DSQRT(fn/fd)/DFLOAT((ig+1)*(ie+1))
	u=0.D0
	DO js=jsi,jsf,2
	   is=js-1
	   lb=(ie+ih+is)/2+2
	   lc=(ig+IF+is)/2+2
	   mb=(ia+id+is)/2+2
	   ly=(ie+ih-is)/2+1
	   my=(ig+IF-is)/2+1
	   ny=(ia-id+is)/2+1
	   ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*
     :        q(js,ny)
	   l0=MAX0(la,lb,lc,ld)+1
	   m0=MAX0(ma,mb,mc,lc)+1
	   n0=MAX0(na,nb,mb,lb)+1
	   l1=MIN0(le+ld,lf+lb,lg+lc)
	   m1=MIN0(me+lc,mf+mb,mg+mc)
	   n1=MIN0(ne+lb,nf+nb,ng+mb)
	   x=0.D0
	   DO l=l0,l1
	      x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
	   ENDDO
	   y=0.D0
	   DO m=m0,m1
	      y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
	   ENDDO
	   z=0.D0
	   DO n=n0,n1
	      z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
	   ENDDO
	   ihx=l1+m1+n1
	   p=1+2*(ihx/2*2-ihx)
	   u=u+p*x*y*z/ud
	ENDDO
	snj=u*fs
      ENDIF
      END







c          ********************************************
c          *         calculates 6j-symbols            *
c          ********************************************


      REAL FUNCTION sjs*8(j_a,j_b,j_c,l_a,l_b,l_c)
      IMPLICIT REAL*8 (a-h,o-z)    
      COMMON/binoco/ q(50,50),kh(200)
      sjs=0.0d0
      ja=j_a + 1
      jb=j_b + 1
      jc=j_c + 1
      la=l_a + 1
      lb=l_b + 1
      lc=l_c + 1
      i=kh(ja+jb-jc+99)+kh(jb+jc-ja+99)+kh(jc+ja-jb+99)+kh(ja+lb-lc+99)
     : +kh(lb+lc-ja+99)+kh(lc+ja-lb+99)+kh(la+jb-lc+99)+kh(jb+lc-la+99)
     : +kh(lc+la-jb+99)+kh(la+lb-jc+99)+kh(lb+jc-la+99)+kh(jc+la-lb+99)
      IF(i.LE.0) THEN
	mt=(j_a+j_b+j_c)/2 + 2
	ma=(j_a+l_b+l_c)/2+ 2
	mb=(l_a+j_b+l_c)/2+ 2
	mc=(l_a+l_b+j_c)/2+ 2
	na=mt-ja
	nb=mt-jb
	nc=mt-jc
	ka=ma-lc
	kb=mb-lc
	kc=mc-jc
	f=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)*q(la,kb)
     :    *q(mc,la+1)*q(la,kc))
	fs=DSQRT(f)/(1.d0*l_a + 1.0d0)
	l0=MAX0(mt,ma,mb,mc)+1
	l1=MIN0(ma+na,mb+nb,mc+nc)
	x=0.D0
	DO l=l0,l1
	   x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
	ENDDO
	sjs=-(1+2*(l1/2*2-l1))*fs*x
      ENDIF
      END


c          ********************************************
c          *         calculates 3j-symbols            *
c          ********************************************



      REAL FUNCTION tjs*8(j_a,j_b,j_c,m_a,m_b,m_c)
      IMPLICIT REAL*8 (a-h,o-z)
      COMMON/binoco/ q(50,50),kh(200)
      eps=1.0d-2
      tjs=0.D0
      ja=(j_a+m_a)/2+1
      ma=(j_a-m_a)/2+1
      jb=(j_b+m_b)/2+1
      mb=(j_b-m_b)/2+1
      jc=(j_c+m_c)/2+1
      mc=(j_c-m_c)/2+1
      la=(j_b+j_c-j_a)/2+1
      lb=(j_c+j_a-j_b)/2+1
      lc=(j_a+j_b-j_c)/2+1
      lt=(j_a+j_b+j_c)/2+1
      ld=MIN0(ja,jb,jc,ma,mb,mc,la,lb,lc)
      IF(((m_a+m_b+m_c).LE.0).AND.(ld.GT.0)) THEN
	ja2=j_a+m_a
	jb2=j_b+m_b
	jc2=j_c+m_c
	i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
	IF(i.EQ.0) then 
	  fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2)
     :      *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
	  k0=MAX0(0,lc-ja,lc-mb)+1
	  k1=MIN0(lc,ma,jb)
	  x=0.D0
	  DO k=k0,k1
	     x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
	  ENDDO
	  ip=k1+lb+jc
	  p=1-2*(ip-ip/2*2)
	  tjs=p*x*DSQRT(fn)
	ENDIF
      ENDIF
      END


c
c
c     ******************************************************************
c     *                factorials for 3j,6j and 9j symbols             *
c     ******************************************************************
c
c


      SUBROUTINE binom
      IMPLICIT REAL*8 (a-h,o-z)    
      COMMON/binoco/ q(50,50),kh(200)
      DO l=1,200
	 kh(l)=1
      ENDDO
      kh(100) =0
      DO l=1,50
	 q(l,1)=1.0d0
	 q(l,l)=1.0d0
	 kh(l+l+100)=0
      ENDDO
      DO l=2,49
	 DO k=2,l
	    q(l+1,k)=q(l,k-1)+q(l,k)
	 ENDDO
      ENDDO

      END



      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=250,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(NMAX),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
 
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit real*8(a-h,o-z)
      parameter(nmax=250)
      DIMENSION A(NP,NP),INDX(NMAX),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
 

