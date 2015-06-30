
c             ph-screening program for all vertices
c             please notify the author, morten hjorth-jensen at
c             mhjensen@ect.unitn.it, if changes are made
c             
c             ****  last update   : may 1996          ****             
c             ****  by            : m. hjorth-jensen  ****
c
c
c

      PROGRAM  ph_screen
      IMPLICIT REAL*8 (A-H,O-Z)

c      CALL read_data(iqbox,nbody,istart,iorder,i_effoper) ! read veff data
      CALL pauli                            ! find all model-space configs
      CALL binom                            ! factorials for 6j and 9j
      CALL qsbox

      END



c
c
c     ******************************************************************
c      *        This subrout reads the data which define the veffint   *
c     ******************************************************************
c
c



      SUBROUTINE read_data(iqbox,nbody,istart,iorder,i_effoper)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/jtzdat/itmin,itmax,jtmin,jtmax
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      COMMON /nlja/na(200),la(200),ja(200),nhbar(153)
      COMMON/print_opt/print
      LOGICAL print
      CHARACTER*70 gmatrix, qboxfile, veffint


 1000 format(2x,a70)
 1003 format(2x,5htmin=,i3,2x,5htmax=,i3,2x,5hjmin=,i3,2x,5hjmax=,i3)
 1001 format(f8.4,3i3,l4)
 1002 format(5f12.4)
 1004 format(5i4,i5,f10.4)

c        *****************************************************
c        *          read in the various variables            *
c        *****************************************************

      READ(5,*)
      READ(5,1000) gmatrix
      READ(5,1000) qboxfile
      READ(5,1000) veffint
      OPEN(UNIT=7,FILE=gmatrix)
      open(unit=8,file=veffint)       ! WRITE veffint on unit 8
      READ(5,*)
      READ(5,*) nnh, lnh, jnh      ! READ qm numbers of last neutron hole
      READ(5,*)
      READ(5,*) nnp, lnp, jnp      ! READ qm numbers of last neutron part
      READ(5,*)
      READ(5,*) nph, lph, jph      ! READ qm numbers of last proton hole
      READ(5,*)
      READ(5,*) npp, lpp, jpp      ! READ qm numbers of last proton part
      READ(5,*)
      READ(5,*) mnnh, mlnh, mjnh   ! READ qm numbers of 1st neutron mspace
      READ(5,*)
      READ(5,*) mnnp, mlnp, mjnp   ! READ qm numbers of last neutron mspace
      READ(5,*)
      READ(5,*) mnph, mlph, mjph   ! READ qm numbers of first proton mspace
      READ(5,*)
      READ(5,*) mnpp, mlpp, mjpp   ! READ qm numbers of last proton mspace



      CALL pspace(nnh,lnh,jnh,nnp,lnp,jnp,nph,lph,jph,npp,lpp,jpp,
     &     mnnh,mlnh,mjnh,mnnp,mlnp,mjnp,mnph,mlph,mjph,
     &     mnpp,mlpp,mjpp)


      IF(np1.GT.nn1) THEN 
        WRITE(6,*) ' warning & more proton orbs than
     &               neutron orbs in the closed shell core'
        STOP
      ENDIF
      READ(5,*)
      READ(5,*) itmin,itmax,jtmin,jtmax     ! min & max tz, min & max j   
      READ(5,*)
      READ(5,*) wstart                      ! starting energy
      READ(5,*)
      READ(5,*) ihomega                     ! #\hbar\omega restriction
      READ(5,*)
      READ(5,*) iorder                      ! order of the interaction
      READ(5,*)
      READ(5,*) nbody                       ! 1-body or 2-body int.
      READ(5,*)
      READ(5,*) i_effoper                   ! effective operator


c        *****************************************************
c        *     here we set up the energy variables to be     *
c        *     used in the lagrangian interpolation          *
c        *                                                   *
c        *     ----->      important  note        <-----     *
c        *                                                   *
c        *     observe that the definition of wcn            *
c        *     is given for ascending   starting energies    *
c        *     if you have a descending series of energies   *
c        *     invert the signs of (2.*i-1)-11)              *
c        *       descending wstart, e.g. -140 -90 -50 -20 -5 *
c        *     remember to make similar changes for the      *
c        *     derivatives as well                           *
c        *                                                   *
c        *     ascending wstart, e.g. -5 -10 -20 -40 -85     *
c        *****************************************************

      do i=1,11                            ! starting energy variable
	 wcn(i)=-1.d0*i+6.d0               ! used in the calc of denos
      enddo                                ! in various diagrams
      READ(5,*)
      READ(5,*)(a(i),i=1,5)                ! READ in starting energies
                                           ! for the g-matrix
      READ(5,*)
      READ(5,*) istart, iqbox
      READ(5,*)
      READ(5,*) mass_n                     ! READ #parts  of closed shell core
      READ(5,*)
      READ(5,*) print                      ! READ print option


c        *****************************************************
c        *     here we set up the proton and neutron sp      *
c        *     quantum numbers and their unperturbed energies*
c        *****************************************************

      hbar=45.d0*dfloat(mass_n)**(-1.d0/3.d0)-   ! sets up the osc.  energy
     &       25.d0*dfloat(mass_n)**(-2.d0/3.d0)  ! defined by total a=mass_n

      emod=wstart/2.d0                           ! the sp part ener of p-space
      shift=(dfloat(nhbar(mnf))+1.5d0)*hbar-emod



      do j = 1, nn2
         jj(j*2)=ja(j)                     ! neutrons are 
         ll(j*2)=la(j)                     ! always even numbers
         nn(j*2)=na(j)                     ! protons are odd
         nshell(j*2)=nhbar(j)
         e(j*2)=(dfloat(nhbar(j))+1.5d0)*hbar-shift
         jj(j*2-1)=ja(j)                    
         ll(j*2-1)=la(j)
         nn(j*2-1)=na(j)
         nshell(j*2-1)=nhbar(j)
         e(j*2-1)=(dfloat(nhbar(j))+1.5d0)*hbar-shift
         itzp(j*2-1)=-1                   ! protons, nucl. phys. def.
         itzp(j*2)=1                      ! neutrons, nucl. phys. def.
      enddo

      WRITE(6,*) '  the single particle data '
      WRITE(6,*) '  #&  n   l 2*j  tz hbar unpert-sp-energy '
      do j=1,nn2*2
         WRITE(6,1004) j,nn(j),ll(j),jj(j),itzp(j),nshell(j),
     &                 e(j)
      enddo
      WRITE(6,*) ' last proton hole:     ', 2*np1-1
      WRITE(6,*) ' last proton particle: ', 2*np2-1
      WRITE(6,*) ' last neutron hole:    ', 2*nn1
      WRITE(6,*) ' last neutron particle:', 2*nn2
      WRITE(6,*) ' first proton mspace:  ', 2*mpf-1
      WRITE(6,*) ' last proton mspace:   ', 2*mpl-1
      WRITE(6,*) ' first neutron mspace: ', 2*mnf
      WRITE(6,*) ' last neutron mspace : ', 2*mnl
      WRITE(6,1003)itmin,itmax,jtmin,jtmax
      WRITE(6,*) ' starting energy: ', wstart
      WRITE(6,*) ' # of harmonic oscillator quanta: ', ihomega 
      WRITE(6,*) ' order of the q-box: ', iorder 
      WRITE(6,*) ' type of interaction, n-body: ', nbody
      WRITE(6,*) ' if i_effoper=1, calc effoper, else only effint'
      WRITE(6,*) ' i_effoper:  ', i_effoper

      WRITE(6,*) ' starting energies for the g-mat: '
      WRITE(6,1002)(a(i),i=1,5)
      WRITE(6,*) ' if istart=0, only q-box, if istart=1, rs pert theory'
      WRITE(6,*) ' and if istart=2 folded diagrams'
      WRITE(6,*) ' istart: ', istart
      WRITE(6,*) ' if iqbox=1, calc q-box, else iqbox=2, READ from tape'
      WRITE(6,*) ' iqbox:  ', iqbox
      WRITE(6,*) ' mass number of closed shell core: ', mass_n
      WRITE(6,*) ' oscillator energy in mev:         ', hbar
      WRITE(6,*) ' the energy shift used:            ', shift



      nn1=2*nn1
      nn2=2*nn2
      np1=2*np1-1
      np2=2*np2-1
      mpf=2*mpf-1
      mpl=2*mpl-1
      mnf=2*mnf
      mnl=2*mnl

      CALL readme                           ! READ in the g-matrix

      END



     
c     *****************************************************************
c     *   sets up h.o. sp basis and the calculations limits in array  *
c     *         common/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl         *
c     *****************************************************************


      SUBROUTINE pspace(nnh,lnh,jnh,nnp,lnp,jnp,nph,lph,jph,npp,lpp,jpp,
     &     mnnh,mlnh,mjnh,mnnp,mlnp,mjnp,mnph,mlph,mjph,
     &     mnpp,mlpp,mjpp)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER ( lim=15)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON /nlja/na(200),la(200),ja(200),nhbar(153)
      DIMENSION ishell(lim)
      DATA nhbar/0,
     &        1,1,
     &        2,2,2,
     &        3,3,3,3,
     &        4,4,4,4,4,
     &        5,5,5,5,5,5,
     &        6,6,6,6,6,6,6,
     &        7,7,7,7,7,7,7,7,
     &        8,8,8,8,8,8,8,8,8,
     &        0,0,1,1,2,2,3,3,4,4,
     &        0,0,1,1,2,2,3,3,4,4,5,
     &        0,0,1,1,2,2,3,3,4,4,5,5,
     &        0,0,1,1,2,2,3,3,4,4,5,5,6,
     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,
     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,
     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,
     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8/



c      data nhbar/0,
c     &        1,1,
c     &        2,2,2,
c     &        3,4,4,4,4,
c     &        5,5,5,5,5,
c     &        6,6,6,6,6,
c     &        6,7,7,7,7,7,7,
c     &        7,8,8,8,8,8,8,8,
c     &        9,9,9,9,9,9,9,9,9,
c     &        0,0,1,1,2,2,3,3,4,4,
c     &        0,0,1,1,2,2,3,3,4,4,5,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,6,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,
c     &        0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8/


      DO i=1,lim
         ishell(i)=i-1                    ! gives 2n+l 
      ENDDO
      icount=0
      DO i=1,lim
         max_nl=ishell(i)                ! max value of 2n+l for this shell
         ipar=MOD(max_nl,2)              ! parity
         n_min=0                         ! min value of principal qm number
         n_max=max_nl/2                  ! max value of principal qm number
         l_min=max_nl-2*n_max            ! lowest orb. mom. for this shell 
         l_max=max_nl-2*n_min            ! highest orb. mom. for this shell
         j_min=2*l_min-1
         j_max=2*l_max+1
         IF(j_min.LE.0) j_min=1          ! make sure that ns_{1/2} is accounted
         DO jpart=j_max,j_min,-2
            icount=icount+1
            ja(icount)=jpart
            l1=(jpart+1)/2
            l2=(jpart-1)/2
            if(mod(l1,2).eq.ipar)la(icount)=l1
            if(mod(l2,2).eq.ipar)la(icount)=l2
            na(icount)=(max_nl-la(icount))/2
         ENDDO
      ENDDO

      DO j=1,icount
         ll=la(j)
         jj=ja(j)
         nn=na(j)
         if((nn.eq.nnh).and.(ll.eq.lnh).and.(jj.eq.jnh)) nn1=j
         if((nn.eq.nnp).and.(ll.eq.lnp).and.(jj.eq.jnp)) nn2=j
         if((nn.eq.nph).and.(ll.eq.lph).and.(jj.eq.jph)) np1=j
         if((nn.eq.npp).and.(ll.eq.lpp).and.(jj.eq.jpp)) np2=j
         if((nn.eq.mnnh).and.(ll.eq.mlnh).and.(jj.eq.mjnh))mnf=j
         if((nn.eq.mnnp).and.(ll.eq.mlnp).and.(jj.eq.mjnp))mnl=j
         if((nn.eq.mnph).and.(ll.eq.mlph).and.(jj.eq.mjph))mpf=j
         if((nn.eq.mnpp).and.(ll.eq.mlpp).and.(jj.eq.mjpp))mpl=j
      ENDDO
      END





     
c     *****************************************************************
c     *        calculates all possible model space   states           *
c     *****************************************************************

      SUBROUTINE pauli(jtot,itz,ipar,nconf,ipba) 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_conf=12000)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION ipba(max_conf)
      LOGICAL triag

      init_m=MIN0(mpf,mnf)
      last_m=MAX0(mpl,mnl)
      nconf=0
      DO ia=init_m,last_m
         la=ll(ia)
         ja=jj(ia)
         IF(itz.NE.0) THEN
            istart=ia
         ELSE
            istart=1
         ENDIF
         DO 10 ib=istart,last_m
            lb=ll(ib)
            jb=jj(ib)
            IF(MOD((la+lb),2).NE.ipar) GOTO 10
            IF(triag(2*jtot,ja,jb)) GOTO 10
            IF(itzp(ia)+itzp(ib).NE.2*itz) goto 10
            IF((ia.EQ.ib).AND.(MOD(jtot,2).NE.0)) GOTO 10
            IF((itzp(ia).EQ.-1).AND.(ia.GT.mpl)) GOTO 10
            IF((itzp(ia).EQ.-1).AND.(ia.LT.mpf)) GOTO 10
            IF((itzp(ia).EQ.1).AND.(ia.GT.mnl)) GOTO 10
            IF((itzp(ia).EQ.1).AND.(ia.LT.mnf)) GOTO 10
            IF((itzp(ib).EQ.-1).AND.(ib.GT.mpl)) GOTO 10
            IF((itzp(ib).EQ.-1).AND.(ib.LT.mpf)) GOTO 10
            IF((itzp(ib).EQ.1).AND.(ib.GT.mnl)) GOTO 10
            IF((itzp(ib).EQ.1).AND.(ib.LT.mnf)) GOTO 10
            IF(itzp(ia).GT.itzp(ib)) GOTO 10
            nconf=nconf+1
            k2=nconf*2
            k1=k2-1
            ipba(k1)=ia
            ipba(k2)=ib
 10      CONTINUE
      ENDDO
      END






      SUBROUTINE qsbox
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_conf=12000)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      COMMON/jtzdat/itmin,itmax,jtmin,jtmax
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      DIMENSION twobd(5),  ipba(max_conf), ans(5),dg1(5)
      COMMON/iterat/iter
      COMMON/renorm_v/ index(200000),vmtx(200000,5)
      COMMON/find_jpt/istart(0:2,0:30,0:1),iend(0:2,0:30,0:1)
      DIMENSION i_temp(200000),temp_v(200000,5),
     &          temp_istart(0:2,0:30,0:1),temp_iend(0:2,0:30,0:1)


      LOGICAL end_iter
      
      end_iter=.FALSE.
      iter=1
      DO WHILE((.NOT.end_iter).AND.(iter.LE.5))
         num=0
         dif=0.D0
         sum=0.0
         DO itz=itmin,itmax,1
            DO jtot=jtmin,jtmax
               DO ipar=0,1
                  CALL pauli(jtot,itz,ipar,nconf,ipba)
                  nn_conf=0
                  DO ibra=1, nconf
                     ja=ipba(2*ibra-1)
                     jb=ipba(2*ibra)                     
                     DO iket=ibra, nconf
                        jc=ipba(2*iket-1)
                        jd=ipba(2*iket)
                        IF(iter.GE.2) THEN
                        CALL renorm(ja,jb,jc,jd,jtot,itz,ipar,ans)
                        ENDIF
                        CALL gettwo(itz,ja,jb,jc,jd,jtot,twobd)
                        CALL pphhmtx(ja,jb,jc,jd,jtot,itz,ipar,dg1)
                        fnorm=1.d0/(dij(ja,jb)*dij(jc,jd))
                        num=num+1
                        nn_conf=nn_conf+1
                        iab=100*ja+jb
                        icd=100*jc+jd
                        IF(iab.GE.icd) THEN
                           jconf=icd*10000+iab
                        ELSE
                           jconf=iab*10000+icd
                        ENDIF
                        i_temp(num)=jconf
                        DO ie=1,5
                           temp_v(num,ie)=twobd(ie)
                           ans(ie)=ans(ie)*fnorm
                           dg1(ie)=dg1(ie)*fnorm
                        ENDDO
                        WRITE(8,1000)itz,ipar,jtot,jconf,
     &                               dg1(1),ans(1),temp_v(num,1)

C     &                              ( vmtx(num,ie),ie=1,5)           
                        sum=sum+twobd(1)**2
                        dif=dif+(ans(1)-twobd(1))**2
                     ENDDO
                  ENDDO
                  IF(nconf.NE.0) THEN
                    temp_istart(itz+1,jtot,ipar)=num+1-nn_conf
                    temp_iend(itz+1,jtot,ipar)=num
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         DO i_tot=1,num
            index(i_tot)=i_temp(i_tot)
            DO ie=1,5
               vmtx(i_tot,ie)=temp_v(i_tot,ie)
            ENDDO
         ENDDO
         DO itz=itmin,itmax,1
            DO jtot=jtmin,jtmax
               DO ipar=0,1
                  istart(itz+1,jtot,ipar)=temp_istart(itz+1,jtot,ipar)
                  iend(itz+1,jtot,ipar)=temp_iend(itz+1,jtot,ipar)
               ENDDO
            ENDDO
         ENDDO
         test=dif/sum
         iter=iter+1
         IF(test.LE.0.0001) end_iter=.TRUE.
         IF((iter.EQ.5).AND.(.NOT.end_iter))WRITE(6,*) 'itermax =5'
         WRITE(8,*) 'End iteration number:', iter, test
      ENDDO

 1000 FORMAT(1X,I2,1X,I2,1X,I2,1X,I10,5X,5F9.5)

      END






c
c
c    *******************************************************************
c    *               selects all two-body contributions                *
c    *******************************************************************
c
c
      SUBROUTINE gettwo(itz,ja,jb,jc,jd,jtot,twobd)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/jtzdat/itmin,itmax,jtmin,jtmax
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION twobd(5),dg1(5),dg2(5),dg4(5)

      COMMON/print_opt/print
      LOGICAL  print

 1000 FORMAT(/2x,17htwobody diagrams ,6i4/)
 1001 FORMAT(/2x,21htwobody contribution=,11f8.4/)
 1002 FORMAT(2x,22h1st-order contribution,2x,11f8.4)
 1003 FORMAT(2x,20hcore-polarization   ,11f8.4)
 1005 FORMAT(2x,20h2nd-order pp-ladder ,11f8.4)

      DO i=1,5
	 twobd(i)=0.0d0
      ENDDO
      ipar=MOD(ll(ja)+ll(jb),2)
      fnorm=1.d0/(dij(ja,jb)*dij(jc,jd))
      CALL pphhmtx(ja,jb,jc,jd,jtot,itz,ipar,dg1)
      CALL select3(itz,ja,jb,jc,jd,jtot,dg2,iorder)
c      CALL diagpp(ja,jb,jc,jd,jtot,dg4)
      DO i=1,5
	 twobd(i)=(dg1(i)+dg2(i)+dg4(i))*fnorm
      ENDDO
      IF(print)THEN
         WRITE(6,1000)itz, ja, jb, jc, jd, jtot
   	 WRITE(6,1002)(dg1(i)*fnorm,i=1,5)
	 WRITE(6,1003)(dg2(i)*fnorm,i=1,5)
	 WRITE(6,1005)(dg4(i)*fnorm,i=1,5)
         WRITE(6,1001)(twobd(i),i=1,5)
      ENDIF

      END



c
c
c    *******************************************************************
c               selects the contributions from the screened 
C               core-pol diagram
c    *******************************************************************
c
c
      SUBROUTINE select3(itz,ja,jb,jc,jd,jtot,dg2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION dg2(5),d2a(5),d2b(5),d2c(5),d2d(5),dia2(5)

      iab=(jj(ja)+jj(jb))/2-jtot+1
      icd=(jj(jc)+jj(jd))/2-jtot+1
      sgabjt=iph(iab)      
      sgcdjt=iph(icd)      
      sgabcd=sgcdjt*sgabjt

      DO i=1,5
         dg2(i)=0.d0
	 d2a(i)=0.d0
         d2b(i)=0.d0
         d2c(i)=0.d0
         d2d(i)=0.d0
      ENDDO
      CALL diag2(ja,jb,jc,jd,jtot,dia2)
      DO i=1,5
         d2a(i)=dia2(i)
      ENDDO
      IF(jc.EQ.jd) THEN
         DO i=1,5
            d2b(i)=sgcdjt*d2a(i)
         ENDDO
      ELSE
         CALL diag2(ja,jb,jd,jc,jtot,dia2)
         DO i=1,5
            d2b(i)=sgcdjt*dia2(i)
         ENDDO
      ENDIF
      IF(ja.EQ.jb) THEN
         DO i=1,5
            d2c(i)=sgabjt*d2b(i)
         ENDDO
      ELSE
         CALL diag2(jb,ja,jd,jc,jtot,dia2)
	 DO i=1,5
            d2c(i)=sgabcd*dia2(i)
         ENDDO
      ENDIF
      IF(ja.EQ.jb)THEN
         DO i=1,5
            d2d(i)=sgabjt*d2a(i)
         ENDDO
      ELSE
         CALL diag2(jb,ja,jc,jd,jtot,dia2)
         DO i=1,5
            d2d(i)=sgabjt*dia2(i)
         ENDDO
      ENDIF
      DO i=1,5
         dg2(i)=d2a(i)+d2b(i)+d2c(i)+d2d(i)
      ENDDO

      END


c
c
c
c    *******************************************************************
c    *                  core-polarization diagram                      *
c    *******************************************************************
c
c
      SUBROUTINE diag2(ja,jb,jc,jd,jtot,dia2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_conf=12000, max_dim=250)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION dia2(5),ans1(5),ans2(5),ipba(max_conf),
     &          f_mtx(5,max_dim,max_dim)
      iad=(jj(ja)+jj(jd))/2+jtot
      DO i=1,5
         dia2(i)=0.d0
      ENDDO
      j1min=IABS((jj(ja)-jj(jc))/2)
      j1max=(jj(ja)+jj(jc))/2
      j2min=IABS((jj(jb)-jj(jd))/2)
      j2max=(jj(jb)+jj(jd))/2
      jphmin=MAX0(j1min,j2min)
      jphmax=MIN0(j1max,j2max)
      IF(jphmin.GT.jphmax)RETURN
      iparph=MOD(ll(ja)+ll(jc),2)
      IF(MOD(ll(jd)+ll(jb),2).NE.iparph) RETURN
      DO 10 jph=jphmin,jphmax
         sixj=sjs(jj(jc),jj(jd),2*jtot,jj(jb),jj(ja),2*jph)
         CALL ph_config(jph,npd,jd,jb,ipba)
         CALL mtx_f(jph,jd,jb,npd,ipba,f_mtx)
         DO 20 nconf1=1,npd
            jh1=ipba(2*nconf1-1)
            jp1=ipba(2*nconf1)
            CALL crossc(ja,jh1,jc,jp1,jph,ans1)
            IF(ans1(2).EQ.0.d0) GOTO 20
            DO 30 nconf2=1,npd
               jh2=ipba(2*nconf2-1)
               jp2=ipba(2*nconf2)
               CALL crossc(jp2,jb,jh2,jd,jph,ans2)
               IF(ans2(2).EQ.0.d0) GOTO 30
               iphx=(jj(jh2)+jj(jp2))/2
               sg=iph(iad+iphx+jj(jh2))*sixj
               DO 40 ie=1,5
                  dia2(ie)=dia2(ie)+sg*ans1(ie)*ans2(ie)
     &                    * f_mtx(ie,nconf2,nconf1)
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE

      END

c
c
c    *******************************************************************
c    *                      second-order pp ladder                     *
c    *******************************************************************
c
c
      SUBROUTINE diagpp(ja,jb,jc,jd,jtot,dpp)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/iterat/iter
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      COMMON/jtzdat/itmin,itmax,jtmin,jtmax
      PARAMETER (max_conf=12000, max_dim=250)
      DIMENSION dpp(5),ans1(5),ans2(5),ipba(max_conf),
     &          gr_mtx(5,max_dim,max_dim)

      DO i=1,5
         dpp(i)=0.d0
      ENDDO
      itz=(itzp(ja)+itzp(jb))/2
      ipar=MOD(ll(jb)+ll(ja),2)
      CALL pp_config(jtot,npd,jc,jd,ipar,itz,ipba)
      CALL mtx_g(jtot,jc,jd,npd,ipba,gr_mtx)
      DO 10 nc1=1,npd
         jp1=ipba(2*nc1-1)
         jp2=ipba(2*nc1)
         IF(iter.EQ.1) THEN
            CALL pphhmtx(ja,jb,jp1,jp2,jtot,itz,ipar,ans1)
         ELSE
            CALL renorm(ja,jb,jp1,jp2,jtot,itz,ipar,ans1)
         ENDIF
         IF(ans1(2).EQ.0.d0) GOTO 10
         DO 20 nc2=1,npd
            jp3=ipba(2*nc2-1)
            jp4=ipba(2*nc2)
            IF(iter.EQ.1) THEN
               CALL pphhmtx(jp3,jp4,jc,jd,jtot,itz,ipar,ans2)
            ELSE
               CALL renorm(jp3,jp4,jc,jd,jtot,itz,ipar,ans2)
            ENDIF
            IF(ans2(2).EQ.0.d0) GOTO 20
            DO l=1,5
               dpp(l)=dpp(l)+0.5*ans1(l)*ans2(l)*gr_mtx(l,nc2,nc1)
            ENDDO
 20      CONTINUE
 10   CONTINUE

      END       






      SUBROUTINE  ph_config(jph,npd,jd,jb,ipba)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_conf=12000)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      LOGICAl triag, dencheck, ptest
      DIMENSION ipba(max_conf)
      npd=0
      DO 10 jh=1,nn1
         jjh=jj(jh)
         IF(ptest(jh,nn1)) GOTO 10
         DO 20 jp=np1,nn2
            IF(ptest(jp,nn2)) GOTO 20
            jjp=jj(jp)
            IF(triag(2*jph,jjp,jjh)) GOTO 20
            IF(MOD(ll(jd)+ll(jb),2).NE.MOD(ll(jp)+ll(jh),2)) GOTO 20
            nshell1=nshell(jh)
            nshell2=nshell(jp)
            idiff=nshell1-nshell2
            IF(dencheck(idiff)) GOTO 20
            npd=npd+1
            k2=npd*2
            k1=k2-1
            ipba(k1)=jh
            ipba(k2)=jp
 20      CONTINUE
 10   CONTINUE

      END



      SUBROUTINE  pp_config(jph,npd,jc,jd,ipar,itz,ipba)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_conf=12000)
      COMMON/ipndat/nn1,nn2,np1,np2,mpf,mpl,mnf,mnl
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      LOGICAl triag, dencheck, ptest
      DIMENSION ipba(max_conf)
      npd=0
      DO 10 jh=np1,nn2
         jjh=jj(jh)
         IF(ptest(jh,nn2)) GOTO 10
         DO 20 jp=np1,nn2
            IF(ptest(jp,nn2)) GOTO 20
            jjp=jj(jp)
            IF(triag(2*jph,jjp,jjh)) GOTO 20
            IF(ipar.NE.MOD(ll(jp)+ll(jh),2)) GOTO 20
            IF(2*itz.NE.(itzp(jp)+itzp(jh))) GOTO 20
            nshell1=nshell(jh)+nshell(jp)
            nshell2=nshell(jd)+nshell(jc)
            idiff=nshell1-nshell2
            IF(dencheck(idiff)) GOTO 20
            npd=npd+1
            k2=npd*2
            k1=k2-1
            ipba(k1)=jh
            ipba(k2)=jp
 20      CONTINUE
 10   CONTINUE

      END







      SUBROUTINE mtx_f(jph,jd,jb,npd,ipba,f_mtx)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER( max_conf=12000,max_dim=250, istart=5)
      COMMON/iterat/iter
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION ipba(max_conf),f_mtx(istart,max_dim,max_dim),
     &          amtx(max_dim,max_dim),ans(istart)
      DO ie=1,istart
         DO nc1=1,npd
            jh2=ipba(2*nc1-1)
            jp2=ipba(2*nc1)
            DO nc2=1,npd
               jh1=ipba(2*nc2-1)
               jp1=ipba(2*nc2)
               CALL crossc(jp1,jh2,jh1,jp2,jph,ans)
               sg=iph(jj(jh1)+(jj(jh1)+jj(jp1))/2-jph)
               amtx(nc2,nc1)=sg*ans(ie)
            ENDDO
            amtx(nc1,nc1)=amtx(nc1,nc1)+(e(jh2)-e(jp2))*iter
         ENDDO
         CALL matinv(amtx,npd,max_dim)
         DO nc1=1,npd
            DO nc2=1,npd
               f_mtx(ie,nc2,nc1)=amtx(nc2,nc1)
            ENDDO
         ENDDO
      ENDDO

      END




      SUBROUTINE mtx_g(jph,jd,jb,npd,ipba,f_mtx)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER( max_conf=12000,max_dim=250, istart=5)
      COMMON/iterat/iter
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION ipba(max_conf),f_mtx(istart,max_dim,max_dim),
     &          amtx(max_dim,max_dim),ans(istart)
      DO ie=1,istart
         DO nc1=1,npd
            jp1=ipba(2*nc1-1)
            jp2=ipba(2*nc1)
            itz=(itzp(jp1)+itzp(jp2))/2
            ipar=MOD(ll(jp1)+ll(jp2),2)
            DO nc2=1,npd
               jp3=ipba(2*nc2-1)
               jp4=ipba(2*nc2)
               IF(iter.EQ.1) THEN
                  CALL pphhmtx(jp3,jp4,jp2,jp1,jtot,itz,ipar,ans)
               ELSE
                  CALL renorm(jp3,jp4,jp2,jp1,jtot,itz,ipar,ans)
               ENDIF
               amtx(nc2,nc1)=-ans(ie)*0.5
            ENDDO
            amtx(nc1,nc1)=amtx(nc1,nc1)-e(jp1)-e(jp2)+a(ie)
         ENDDO
         CALL matinv(amtx,npd,max_dim)
         DO nc1=1,npd
            DO nc2=1,npd
               f_mtx(ie,nc2,nc1)=amtx(nc2,nc1)
            ENDDO
         ENDDO
      ENDDO

      END




C     *********************************************************
C            Routines to do mtx inversion, from Numerical
C            Recepies, Teukolsky et al. Routines included
C            below are MATINV, LUDCMP and LUBKSB. See chap 2
C            of Numerical Recepies for further details
C            Note well the max dim set to nmax
C     *********************************************************


      SUBROUTINE matinv(a,n,np)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(nmax=250)
      DIMENSION a(np,np),y(nmax,nmax),indx(nmax)

      DO i=1,n
         DO j=1,n
            y(i,j)=0.
         ENDDO
         y(i,i)=1.
      ENDDO
      CALL  ludcmp(a,n,np,indx,d)
      DO j=1,n
         call lubksb(a,n,np,indx,y(1,j))
      ENDDO
      DO i=1,n
         DO j=1,n
            a(i,j)=y(i,j)
         ENDDO
      ENDDO

      END
 



C
C      **************************************************************
C      *   CALCULATES THE CROSS-COUPLED  MATRIX ELEMENTS            *
C      **************************************************************
C


      SUBROUTINE crossc(ja,jb,jc,jd,jtot,cross)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/iterat/iter
      PARAMETER (istart=5)    
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION cross(istart),ans(istart)

      DO i=1,istart
         cross(i)=0.0d0
      ENDDO
      fnorm=DSQRT(2.d0*jtot+1.d0)
      IF(iph(ll(ja)+ll(jb)).NE.iph(ll(jc)+ll(jd))) RETURN
      IF((itzp(ja)+itzp(jb)).NE.(itzp(jc)+itzp(jd))) RETURN
      itz=(itzp(ja)+itzp(jb))/2
      ipar=MOD(ll(ja)+ll(jb),2)
      iad=(jj(ja)+jj(jd))/2+jtot
      jbra_min=IABS(jj(ja)-jj(jb))
      jbra_max=jj(ja)+jj(jb)
      jket_min=IABS(jj(jc)-jj(jd))
      jket_max=jj(jc)+jj(jd)
      jmax=MIN0(jbra_max,jket_max)/2
      jmin=MAX0(jbra_min,jket_min)/2
      IF(jmin.GT.jmax)RETURN

      DO 10 jt=jmin,jmax
         IF((ja.EQ.jb).AND.(MOD(jt,2).NE.0)) GOTO 10
         IF((jc.EQ.jd).AND.(MOD(jt,2).NE.0)) GOTO 10      
         sgbdjt=DFLOAT(iph(iad+jt))
         xlshat=(2.d0*jt+1.d0)*sgbdjt*fnorm
         sixj=sjs(jj(jc),jj(ja),2*jtot,jj(jb),jj(jd),2*jt)
         IF(iter.EQ.1) THEN
            CALL pphhmtx(ja,jb,jc,jd,jt,itz,ipar,ans)
         ELSE
            CALL renorm(ja,jb,jc,jd,jt,itz,ipar,ans)
         ENDIF
         DO l=1,5
            cross(l)=cross(l)+ans(l)*sixj*xlshat
         ENDDO
 10   CONTINUE

      END




C            ************************************************************
C            *             SUBROUTINE TO READ IN THE G-MATRIX           *
C            ************************************************************


      SUBROUTINE readme
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION pmtx(5)
      LOGICAL tester
      INTEGER it,ip,jt,t_temp,p_temp,j_temp,i_count,i_energy
      INTEGER nummtxel
      INTEGER*4 iconf,index
      COMMON/one/ index(200000),vmtx(200000,5)
      COMMON/locateg/istart(2,2,20),iend(2,2,20)
 1000 FORMAT(I10)
 1001 FORMAT(1X,I1,1X,I1,1X,I2,1X,I10,5X,5F9.5)
      DO i=1,2,1
         DO j=1,2,1
            DO l=1,20,1
               istart(i,j,l)=0 
               iend(i,j,l)=0
            ENDDO
         ENDDO
      ENDDO
      READ(7,1000)nummtxel
      istart(1,1,1)=1
      DO i_count=1,nummtxel,1
         READ(7,1001)it,ip,jt,iconf,(pmtx(i),i=1,5)
         it=it+1
         ip=ip+1
         jt=jt/2+1
         index(i_count)=iconf
         DO i_energy=1,5
            vmtx(i_count,i_energy)=pmtx(i_energy)
         ENDDO
         IF(i_count.EQ.nummtxel)THEN
           iend(it,ip,jt)=i_count
         ENDIF
         tester=(t_temp.NE.it).OR.(p_temp.NE.ip).OR.(j_temp.NE.jt)
         IF(tester.AND.i_counT.GE.2)THEN
           istart(it,ip,jt)=i_count
           iend(t_temp,p_temp,j_temp)=i_count-1
         ENDIF
         t_temp=it
         p_temp=ip
         j_temp=jt
      ENDDO

      END

C                *******************************************************
C                *     ROUTINE TO SET UP THE PN G-MATRIX FROM          *
C                *     THAT IN THE ISOSPIN FORMALISM                   *
C                *******************************************************

      SUBROUTINE pphhmtx(ja,jb,jc,jd,jt,itz,ipar,gmtpn)
      IMPLICIT REAL*8 (A-H,O-Z)    
      LOGICAL triag
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      DIMENSION gmtpn(5), gmtx0(5), gmtx1(5)
      DO i=1,5
         gmtpn(i)=0.d0
         gmtx0(i)=0.d0
         gmtx1(i)=0.d0
      ENDDO
      IF((ja.EQ.jb).AND.(MOD(jt,2).NE.0)) RETURN       
      IF((jc.EQ.jd).AND.(MOD(jt,2).NE.0)) RETURN       
      IF(triag(jj(ja),jj(jb),2*jt)) RETURN
      IF(triag(jj(jc),jj(jd),2*jt)) RETURN
      fbra=1.D0/DSQRT(2.D0)
      fket=1.D0/DSQRT(2.D0)
      ia=ja/2+MOD(ja,2)
      ib=jb/2+MOD(jb,2)
      ic=jc/2+MOD(jc,2)
      id=jd/2+MOD(jd,2)

      IF(IABS(itz).EQ.1) THEN
        CALL findg(1,ia,ib,ic,id,jt,gmtpn)
      ELSEIF(itz.EQ.0) THEN
        IF(itzp(ja).LT.itzp(jb)) fbra=-fbra
        IF(itzp(jc).LT.itzp(jd)) fket=-fket
        CALL findg(0,ia,ib,ic,id,jt,gmtx0)
        CALL findg(1,ia,ib,ic,id,jt,gmtx1)
        DO ie=1,5
           gmtpn(ie)=0.5d0*gmtx1(ie)+ 
     :               gmtx0(ie)*fbra*fket
        ENDDO
      ENDIF

      END





C                *******************************************************
C                *             G-MTX SEARCH ROUTINE                    *
C                *             I   IB.GE.IA                            *
C                *             II  ID.GE.IC                            *
C                *             III IB.GE.ID                            *
C                *             IV  IB.EQ.ID THEN IA.GE.IC              *
C                *     RETURNS UNNORMALIZED MTX EL                     *
C                *******************************************************


      SUBROUTINE FINDG(ITT,JAA,JBB,JCC,JDD,JT,GMTXEL)
      IMPLICIT REAL*8 (A-H,O-Z)    
      INTEGER*4 ICONF,INDEX
      COMMON/LOCATEG/ISTART(2,2,20),IEND(2,2,20)
      COMMON /NLJA/NA(200),LA(200),JA(200),NHBAR(153)
      COMMON/ONE/INDEX(200000),VMTX(200000,5)
      LOGICAL SEARCHER
      DIMENSION GMTXEL(5)
 1000 FORMAT(2X,13HNO MATRIX FOR,2X,2I4,2X,I11)
      DO 10 I=1,5
         GMTXEL(I)=0.0D0
 10   CONTINUE
      FNORM=DIJ(JAA,JBB)*DIJ(JCC,JDD)
      IPAR=MOD((LA(JAA)+LA(JBB)),2)
      ICHECK=-IPH(ITT+JT)
      IF((JAA.EQ.JBB).AND.(ICHECK.LT.0).OR.(JCC.EQ.JDD).AND.
     :   (ICHECK.LT.0)) RETURN
      IF(JAA.LE.JBB) THEN
         KAB=2
         IA=JAA
         IB=JBB
      ELSE
         IA=JBB
         IB=JAA
         KAB=(JA(JAA)+JA(JBB))/2 +JT+ITT
      ENDIF
      IF(JCC.LE.JDD) THEN
         IC=JCC
         ID=JDD
         KCD=2
      ELSE
         IC=JDD
         ID=JCC
         KCD=(JA(JCC)+JA(JDD))/2+JT+ITT
      ENDIF
      IF((ID.LT.IB).OR.((ID.EQ.IB).AND.(IC.LT.IA))) THEN
         ICD=100*IA+IB
         IAB=100*IC+ID
      ELSE
         IAB=100*IA+IB
         ICD=100*IC+ID
      ENDIF
      ICONF=10000*ICD+IAB
      LOW=ISTART(ITT+1,IPAR+1,JT+1)
      IHIGH=IEND(ITT+1,IPAR+1,JT+1)
      SEARCHER=.TRUE.
      IS=0
      IFIND=0

C     *************************************************************
C     I            BINARY SEARCH FOR MATRIX ELEMENTS              I
C     *************************************************************

       DO WHILE(SEARCHER)
          MID=(LOW+IHIGH)/2
          IS=IS+1
          IF(IS.GE.100)THEN
            SEARCHER=.FALSE.
          ELSEIF(ICONF-INDEX(MID).LT.0)THEN
            IHIGH=MID-1
            SEARCHER=.TRUE.
          ELSEIF(ICONF-INDEX(MID).GT.0)THEN
            LOW=MID+1
            SEARCHER=.TRUE.
          ELSEIF(ICONF-INDEX(MID).EQ.0)THEN
            IFIND=MID
            SEARCHER=.FALSE.
          ENDIF
       ENDDO
       IF(IFIND.EQ.0)THEN
          RETURN
       ENDIF
       DO 20 I=1,5
          GMTXEL(I)=VMTX(IFIND,I)*DFLOAT(IPH(KAB+KCD))*FNORM
 20    CONTINUE

       END



C                *******************************************************
C                *     RETURNS UNNORMALIZED MTX EL                     *
C                *******************************************************


      SUBROUTINE renorm(ja,jb,jc,jd,jt,itz,ipar,gmtxel)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/renorm_v/ index(200000),vmtx(200000,5)
      COMMON/find_jpt/istart(0:2,0:30,0:1),iend(0:2,0:30,0:1)
      COMMON /nlj/nn(200),ll(200),jj(200),itzp(200),nshell(200)
      LOGICAL searcher
      DIMENSION gmtxel(5)

      DO i=1,5
         gmtxel(i)=0.0D0
      ENDDO
      fnorm=dij(ja,jb)*dij(jc,jd)
      IF(IABS(itz).EQ.1) THEN
         IF(ja.LE.jb) THEN
            ia=ja
            ib=jb
            kab=2
         ELSE
            ia=jb
            ib=ja
            kab=(jj(ja)+jj(jb))/2-jt+1
         ENDIF
         IF(jc.LE.jd) THEN
            ic=jc
            id=jd
            kcd=2
         ELSE
            ic=jd
            id=jc
            kcd=(jj(jd)+jj(jc))/2-jt+1
         ENDIF
      ELSEIF(itz.EQ.0) THEN
         IF(itzp(ja).LT.itzp(jb)) THEN
            ia=ja
            ib=jb
            kab=2
         ELSEIF(itzp(ja).GT.itzp(jb)) THEN
            ia=jb
            ib=ja
            kab=(jj(ja)+jj(jb))/2-jt+1
         ENDIF
         IF(itzp(jc).LT.itzp(jd)) THEN 
            ic=jc
            id=jd
            kcd=2
         ELSEIF(itzp(jc).GT.itzp(jd)) THEN
            ic=jd
            id=jc
            kcd=(jj(jd)+jj(jc))/2-jt+1
         ENDIF
      ENDIF
      iab=100*ia+ib
      icd=100*ic+id
      IF(iab.LE.icd) THEN
         iconf=10000*iab+icd
      ELSE
         iconf=10000*icd+iab
      ENDIF
      low=istart(itz+1,jt,ipar)
      ihigh=iend(itz+1,jt,ipar)
      searcher=.TRUE.
      is=0
      ifind=0
      DO WHILE(searcher)
         mid=(low+ihigh)/2
         is=is+1
         IF(is.GE.100)THEN
            searcher=.FALSE.
         ELSEIF(iconf-index(mid).LT.0)THEN
            ihigh=mid-1
            searcher=.TRUE.
         ELSEIF(iconf-index(mid).GT.0)THEN
            low=mid+1
            searcher=.TRUE.
         ELSEIF(iconf-index(mid).EQ.0)THEN
            ifind=mid
            searcher=.FALSE.
         ENDIF
      ENDDO
      IF(ifind.EQ.0)THEN
         RETURN
      ENDIF
      DO i=1,5
         gmtxel(i)=vmtx(ifind,i)*DFLOAT(iph(kab+kcd))*fnorm
      ENDDO
      
      END




C
C             *********************************************
C             *    FUNCTION TO CHECK # OSC. EXCITATIONS   *
C             *********************************************
C

      LOGICAL FUNCTION dencheck(i)
      IMPLICIT REAL*8 (A-H,O-Z)    
      COMMON/energy/a(5),wstart,wcn(11),e(200),ihomega
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
 

