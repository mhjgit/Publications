       function pot(r)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      radial part of the potential                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit real*8 (a-h,o-z)
        parameter (nwavv0=100)
        dimension v(8)
        dimension rx(nwavv0),px(nwavv0,31) 
        common/potsel/ipot,np
c        vc(l)=v(1)
c        vs(l)=v(3)
c        vt(l)=v(5)
c        vls(l)=v(7)
c        vct(l)=v(2)
c        vst(l)=v(4)
c        vtt(l)=v(6)
c        vlst(l)=v(8)
c         reid 8
c         lpot=1
c         call potreid(lpot,r,v)
c         argonne 8
          las=1
          lae=1
         rx(1)=r
c         call arpot14(las,lae,rx,px) 
c         do i=1,8
c         v(i)=px(1,i)
c         end do
         call av8pop(r,v)
         if(np.eq.1) pot=v(1)-v(2)-v(3)+v(4)
         if(np.eq.2) pot=-4.d0*v(4)
         if(np.eq.3) pot=2.d0*v(3)-2.d0*v(4)
         if(np.eq.4) pot=-2.d0*v(2)+2.d0*v(4)
         if(np.eq.5) pot=v(7)
         if(np.eq.6) pot=v(8)
         if(np.eq.7) pot=v(5)
         if(np.eq.8) pot=v(6)
         end


      subroutine arpot14(las,lae,rx,px) 
c
c     Argonne V6 potential: Use only the first six terms
c
      implicit real*8(a-h,o-z)
      parameter (nwavv0=100)
      parameter ( one=1.e0,two=2.e0,three=3.e0,four=4.e0,five=5.e0,     
     &  six=6.e0,seven=7.e0,eight=8.e0,xnine=9.e0,ten=10.e0,            
     &  eleven=11.e0,twelve=12.e0,zero=0.e0,                            
     &  half=0.5e0,third=1.e0/3.e0,fourth=0.25e0,fifth=0.2e0,           
     &  sixth=1.e0/6.e0,tenth=0.1e0,twelfth=1.e0/12.e0 )                
      parameter (ione=1,itwo=2) 
      dimension rx(nwavv0),px(nwavv0,31) 
      dimension x(nwavv0),rcut(nwavv0),ypi(nwavv0),tpi(nwavv0),         
     & tpi2(nwavv0),w(nwavv0)                                           
      do 10 la=las,lae 
      x(la)=138.03e0*rx(la)/197.33e0 
      www=two*rx(la)*rx(la)
      if(www.gt.300.d0) then
      rcut(la)=one
      else
      rcut(la)=one-exp(-two*rx(la)*rx(la)) 
      endif
      ypi(la)=exp(-x(la))*rcut(la)/x(la) 
      tpi(la)=(one+three/x(la)+three/x(la)**2)*ypi(la)*rcut(la) 
      tpi2(la)=tpi(la)*tpi(la) 
      w(la)=one/(one+exp((rx(la)-.5e0)/.2e0)) 
   10 continue 
      do 20 la=las,lae 
      px(la,1)=-4.801125e0*tpi2(la)+2061.5625e0*w(la)
      px(la,2)=.798925e0*tpi2(la)-477.3125e0*w(la)
      px(la,3)=1.189325e0*tpi2(la)-502.3125e0*w(la)
      px(la,4)=3.72681e0*ypi(la)+.182875e0*tpi2(la)+97.0625e0*w(la)
      px(la,5)=(-.1575e0*tpi2(la)+108.75e0*w(la))
      px(la,6)=(3.72681e0*tpi(la)-.7525e0*tpi2(la)+297.25e0*w(la))
      px(la,7)=(.5625e0*tpi2(la)-719.75e0*w(la))
      px(la,8)=(.0475e0*tpi2(la)-159.25e0*w(la))
      px(la,9)=0.e0 
      px(la,10)=0.e0 
      px(la,11)=0.e0 
      px(la,12)=0.e0 
      px(la,13)=0.e0 
      px(la,14)=0.e0 
   20 continue 
      return 
      END                                           

      subroutine arpot8(las,lae,rx,px) 
c
c     Argonne V8 potential
c
      implicit real*8(a-h,o-z)
      parameter (nwavv0=100)
      parameter ( one=1.e0,two=2.e0,three=3.e0,four=4.e0,five=5.e0,     
     &  six=6.e0,seven=7.e0,eight=8.e0,xnine=9.e0,ten=10.e0,            
     &  eleven=11.e0,twelve=12.e0,zero=0.e0,                            
     &  half=0.5e0,third=1.e0/3.e0,fourth=0.25e0,fifth=0.2e0,           
     &  sixth=1.e0/6.e0,tenth=0.1e0,twelfth=1.e0/12.e0 )                
      parameter (ione=1,itwo=2) 
      dimension rx(nwavv0),px(nwavv0,31) 
      dimension x(nwavv0),rcut(nwavv0),ypi(nwavv0),tpi(nwavv0),         
     & tpi2(nwavv0),w(nwavv0)                                           
c
      do 10 la=las,lae 
      x(la)=138.03e0*rx(la)/197.33e0 
      rcut(la)=one-exp(-two*rx(la)*rx(la)) 
      ypi(la)=exp(-x(la))*rcut(la)/x(la) 
      tpi(la)=(one+three/x(la)+three/x(la)**2)*ypi(la)*rcut(la) 
      tpi2(la)=tpi(la)*tpi(la) 
      w(la)=one/(one+exp((rx(la)-.5e0)/.2e0)) 
   10 continue 
      do 20 la=las,lae 
      px(la,1)=                 -5.263625*tpi2(la)+2415.9375  *w(la) 
      px(la,2)=                  0.541425*tpi2(la)- 298.6875  *w(la) 
      px(la,3)=                  0.931825*tpi2(la)- 323.6875  *w(la) 
      px(la,4)=3.72681e0*ypi(la)+0.200375*tpi2(la)+  96.104167*w(la) 
      px(la,5)=                  0.011250*tpi2(la)-  58.75    *w(la) 
      px(la,6)=3.72681e0*tpi(la)-0.696250*tpi2(la)+ 241.41667 *w(la) 
      px(la,7)=                  0.937500*tpi2(la)- 942.75    *w(la) 
      px(la,8)=                 -0.057500*tpi2(la)- 204.25    *w(la) 
      px(la,9)=0.e0 
      px(la,10)=0.e0 
      px(la,11)=0.e0 
      px(la,12)=0.e0 
      px(la,13)=0.e0 
      px(la,14)=0.e0 
   20 continue 
      return 
      END                                           



c *id* av8op **********************************************************
c subroutine for strong interaction part of argonne v8' potential 
c in operator format
c calls subroutine consts
c ----------------------------------------------------------------------
c arguments for av8ppot
c r:   separation in fm
c vnn: output potential in MeV (8 component array)
c ----------------------------------------------------------------------
c order of operators l in vnn(l):
c l:    1=1                              2=t1.t2
c       3=s1.s2                          4=(s1.s2)(t1.t2)
c       5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
c       7=L.S                            8=L.S(t1.t2)
c where s1=sigma_1, t1=tau_1, etc.
c ----------------------------------------------------------------------
      subroutine av8pop(r,vnn)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension vnn(8)
      real*8 mpi0,mpic,mp,mn
      real*8 mpi,mu0,muc,mu
      data small/1e-4/
      do 5 l=1,8
        vnn(l)=0
    5 continue
      call consts(hc,mpi0,mpic,mp,mn,alpha)
      mpi=(mpi0+2.*mpic)/3.
      mu0=mpi0/hc
      muc=mpic/hc
      mu=mpi/hc
      fsq=.075
      cpi=2.1
      rws=.5
      aiws=5.
      x=mu*r
      x0=mu0*r
      xc=muc*r
      if (r.le.small) then
        tpi=3*cpi**2*r/mu**3
        ypi0=(mpi0/mpic)**2*(mpi0/3)*cpi*r/mu0
        tpi0=3*cpi*ypi0/mu0**2
        ypic=(mpic/3)*cpi*r/muc
        tpic=3*cpi*ypic/muc**2
      else
        rcut=1-exp(-cpi*r*r)
        ypi=exp(-x)*rcut/x
        tpi=(1+(3+3/x)/x)*ypi*rcut
        ypi0=(mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        tpi0=(1+(3+3/x0)/x0)*ypi0*rcut
        ypic=(mpic/3)*exp(-xc)*rcut/xc
        tpic=(1+(3+3/xc)/xc)*ypic*rcut
      end if
      ypi0=fsq*ypi0
      ypic=fsq*ypic
      tpi0=fsq*tpi0
      tpic=fsq*tpic
      tpi2=tpi*tpi
      ws=1/(1+exp((r-rws)*aiws))
      ws0=1/(1+exp(-rws*aiws))
      wsp=ws*(1+aiws*exp(-rws*aiws)*ws0*r)
      wsx=ws*x
      wsx2=wsx*x
      dypi00=(mpi0/mpic)**2*(mpi0/3)*cpi/mu0
      dypic0=(mpic/3)*cpi/muc
      ypi0p=ypi0-fsq*dypi00*ws*r/ws0
      ypicp=ypic-fsq*dypic0*ws*r/ws0
      ypi=(ypi0+2*ypic)/3
      tpi=(tpi0+2*tpic)/3
      p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p
      p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp
      p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p
      pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic
      pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      pls1=    -.62697*tpi2 -570.5571*wsp +819.1222*wsx2
      pl211=    .06709*tpi2 +342.0669*wsp -615.2339*wsx2
      pls21=    .74129*tpi2   +9.3418*wsp -376.4384*wsx2
      p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp
      pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic
      pls0=     .10180*tpi2  +86.0658*wsp -356.5175*wsx2
      pl210=   -.13201*tpi2 +253.4350*wsp   -1.0076*wsx2
      pls20=    .07357*tpi2 -217.5791*wsp  +18.3935*wsx2
      p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p
      p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
      p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p
      pl201=    .12472*tpi2  +16.7780*wsp
      p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp)
      pl200=   -.31452*tpi2 +217.4559*wsp
      p11=(p11pp+p11nn+p11np)/3
      p11cd=(.5*(p11pp+p11nn)-p11np)/6
      p11cs=(p11pp-p11nn)/4
      pt1=(pt1pp+pt1nn+pt1np)/3
      pt1cd=(.5*(pt1pp+pt1nn)-pt1np)/6
      pt1cs=(pt1pp-pt1nn)/4
      p01=(p01pp+p01nn+p01np)/3
      p01cd=(.5*(p01pp+p01nn)-p01np)/6
      p01cs=(p01pp-p01nn)/4
c -----------
      p00=p00+2*pl200
      pls0=pls0-2*pl210-3*pls20
      p11=p11+2*pl211+4*pls21/3
      pt1=pt1-5*pls21/12
      pls1=pls1-.5*pls21
c -----------
      vnn(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vnn(2)=.0625*(3*p11-3*p10  +p01-p00)
      vnn(3)=.0625*(3*p11  +p10-3*p01-p00)
      vnn(4)=.0625*(  p11  -p10  -p01+p00)
      vnn(5)=.25*(3*pt1+pt0)
      vnn(6)=.25*(  pt1-pt0)
      vnn(7)=.25*(3*pls1+pls0)
      vnn(8)=.25*(  pls1-pls0)
      return
      end
c *id* empot ***********************************************************
c subroutine for Coulomb part of Argonne v8' potential
c calls subroutine consts
c ----------------------------------------------------------------------
c arguments for empot
c r:    input separation in fm
c vem:  output potential in MeV
c ----------------------------------------------------------------------
c C1 = one-photon-exchange Coulomb with form factor
c ----------------------------------------------------------------------
      subroutine empot(r,vem)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 mpi0,mpic,mp,mn,mup,mun
      data small/1e-5/
      call consts(hc,mpi0,mpic,mp,mn,alpha)
      b=4.27
      br=b*r
      if (r.lt.small) then
        fcoulr=5*b/16
      else
        fcoulr=(1-(1+11*br/16+3*br**2/16+br**3/48)*exp(-br))/r
      end if
      vem=alpha*hc*fcoulr
      return
      end
c *id* consts **********************************************************
c subroutine for constants in av8' potential
c ----------------------------------------------------------------------
c arguments for consts
c hc:    output value for hbar*c (MeV-fm)
c mpi0:    "      "    "  neutral pion mass (MeV)
c mpic:    "      "    "  charged pion mass (MeV)
c mp:      "      "    "  proton mass (MeV)
c mn:      "      "    "  neutron mass (MeV)
c ----------------------------------------------------------------------
      subroutine consts(hc,mpi0,mpic,mp,mn,alpha)
      real*8 hc,mpi0,mpic,mp,mn
      hc=197.327053
      mpi0=134.9739
      mpic=139.5675
      mp=938.27231
      mn=939.56563
      alpha=1./137.035989
      return
      end
