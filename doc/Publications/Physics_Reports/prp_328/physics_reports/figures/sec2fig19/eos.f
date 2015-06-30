           implicit real*8(a-h,o-z)
           parameter (e=15.8, s=32., gamma=0.6)
           parameter (delta=0., rho0=0.16)

           a0=2./3/acos(-1.)/acos(-1.)*((4.*s/197.327)**3)
           a0=a0/rho0
           do i=1,200
              rho=0.01*i
              x=rho/rho0
              p=1./3./a0/x
              q=-1./a0/2./x
c              fu=-q+sqrt(q**2+p**3)
c              fv=-q-sqrt(q**2+p**3)
c              u=fu/abs(fu)*(abs(fu)**0.333333333)
c              v=fv/abs(fv)*(abs(fv)**0.333333333)
c              chi=u+v
              phi=atan(sqrt(p**3)/q)
              psi=-(abs(tan(phi/2.)))**0.33333333
              psi=atan(psi)
              chi=-2.*sqrt(p)/tan(2.*psi)
              xp=(1.-chi)*0.5
              erg=e*x*(x-2.-delta)/(1.+delta*x)+
     :            s*(x**gamma)*(chi**2)
              WRITE(6,'(E12.6,2x,e12.6,2x,e12.6)')rho,erg
           enddo
           end
