           implicit real*8(a-h,o-z)
           parameter (e=15.8, s=30., gamma=0.5)
           parameter (delta=0.2, rho0=0.16)

           a0=0.095
           do i=1,200
              rho=0.01*i
              x=rho/rho0
              p=1./3./a0/x
              q=-1./a0/2./x
              fu=-q+sqrt(q**2+p**3)
              fv=-q-sqrt(q**2+p**3)
              u=fu/abs(fu)*(abs(fu)**0.333333333)
              v=fv/abs(fv)*(abs(fv)**0.333333333)
c              chi=u+v
               chi=0.
c              phi=atan(sqrt(p**3)/q)
c              psi=-(abs(tan(phi/2.)))**0.33333333
c              psi=atan(psi)
c              chi=-2.*sqrt(p)/tan(2.*psi)
c              xp=0.5*(1.-chi)
              erg=e*x*(x-2.-delta)/(1.+delta*x)+
     :            s*(x**gamma)*(chi**2)
              WRITE(6,'(E12.6,2x,e12.6,2x,e12.6)')rho,erg
C,xp
           enddo
           end
