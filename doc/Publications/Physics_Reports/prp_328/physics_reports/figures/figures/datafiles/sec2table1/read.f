         implicit real*8(a-h,o-z)
          pi = 3.1415926535897932
         open(unit=6,file='n.dat')
         open(unit=7,file='p.dat')
         open(unit=8,file='e.dat')
         open(unit=9,file='mu.dat')
         read(5,*) 
         do i=1, 56
            read(5,*) rho,xp,un,up,ue,a,b,c,d
            xn= 1.-xp
            rhon=xn*rho
            rhop=xp*rho
            xp=(3.*pi*pi*rhop)**(1./3.)
            xn=(3.*pi*pi*rhon)**(1./3.)
            xe=ue/197.327
            if (ue**2 - 105.7**2 .ge. 0. ) then
               xmu=sqrt(ue**2 - 105.7**2)/197.327
            endif
            write(6,'(2e12.4)') rho, xn
            write(7,'(2e12.4)') rho, xp
            write(8,'(2e12.4)') rho, xe
            write(9,'(2e12.4)') rho, xmu
          enddo
          end
