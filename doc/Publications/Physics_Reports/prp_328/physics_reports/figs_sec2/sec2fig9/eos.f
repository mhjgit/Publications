         implicit real*8(a-h,o-z)
          pi = 3.1415926535897932
         read(5,*) n1, n2
        write(6,*) n1, n2
         do i=1, n1
            read(5,*) dens, e
            xp=(3.*pi*pi*dens)**(1./3.)
            xkin=3.*((xp*197.327)**2)/938.926/10.
            write(6,*) dens, e+xkin
          enddo
          end

