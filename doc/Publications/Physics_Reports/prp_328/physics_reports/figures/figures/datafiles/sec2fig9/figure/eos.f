         implicit real*8(a-h,o-z)
         read(5,*) n1
         do i=1, n1
            read(5,*) dens, e, p, xk, vc, gamma
            write(6,*) dens, p
            write(7,*) dens, xk
            write(8,*) dens, vc
            write(9,*) dens, gamma
          enddo
          end

