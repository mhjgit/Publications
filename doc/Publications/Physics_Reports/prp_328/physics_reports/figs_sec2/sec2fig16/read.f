         implicit real*8(a-h,o-z)
         do i=1, 19
            read(5,*) xkf,e
             rho=(xkf**3)*2./3./acos(-1.)/acos(-1.)
             write(6,*) rho,e
          enddo
          end

