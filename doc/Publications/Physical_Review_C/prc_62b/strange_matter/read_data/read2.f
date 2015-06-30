             implicit real*8(a-h,o-z)
             xsp=0.
             xs0=0.
             do i=1,288
                read(5,*) dens, xn,xp,xsm,xl,x1, x2, x3, x4, x5, x6, x7,
     &                    x8, x9, x10, x11, x12, x13, x14, x15
                write(6,'(22(E12.6,2x))') dens, xn,xp,xsm,xl,xs0,xsp, 
     &                     x1, x2, x3, x4, x5, x6, x7,
     &                    x8, x9, x10, x11, x12, x13, x14, x15

             enddo
             end

