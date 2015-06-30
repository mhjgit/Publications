             implicit real*8(a-h,o-z)
             xsp=0.
             xl=0.
             xsm=0.
             xs0=0.
             pi=acos(-1.)
             do i=1, 33
                read(5,*) dens, frac, x1, x2, x3, x4, x5, x6, x7,
     &                    x8, x9, x10, x11, x12, x13, x14, x15
                xp=((pi**2)*dens*frac*3)**(1./3.)
                xn=((pi**2)*dens*(1.-frac)*3)**(1./3.)
                write(6,'(22(E12.6,2x))') dens, xn,xp,xsm,xl,xs0,xsp, 
     &                     x1, x2, x3, x4, x5, x6, x7,
     &                    x8, x9, x10, x11, x12, x13, x14, x15

             enddo
             end

