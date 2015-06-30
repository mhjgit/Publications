             implicit real*8(a-h,o-z)

             pi=acos(-1.)
             do i=1, 120
                read(5,*) dens, x1, x2, x3, x4, x5, x6, x7
                xn=((pi**2)*x1*3)**(1./3.)
                xp=((pi**2)*x2*3)**(1./3.)
                xs=((pi**2)*x5*3)**(1./3.)
                xl=((pi**2)*x6*3)**(1./3.)
                write(6,'(5(E12.6,2x))') dens, xn, xp,xs,xl 
             enddo
             end
