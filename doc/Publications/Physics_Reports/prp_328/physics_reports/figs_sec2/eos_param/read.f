          implicit real*8(a-h,o-z)
          open(unit=7,file='v18pnm.dat')
          open(unit=8,file='v18snm.dat')

          write(*,*) 'number of data points:'
          read(*,*) np
          read(7,*) nnn, lll
          read(8,*) nnn, lll
          do n=1,np
             read(7,*) dens, xpnm
             read(8,*) dens, xsnm
             write(6,*) dens, xpnm-xsnm
          enddo
          end
