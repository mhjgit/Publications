             implicit real*8(a-h,o-z)
             open(unit=6,file='bhf_n.dat')
             open(unit=7,file='bhf_p.dat')
             open(unit=8,file='bhf_e.dat')
             open(unit=9,file='bhf_mu.dat')
             open(unit=10,file='bhf_sm.dat')
             open(unit=11,file='bhf_l.dat')

             do i=1,120
                read(5,*) dens, x1, x2, x3, x4, x5, x6, x7
                write(6,'(E12.6,2x,e12.6)') dens, x1
                write(7,'(E12.6,2x,e12.6)') dens, x2
                write(8,'(E12.6,2x,e12.6)') dens, x3
                write(9,'(E12.6,2x,e12.6)') dens, x4
                write(10,'(E12.6,2x,e12.6)') dens, x5
                write(11,'(E12.6,2x,e12.6)') dens, x6
             enddo
             end

