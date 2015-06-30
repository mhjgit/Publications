             implicit real*8(a-h,o-z)
             do i=1,120
                read(5,*) dens, x1, x2, x3, x4
                write(6,'(E12.6,2x,e12.6)') dens, x4
             enddo
             end

