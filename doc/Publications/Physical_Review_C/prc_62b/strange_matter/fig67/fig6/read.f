             implicit real*8(a-h,o-z)
             do i=1,50
                read(5,*) x1, x2, x3, x4, x5, x6
                write(6,'(E12.6,2x,e12.6)') x2, x3
             enddo
             end

