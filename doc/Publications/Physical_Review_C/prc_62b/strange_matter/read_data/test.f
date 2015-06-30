             implicit real*8(a-h,o-z)
             dimension xnn(12), xll(12), xpp(12), xss(12)
             do i=1, 12
                read(5,*) dens, xnn(i),xpp(i),xss(i),xll(i)

             enddo
             write(6,'(3(E12.6,2x))') (xnn(i),i=1,12) 
            write(6,'(3(E12.6,2x))') (xpp(i),i=1,12) 
            write(6,'(3(E12.6,2x))') (xss(i),i=1,12) 
            write(6,'(3(E12.6,2x))') (xll(i),i=1,12) 
             end
