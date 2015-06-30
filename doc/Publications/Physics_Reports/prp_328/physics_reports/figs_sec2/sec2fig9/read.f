         implicit real*8(a-h,o-z)
          pi = 3.1415926535897932
         open(unit=6,file='ebonn.dat')
         open(unit=7,file='pbonn.dat')
         open(unit=8,file='csbonn.dat')
         open(unit=9,file='gammabonn.dat')
         do i=1, 55
            read(5,*) a,b,c,d,e,f,g
            write(6,'(2e12.4)') a, b-938.926*a
            write(7,'(2e12.4)') a,d
            write(8,'(2e12.4)') a, f
            write(9,'(2e12.4)') a, g
          enddo
          end
