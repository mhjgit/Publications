         implicit real*8(a-h,o-z)
          do i=1,500
             read(5,*) a,b,c,d,e
             write(6,*) a, abs(b-c)
          enddo
          end