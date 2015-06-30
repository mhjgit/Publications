         implicit real*8(a-h,o-z)
         open(7,file='s11c')
         open(8,file='s12c')
          do i=1,800
             read(7,*) a,b
             read(8,*) c,d
             write(6,*) a,abs(b-d)
          enddo
          end
