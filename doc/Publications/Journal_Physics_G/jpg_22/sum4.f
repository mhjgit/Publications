      implicit real*8(a-h,o-z)
      sum1=0.0d0
      sum0=0.0d0
      i=0
      read(*,*) jmx
      do j0=1,jmx+1
         sum0=sum0+(j0*2+1)
      enddo
      do j1=0,jmx
         sum1=sum1+3*(j1*2+1)
      enddo
      write(*,*) sum0+sum1,sum0,sum1
      end

