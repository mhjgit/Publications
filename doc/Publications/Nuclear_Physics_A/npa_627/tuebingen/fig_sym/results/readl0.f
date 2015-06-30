      program read_files
      implicit none
      integer i
      real*8 a(100), b(100), c(100), d(100) 
      read(5,*)
      read(5,*) (a(i),i=1,19)
      read(5,*)
      read(5,*) (b(i),i=1,19)
      read(5,*)
      read(5,*) (c(i),i=1,19)
      read(5,*)
      read(5,*) (d(i),i=1,19)
      do i=1,19
         write(6,1000) a(i), d(i)-b(i)-c(i)
      enddo
 1000 format(2e12.6)
      end
      

