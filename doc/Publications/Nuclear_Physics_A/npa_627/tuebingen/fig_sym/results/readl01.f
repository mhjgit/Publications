      program read_files
      implicit none
      integer i
      real*8 a(100), b(100), c(100), d(100),e(100) ,f(100) 
     &       ,ff(100) ,g(100) ,gg(100) ,h(100) ,hh(100) 
      read(5,*)
      read(5,*) (a(i),i=1,19)
c      write(6,*) (a(i),i=1,19)
      read(5,*)
      read(5,*) (b(i),i=1,19)
c      write(6,*) (b(i),i=1,19)
      read(5,*)
      read(5,*) (c(i),i=1,19)
c      write(6,*) (c(i),i=1,19)
      read(5,*)
      read(5,*) (d(i),i=1,19)
c      write(6,*) (d(i),i=1,19)
      read(5,*)
      read(5,*) (e(i),i=1,19)
c      write(6,*) (e(i),i=1,19)
      read(5,*)
      read(5,*) (f(i),i=1,19)
c      write(6,*) (f(i),i=1,19)
      read(5,*)
      read(5,*) (ff(i),i=1,19)
c      write(6,*) (ff(i),i=1,19)
      read(5,*)
      read(5,*) (g(i),i=1,19)
c      write(6,*) (g(i),i=1,19)
      read(5,*)
      read(5,*) (gg(i),i=1,19)
c      write(6,*) (gg(i),i=1,19)
      read(5,*)
      read(5,*) (h(i),i=1,19)
c      write(6,*) (h(i),i=1,19)
      read(5,*)
      read(5,*) (hh(i),i=1,19)
c      write(6,*) (hh(i),i=1,19)
      do i=1,19
         write(6,1000) a(i),g(i)+gg(i)+h(i)+hh(i)-f(i)-ff(i) 
     &                 -e(i)-d(i)-b(i)-c(i)
      enddo
 1000 format(2e12.6)
      end
      

