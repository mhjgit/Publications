       implicit real*8 (a-h,o-z)
       dimension dens(60), e1(60),e2(60), e3(60)


        read(5,*) (dens(i),i=1,19)  
        read(5,*) (e1(i),i=1,19)
        read(5,*) (e2(i),i=1,19)

       do i=1,19
          write(6,*) dens(i) , e1(i)
c           write(6,*) dens(i) , e1(i)-e2(i)-e3(i)
       enddo
       end



