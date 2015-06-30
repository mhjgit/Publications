       implicit real*8 (a-h,o-z)
       dimension dens(60), eos(60),prot(60)

       do i=1,56
        read(5,*) dens(i),prot(i),b,c,d,e,f,g,eos(i)
        write(6,*) dens(i),prot(i)
       enddo 

       end


