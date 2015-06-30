      implicit real*8(a-h,o-z)
      sum=0.0d0
      sum1=0.0d0
      sum0=0.0d0
      sumt0=0.d0
      sumt=0.d0
      sumt1=0.d0
      i=0
    1 read(5,100,end=2) n11,n12,n21,n22,n31,n32,n41,n42,it
     1 ,ij,a 
  100 format(10i3,f9.4)
      i=i+1
      sum=sum+((iabs(it)+1)*(ij+1)*a)
      sumt=sumt+(iabs(it)+1)*(ij+1)
      if(it.eq.0) then
        sum0=sum0+((iabs(it)+1)*(ij+1)*a)
        sumt0=sumt0+(iabs(it)+1)*(ij+1)
      elseif(iabs(it).eq.2) then
        sum1=sum1+((iabs(it)+1)*(ij+1)*a)
        sumt1=sumt1+(iabs(it)+1)*(ij+1)
      endif
      goto 1
    2 continue
      write(*,*) i,sum/sumt,sum0/sumt0,sum1/sumt1
      write(*,*) i,sumt,sumt0,sumt1
      write(*,*) i,sum,sum0,sum1
      end


