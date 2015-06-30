      implicit real*8(a-h,o-z)
      sum=0.0d0
      sum1=0.0d0
      sum0=0.0d0
      i=0
    1 read(5,100,end=2) n11,n12,n21,n22,n31,n32,n41,n42,it
     1 ,ij,a 
  100 format(10i3,f9.4)
      i=i+1
      sum=sum+a
      if(it.eq.0) sum0=sum0+a
      if(iabs(it).eq.2) sum1=sum1+a
      goto 1
    2 continue
      write(*,*) i,sum,sum0,sum1

      end



