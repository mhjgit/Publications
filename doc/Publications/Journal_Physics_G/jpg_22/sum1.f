      implicit real*8(a-h,o-z)
      sum=0.0d0
      sum1=0.0d0
      sum0=0.0d0
      sumt0=0.d0
      sumt=0.d0
      sumt1=0.d0
      i=0
      i0=0
      i1=0
      x=1000.D0
    1 read(5,100,end=2) n11,n12,n21,n22,n31,n32,n41,n42,it
     1 ,ij,a 
  100 format(10i3,f9.4)
      if(n11.ne.n31.or.n21.ne.n41) goto 1
      if(n12.ne.n32.or.n22.ne.n42) goto 1 
      i=i+1
      sum=sum+((iabs(it)+1)*(ij+1)*a)
      sumt=sumt+(iabs(it)+1)*(ij+1)
      if(it.eq.0) then
        sum0=sum0+((iabs(it)+1)*(ij+1)*a)
        sumt0=sumt0+(iabs(it)+1)*(ij+1)
        i0=i0+1
      elseif(iabs(it).eq.2) then
        sum1=sum1+((iabs(it)+1)*(ij+1)*a)
        sumt1=sumt1+(iabs(it)+1)*(ij+1)
        i1=i1+1
      endif
      goto 1
    2 continue
      write(*,*) i, i0, i1      
      if(i0.gt.0) THEN

      write(*,*)x*sum/sumt/DFLOAT(i),x*sum0/sumt0/DFLOAT(i0),
     &          x*sum1/sumt1/DFLOAT(i1)
      else
      write(*,*)x*sum/sumt/DFLOAT(i),x*sum1/sumt1/DFLOAT(i1)
      endif
      end


