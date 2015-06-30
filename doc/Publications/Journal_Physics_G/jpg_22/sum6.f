      implicit real*8(a-h,o-z)
      character*100 file1
      sum=0.0d0
      sum1=0.0d0
      sum0=0.0d0
      sumt0=0.d0
      sumt=0.d0
      sumt1=0.d0

      i=0
      write(*,*) ' quantum numbers:'
      read(*,*) n, j
      write(*,*) 'data-file'
      read(*,*) file1
      open(unit=5,file=file1)
    1 read(5,100,end=2) na,ja,nb,jb,nc,jc,nd,jd,it
     1 ,ij,a 
  100 format(10i3,f9.4)
      if((na.ne.n).or.(nb.ne.n).or.(nc.ne.n).or.(nd.ne.n)) goto 1
      if((ja.ne.j).or.(jb.ne.j).or.(jc.ne.j).or.(jd.ne.j)) goto 1
      i=i+1
      sum=sum+(iabs(it)+1)*(ij+1)*a
      sumt=sumt+(iabs(it)+1)*(ij+1)
      if(it.eq.0) then
        sum0=sum0+(iabs(it)+1)*(ij+1)*a
        sumt0=sumt0+(iabs(it)+1)*(ij+1)
      elseif(iabs(it).eq.2) then
        sum1=sum1+(iabs(it)+1)*(ij+1)*a
        sumt1=sumt1+(iabs(it)+1)*(ij+1)
      endif
      goto 1
    2 continue
      write(*,*) i,sum/sumt,sum1/sumt1
      end

