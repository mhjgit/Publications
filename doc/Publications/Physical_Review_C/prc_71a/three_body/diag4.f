      program diagonalize
      implicit double precision (a-h,o-z)
      dimension arrin(4,4),eval(4),wk(11),xqw(13,2,2),xqwp(13,2,2)
     &,fac(14),arrin3(3,3),eval3(3),wk3(8)
      common/qv/xq(15,2,2),veff(2,2),heff(2,2),eigval(2)
      character*1 value,both,upper,lower
      data value,both,upper,lower/'N','V','U','L'/
      fac(1)=1.0d0
      do 34 i=2,14
 34   fac(i)=fac(i-1)*dfloat(i)
      rt3=dsqrt(3.d0)
      v=-0.4d0
      w=-0.4d0
      u=-0.15d0
      write(14,102)v,w,u
 102  format(3h v=,f7.3,4h  w=,f7.3,4h  u=,f7.3)
      write(14,104)
      arrin(1,1)=-1.5d0
      arrin(1,2)=rt3*u
      arrin(1,3)=rt3*v
      arrin(1,4)=0.d0
      arrin(2,2)=-0.5d0 +w*2.d0
      arrin(2,3)=u*2.d0
      arrin(2,4)=rt3*v
      arrin(3,3)=0.5d0 +w*2.d0
      arrin(3,4)=rt3*u
      arrin(4,4)=1.5d0
      do 2 i=1,3
      do 2 j=i+1,4
 2    arrin(j,i)=arrin(i,j)
      do 3 i=1,4
 3    write(14,105)(arrin(i,j),j=1,4)
      write(14,104)
c diagonalizes  4x4 array arrin gives eval & evectors in arrin
c      call dsyev(both,upper,4,arrin,4,eval,wk,11,infout)
      do 1 i=1,4
 1    write(14,100)eval(i),(arrin(i,j),j=1,4)
      write(14,101)infout
c******************************************************************
c              Effective interaction 2+3-body
c******************************************************************
      write(14,104)
      write(14,107)
      do 10 i=1,13
      xom=-0.0175d0+dfloat(i)*0.0025d0
      xqw(i,1,1)=3.d0*v*v/(xom-2.d0*(1.d0+w)-3.d0*u*u/(xom-3.d0))
      xqw(i,2,1)=rt3*u
     &+2.d0*rt3*u*v/(xom-2.d0*(1.d0+w)-3.d0*u*u/(xom-3.d0))
     &+3.d0*rt3*u*v*v/((xom-3.d0)*(xom-2.d0*(1.d0+w))-3.d0*u*u)
      xqw(i,1,2)=xqw(i,2,1)
      xqw(i,2,2)=1.d0+2.d0*w
     &+4.d0*u*u/(xom-2.d0*(1.d0+w)-3.d0*u*u/(xom-3.d0))
     &+3.d0*v*v/(xom-3.d0-3.d0*u*u/(xom-2.d0*(1.d0+w)))
     &+12.d0*v*u*u/((xom-2.d0*(1.d0+w))*(xom-3.d0)-3.d0*u*u)
 10   continue
      do 12 j=1,2
      do 12 k=1,2
 12   xq(1,j,k)=xqw(7,j,k)
c     derivatives of the Q-box
      do 11 i=1,7
      n=12-i
      do 13 k=1,2
      do 13 m=1,2
      do 14 j=i,n
 14   xqwp(j+1,k,m)=(xqw(j+2,k,m)-xqw(j,k,m))/0.005d0
      do 15 j=i,n
 15   xqw(j+1,k,m)=xqwp(j+1,k,m)
      xq(i+1,k,m)=xqw(7,k,m)/fac(i)
 13   continue
 11   continue
c     0 folds
      do 16 i=1,2
      do 16 j=1,2
 16   heff(i,j)=xq(1,i,j)
      heff(1,1)=heff(1,1)-1.5d0
      heff(2,2)=heff(2,2)-1.5d0
      call roots
      write(14,106)0,eigval(1),eigval(2)
c     1-6 folds
      do 17 i=1,7
      do 18 j=1,2
      do 18 k=1,2
 18   veff(j,k)=xq(1,j,k)
      call folds(i)
      heff(1,1)=heff(1,1)-1.5d0
      heff(2,2)=heff(2,2)-1.5d0
      call roots
 17   write(14,106)i,eigval(1),eigval(2)
c******************************************************************
c              Effective interaction 2-body; test 2-body
c******************************************************************
      write(14,104)
      write(14,108)
      rt2=dsqrt(2.d0)
      arrin3(1,1)=-1.0d0
      arrin3(1,2)=u/rt2
      arrin3(1,3)=v
      arrin3(2,1)=u/rt2
      arrin3(2,2)=w
      arrin3(2,3)=u/rt2
      arrin3(3,1)=v
      arrin3(3,2)=u/rt2
      arrin3(3,3)=1.0d0
      do 23 i=1,3
 23   write(14,105)(arrin3(i,j),j=1,3)
      write(14,104)
c diagonalizes  3x3 array arrin gives eval & evectors in arrin
c      call dsyev(both,upper,3,arrin3,3,eval3,wk3,8,infout)
      do 21 i=1,3
 21   write(14,100)eval3(i),(arrin3(i,j),j=1,3)
      write(14,101)infout
      write(14,104)
      do 30 i=1,15
      z=2.d0**i
      xq(i,1,1)=-v*v/z
      xq(i,2,1)=-v*u/(rt2*z)
      xq(i,1,2)=xq(i,2,1)
      xq(i,2,2)=-u*u/(2.d0*z)
 30   continue
      xq(1,1,2)=u/rt2+xq(1,1,2)
      xq(1,2,1)=u/rt2+xq(1,2,1)
      xq(1,2,2)=1.d0+w+xq(1,2,2)
c     0 folds
      do 26 i=1,2
      do 26 j=1,2
 26   heff(i,j)=xq(1,i,j)
      heff(1,1)=heff(1,1)-1.0d0
      heff(2,2)=heff(2,2)-1.0d0
      call roots
      write(14,106)0,eigval(1),eigval(2)
c     1-14 folds
      do 27 i=1,14
      do 28 j=1,2
      do 28 k=1,2
 28   veff(j,k)=xq(1,j,k)
      call folds(i)
      heff(1,1)=heff(1,1)-1.0d0
      heff(2,2)=heff(2,2)-1.0d0
      call roots
 27   write(14,106)i,eigval(1),eigval(2)
c******************************************************************
c              Effective interaction 2-body applied to 3-body
c******************************************************************
      write(14,104)
      write(14,109)
      do 40 i=1,15
      z=2.d0**i
      y=3.d0**I
      xq(i,1,1)=-3.d0*v*v/z
      xq(i,2,1)=-rt3*v*u/z
      xq(i,1,2)=xq(i,2,1)
      xq(i,2,2)=-u*u/z -v*v/y
 40   continue
      xq(1,1,2)=u*rt3+xq(1,1,2)
      xq(1,2,1)=u*rt3+xq(1,2,1)
      xq(1,2,2)=1.d0+2.d0*w+xq(1,2,2)
c     0 folds
      do 46 i=1,2
      do 46 j=1,2
 46   heff(i,j)=xq(1,i,j)
      heff(1,1)=heff(1,1)-1.5d0
      heff(2,2)=heff(2,2)-1.5d0
      call roots
      write(14,106)0,eigval(1),eigval(2)
c     1-14 folds
      do 47 i=1,14
      do 48 j=1,2
      do 48 k=1,2
 48   veff(j,k)=xq(1,j,k)
      call folds(i)
      heff(1,1)=heff(1,1)-1.5d0
      heff(2,2)=heff(2,2)-1.5d0
      call roots
 47   write(14,106)i,eigval(1),eigval(2)
 100  format(1h ,5(1pe15.6))
 105  format(1h ,4(1pe15.6))
 101  format(1h ,i3)
 104  format(1h )
 106  format(7h folds=,i2,2(1pe15.6))
 107  format(22h full eff int 2+3 body)
 108  format(40h two body eff int test for two body case)
 109  format(38h 2 body eff int applied to 3 body case)
      stop
      end
c**************************************************************
      subroutine roots
c     roots of a 2x2 matrix
      implicit double precision (a-h,o-z)
      common/qv/xq(15,2,2),veff(2,2),heff(2,2),eigval(2)
      x=(heff(1,1)-heff(2,2))**2+4.d0*heff(1,2)*heff(2,1)
      if (x.lt.0.d0) then
      eigval(1)=0.d0
      eigval(2)=0.d0
      write(14,100)((heff(i,j),j=1,2),i=1,2)
 100  format(14h complex roots,4(1pe15.6))
      else
      x=dsqrt(x)
      y=heff(1,1)+heff(2,2)
      eigval(1)=(y-x)/2.d0
      eigval(2)=(y+x)/2.d0
      endif
      return 
      end
c**************************************************************
      subroutine folds(n)
      implicit double precision (a-h,o-z)
      common/qv/xq(15,2,2),veff(2,2),heff(2,2),eigval(2)
      dimension hold(2,2),pwv(14,2,2)
c     powers of Q-box
 6    do 3 k=1,n
      do 10 i=1,2
      do 10 j=1,2
      if (k.eq.1)then
      pwv(1,i,j)=veff(i,j)
      else 
      pwv(k,i,j)=0.d0
      do 2 m=1,2
 2    pwv(k,i,j)=pwv(k,i,j)+pwv(k-1,i,m)*veff(m,j)
      endif
 10   continue
 3    continue
c     effective interaction with folds
      do 12 i=1,2
      do 12 j=1,2
      hold(i,j)=xq(1,i,j)
      do 1 k=1,n
      do 1 m=1,2
 1    hold(i,j)=hold(i,j)+xq(k+1,i,m)*pwv(k,m,j)
 12   continue
      y=0.d0
      do 5 i=1,2
      do 5 j=1,2
      y=y+dabs(veff(i,j)-hold(i,j))
      heff(i,j)=hold(i,j)
      

 5    veff(i,j)=hold(i,j)
      if(y.gt.1.d-7)go to 6
      do i=1,2
         write(14,'(2(4X,E12.6))') (heff(j,i),j = 1, 2)
      enddo  
      return
      end
