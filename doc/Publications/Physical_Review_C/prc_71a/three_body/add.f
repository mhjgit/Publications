      program diagonalize
      implicit double precision (a-h,o-z)
      dimension arrin(4,4),eval(4),wk(11),xqw(13,2,2),xqwp(13,2,2)
     &     ,fac(14),arrin3(3,3),eval3(3),wk3(8)
      common/qv/xq(15,2,2),veff(2,2),heff(2,2),eigval(2)


      fac(1)=1.0d0
      do 34 i=2,14
 34      fac(i)=fac(i-1)*dfloat(i)

         do j=1,2
            do k=1,2
               xq(1,j,k)=xqw(7,j,k)
            enddo
         enddo

         do i=1,6
            n=12-i
            do k=1,2
               do m=1,2
                  do j=i,n
                     xqwp(j+1,k,m)=(xqw(j+2,k,m)-xqw(j,k,m))/0.005d0
                  enddo
               enddo
            enddo
            do j=i,n
               xqw(j+1,k,m)=xqwp(j+1,k,m)
            enddo
            xq(i+1,k,m)=xqw(7,k,m)/fac(i)
         enddo

c     1-6 folds
      do i=1,6
         do j=1,2
            do k=1,2
               veff(j,k)=xq(1,j,k)
            enddo
         enddo
         call folds(i)
      enddo      















c******************************************************************
c     Effective interaction 2-body; test 2-body
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
 23               write(14,105)(arrin3(i,j),j=1,3)
                  write(14,104)
c     diagonalizes  3x3 array arrin gives eval & evectors in arrin
c     call dsyev(both,upper,3,arrin3,3,eval3,wk3,8,infout)
                  do 21 i=1,3
 21                  write(14,100)eval3(i),(arrin3(i,j),j=1,3)
                     write(14,101)infout
                     write(14,104)
                     do 30 i=1,15
                        z=2.d0**i
                        xq(i,1,1)=-v*v/z
                        xq(i,2,1)=-v*u/(rt2*z)
                        xq(i,1,2)=xq(i,2,1)
                        xq(i,2,2)=-u*u/(2.d0*z)
 30                  continue
                     xq(1,1,2)=u/rt2+xq(1,1,2)
                     xq(1,2,1)=u/rt2+xq(1,2,1)
                     xq(1,2,2)=1.d0+w+xq(1,2,2)
c     0 folds
                     do 26 i=1,2
                        do 26 j=1,2
 26                        heff(i,j)=xq(1,i,j)
                           heff(1,1)=heff(1,1)-1.0d0
                           heff(2,2)=heff(2,2)-1.0d0
                           call roots
                           write(14,106)0,eigval(1),eigval(2)
c     1-14 folds
                           do 27 i=1,14
                              do 28 j=1,2
                                 do 28 k=1,2
 28                                 veff(j,k)=xq(1,j,k)
                                    call folds(i)
                                    heff(1,1)=heff(1,1)-1.0d0
                                    heff(2,2)=heff(2,2)-1.0d0
                                    call roots
 27                                 write(14,106)i,eigval(1),eigval(2)
c******************************************************************
c     Effective interaction 2-body applied to 3-body
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
 40                                 continue
                                    xq(1,1,2)=u*rt3+xq(1,1,2)
                                    xq(1,2,1)=u*rt3+xq(1,2,1)
                                    xq(1,2,2)=1.d0+2.d0*w+xq(1,2,2)
c     0 folds
                                    do 46 i=1,2
                                       do 46 j=1,2
 46                                       heff(i,j)=xq(1,i,j)
                                          heff(1,1)=heff(1,1)-1.5d0
                                          heff(2,2)=heff(2,2)-1.5d0
                                          call roots
                                          write(14,106)0,eigval(1),eigval(2)
c     1-14 folds
                                          do 47 i=1,14
                                             do 48 j=1,2
                                                do 48 k=1,2
 48                                                veff(j,k)=xq(1,j,k)
                                                   call folds(i)
                                                   heff(1,1)=heff(1,1)-1.5d0
                                                   heff(2,2)=heff(2,2)-1.5d0
                                                   call roots
 47                                                write(14,106)i,eigval(1),eigval(2)
 100                                               format(1h ,5(1pe15.6))
 105                                               format(1h ,4(1pe15.6))
 101                                               format(1h ,i3)
 104                                               format(1h )
 106                                               format(7h folds=,i2,2(1pe15.6))
 107                                               format(22h full eff int 2+3 body)
 108                                               format(40h two body eff int test for two body case)
 109                                               format(38h 2 body eff int applied to 3 body case)
                                                   stop
                                                   end







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
      return
      end
