
   
      program thermo
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (max_dim=10000)
      DIMENSION eodd(max_dim), eeven(max_dim), ee(max_dim), yp(max_dim)
      DIMENSION sodd(max_dim), seven(max_dim)
      DIMENSION ss1(max_dim), ss2(max_dim)


 1000 FORMAT(6D12.4)
 1001 FORMAT(6D14.6)

      open(unit=7,file='t11100')
      open(unit=8,file='t12100')

      do i=1, 500
         ee(i)=i*0.1
      enddo


      nodd=1000
      do i=1, nodd
         read(7,*) xx, yy, zz,sodd(i),eodd(i)
      enddo

      neven=1000
      do i=1, neven
         read(8,*) xx,yy, zz,seven(i),eeven(i)
      enddo


      yp1=0.
      ypn=0.
      call spline(eodd,sodd,nodd,yp1,ypn,yp)      
      do i=1,500
         call splint(eodd,sodd,yp,nodd,ee(i),eee)
         ss1(i)=eee
      enddo



      yp1=0.
      ypn=0.
      call spline(eeven,seven,neven,yp1,ypn,yp)      
      do i=1,500
         call splint(eeven,seven,yp,neven,ee(i),eee)
         ss2(i)=eee
      enddo

      do i=1,500
         write(6,*) ee(i), ss2(i), ss1(i),  &
         ss2(i)/abs((ss2(i)-ss1(i))),ss1(i)/abs((ss1(i)-ss2(i)))
      enddo
      end


! takes as input x[1,..,n] and y[1,..,n] containing a tabulation
! y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) 
! together with yp_1 and yp2 for first derivatives  f(x) at x_0 
! and x_(n-1), respectively. Then the
! function returns y2[1,..,n] which contains the second 
! derivatives of f(x_i)at each point x_i. If yp1 and/or yp2 
! is larger than the constant INFINITY the function will 
! put corresponding second derivatives to zero.

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT NONE
      INTEGER :: i, k, n
      DOUBLE PRECISION, DIMENSION(n) :: x, y, y2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: u
      DOUBLE PRECISION :: p, qn, sig, un, ypn, yp1 

      ALLOCATE ( u (n) )
      IF (yp1 > .99E30) THEN
         y2(1)=0.
         u(1)=0.
      ELSE
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF
      DO i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO 
      IF (ypn > .99E30) THEN
         qn=0.
         un=0.
      ELSE
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      DO k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      ENDDO
      DEALLOCATE ( u )

      END SUBROUTINE spline


! takes xa[1,..,n] and y[1,..,n] which tabulates a function 
! (with the xa[i]'s in order) and given ya[0,..,n - 1], 
! which is the output from function spline() and with 
! given value of x returns a cubic--spline interpolation value y.

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      IMPLICIT NONE
      INTEGER :: k, n, klo, khi
      DOUBLE PRECISION, DIMENSION(n) :: xa, ya, y2a
      DOUBLE PRECISION ::  x, y, h, b, a

      klo=1
      khi=n
      DO WHILE  (khi-klo > 1)
         k=(khi+klo)/2
         IF(xa(k)  > x)THEN
            khi=k
         ELSE
            klo=k
         ENDIF
      ENDDO
      h=xa(khi)-xa(klo)
      IF (h == 0.) WRITE (6,*) 'Bad XA input in SPLINT.F.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
 
      END SUBROUTINE splint


