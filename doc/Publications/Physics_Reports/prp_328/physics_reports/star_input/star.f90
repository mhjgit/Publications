!            Program to calculate the mass (SM), 
!            metric function, angular velocity, moment of inertia and 
!            radius (r) for a slowly rotating neutron star, with a  
!            user supplied EOS. 
!            For more details about the equations used see Hartle,
!            ApJ 150 (1967) p. 1005
!            Only monopole terms included here.
!            This is the Fortran 90 version.
!            Free format is employed, compile with e.g.,
!            f90 -c -free star.f and link with f90 -o prog.exe star.o
!
!
!            Coded by : M. Hjorth-Jensen            
!                       Dept. of Physics, Univ. of  
!                       Oslo, Norway, june/1991     
!                       based on white dwarf prog   
!                       whited.c of march 91
!            Language : Fortran 90                              
!            Last Fortran 90 update : 28/des/1997
!            Programs available at: http://www.uio.no/~mhjensen
!
!            You should have the energy per particle as function of 
!            particle density n, where n has units fm-3.
!            Note that the EOS you supply and which is read in
!            by SUBROUTINE eosfit MUST be in units of MeVfm^-3.
!            The dimensions used internally in the program are:
!            [rho] = MeVfm-3  energy per baryon 
!            user,
!            [n]   = 1/fm^3    particle density                     
!            [P]   = MeV/fm^3  pressure                       
!            [sm]  = Kg        mass of neutron star                 
!            [r]   = m         radius of neutron star           
!
!            Recall also that the constants defined 
!            below have the following dim           
!            [hbarc]=197.327 MeVfm
!            [g] = m^3Kg^-1s^-2 gravitational const 
!                = 6.67259*10^-11 m^3Kg^-1s^-2      
!                  or 6.67259 MeV^-2\hbarc/c^2
!            [sm_solar]=1.989D+30  kg
!            The solar mass is converted into MeV. Then
!            one has the relation
!                   1 kg=1.D+30/1.78266270D0  1 MeV/c^2
!            All taken from the Particle Data Group 
!            (1994)                                 
!            The program runs with dim-less variables
!            P=P/central_energy_density
!            rho=rho/central_energy_density
!            with central_energy_density the central energy-density of the star
!            r=r/r0 with r0=1.D+19 fm= 10km
!            m=m/m_solar
!            The dimension less constants used are
!            
!            const_1a=r0**3/sm_solar/speedoflight**2/central_energy_density
!            const_2=sm_solar*g/r0/speedoflight**2=const_2
!
!
!          Observables given by print_output:
!          Prints out central density in fm^-3, radius in km,
!          neutron star mass in terms of solar masses, angular velocity
!          at the surface (in units of s^-1 10km) and Moment of inertia
!          (in units of solar masses * (10km)^2)

!
!     Global constants used to render equations dim-less
!     All variables except const_1 (defined in main program
!     for every new value of the central density ) are kept
!     fixed during the execution of the program
      MODULE constants 
         DOUBLE PRECISION, PUBLIC :: const_1
!     10^19 fm = 10 km,typical value of R
         DOUBLE PRECISION, PUBLIC , PARAMETER ::  r0=1.E+19 
!     1 MeV/c^2=1.78266270D-30 kg
         DOUBLE PRECISION, PUBLIC, PARAMETER ::  xkg=1.E+30/1.78266270
!     \hbar  *  speed of light
         DOUBLE PRECISION, PUBLIC, PARAMETER :: hbarc=197.327
!     gravitational constant in fmMeV^-1
         DOUBLE PRECISION, PUBLIC, PARAMETER :: g=6.67259D-45*197.327
!     4*\pi                  
         DOUBLE PRECISION, PUBLIC, PARAMETER :: fourpi=12.56637061 
!     solar mass in MeV/c^2
         DOUBLE PRECISION, PUBLIC, PARAMETER :: sm_solar=1.989E+30*xkg    
!     dimless constant in mass equation           
         DOUBLE PRECISION, PUBLIC, PARAMETER :: c1=fourpi*(r0**3)/sm_solar
!     dimless constant in TOV equation         
         DOUBLE PRECISION, PUBLIC, PARAMETER :: const_2=sm_solar*g/r0
!     number of central densites chosen
         INTEGER, PUBLIC, PARAMETER :: number_central_densities=50
!     the differential equations of Hartle, excluding the quadropole terms
         INTEGER, PUBLIC , PARAMETER :: number_differential_eqs= 5
!     number of runge-kutta iterations
         INTEGER, PUBLIC , PARAMETER :: max_rk4_steps=100000
!     Dimensionless step in RK4 procedure
         DOUBLE PRECISION, PUBLIC, PARAMETER :: diff_eq_step = 0.001
!     Density from which one wants the central density to be calculated 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: first_density = 0.3
!     Last central density 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: last_density = 2.0
      END MODULE constants 


!
!     This module contains the parametrization of the EOS as
!     a polynomial in density. The number of terms kept in the
!     polynomial expansion is given by number_terms.
!     The only global variables which do change under the 
!     various runge-kutta iterations are 
!     pressure and central_energy_density
      MODULE eos
         USE constants
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: polynom_terms
         INTEGER :: number_terms 
         DOUBLE PRECISION :: pressure, central_energy_density
         CONTAINS
!
!     rho: energy per particle in units of MeV/fm^3
!     rho=\sum_{i=1}^{number of polynoms} a_i*density^{(i-1)/3}
!
            FUNCTION rho(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: rho
            DOUBLE PRECISION,  INTENT(IN) :: x
            INTEGER :: i
            rho=polynom_terms(1)
            DO i=2,number_terms
               rho=rho+polynom_terms(i)*(x**(FLOAT(i-1)/3.d0))
            ENDDO

            END FUNCTION rho
!
!     pressure in units of MeV/fm^3
!     pressure = density*d rho/d density - rho
! 
            FUNCTION press(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: press
            DOUBLE PRECISION, INTENT(IN)  :: x
            INTEGER :: i
            press=-rho(x)
            DO i=2,number_terms
               press=press+(FLOAT(i-1)/3.0*polynom_terms(i)*(x**(FLOAT(i-1)/3.)))
            ENDDO
            press=press-pressure
      
            END FUNCTION press

!
!     derivative dp/dn
!
            FUNCTION dpdn(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: dpdn
            DOUBLE PRECISION,  INTENT(IN) :: x
            INTEGER :: i
            dpdn=0.
            DO i=1,number_terms
               dpdn=dpdn+(((FLOAT(i-1)/3.0)**2)*polynom_terms(i)* &
               (x**(FLOAT(i-4)/3.)) &
               -FLOAT(i-1)/3.0*polynom_terms(i)*(x**(FLOAT(i-4)/3.)))
            ENDDO

            END FUNCTION dpdn
!
!     dP/dr, TOV equation, dimensionless
!
             FUNCTION tov(r,e_rho,y)
             IMPLICIT NONE
             DOUBLE PRECISION :: e_rho, tov, r 
             DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN)  :: y
             tov=-const_2*(e_rho+y(1))* &
                 (y(2)+const_1*(r**3)*y(1))/(r*r-2*const_2*r*y(2))
             END FUNCTION tov 
!
!     With given pressure y(1), finds the relevant energy density
!
             FUNCTION energy_density(y,x)
             IMPLICIT NONE
             DOUBLE PRECISION, PARAMETER :: tol=1.d-8
             DOUBLE PRECISION :: energy_density, x, zbrent
             DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN)  :: y
             pressure=y(1)*central_energy_density
             x=zbrent(0.D0,last_density,tol)
             energy_density=rho(x)/central_energy_density
             END FUNCTION  energy_density

      END MODULE eos 

!
!     Main program starts here
!
      PROGRAM  neutron_star
      USE eos
      USE constants 
      IMPLICIT NONE
      INTEGER :: j, i
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: cen_dens, &
                        radius, smass, inertia, omega
      DOUBLE PRECISION, DIMENSION (number_differential_eqs) :: ynew, y, ders
      DOUBLE PRECISION ::  star_radius, ync, density_step, ang_mom, e_rho,xx


!     Read in equation of state and make fit of it
      CALL eosfit
!     steps in choice of central densities
      density_step=(last_density-first_density)/number_central_densities   
      WRITE( 6,*) 'central density, radius, mass, Omega, I'
!     loop over central densities
      DO i=1, number_central_densities                         
!     Preparation of starting values
!     central density in fm^-3
         ync=first_density+density_step*i   
!     central energy density in units MeV/fm^3        
         central_energy_density=rho(ync)            
!     central pressure, in units of MeV/fm^-3
         pressure=0.
         pressure=press(ync)
!     dimless const_1
         const_1=c1*central_energy_density                  
!     Dimensionless pressure at the centre (function of central dens)
         y(1)=pressure/central_energy_density
!     Dimensionless start radius
         star_radius=diff_eq_step/10.
!     Dimensionless mass of the star at the centre 
         y(2)=const_1*(star_radius**3)/3.              
!     Metric function at the centre of the star
         y(3)=0.
!     Start values of the function w (y(5)) and u (y(4)) used in the calculation of
!     the angular velocity Omega and the derivative dudr of u
         y(4)=0.
         y(5)=0.1
!     Initial value of monopole function m_0
!         y(6)=fourpi/15.*const_1*(1.+y(1))*(2.+1/dpdn(ync))* &
!              (y(5)**2)*(star_radius**5)
!     Initial value of monople function p_0
!         y(7)=1./3.*(y(5)*star_radius)**2
         ders=0.
!     Start of Runge-Kutta solution to differential equations
!     number of RK4 steps given by variable j
         j=0                               
         DO WHILE((pressure > 0.).AND.(j <= max_rk4_steps))
            j=j+1
!     new dimless radius 
            star_radius=star_radius+diff_eq_step
!     start values for derivatives
            CALL derivatives(star_radius,y,ders)        
!     RK4 procedure
            CALL runge_kutta_4(star_radius,y,ynew,ders)
!     new values for y(1) (pressure), y(2) (mass) and y(3) (metric function)
            y=ynew                   
!     new pressure, now in units of MeV/fm^-3
            pressure=y(1)*central_energy_density
!     angular velocity in units of s^-1 * 10 km
            omega(i)=(y(5)+y(4)/3./(star_radius**3))
         ENDDO
         IF ((pressure > 0.).AND.(j >= max_rk4_steps)) THEN
            WRITE (6,*) ' no convergence for density:', ync
            WRITE (6,*) ' pressure and j:', pressure, j
         ENDIF
!     Calculate the angular momentum J, dimensionless
         ang_mom=(star_radius**4)*ders(5)/6./const_1/const_2
!     central density in units of fm^-3
!        cen_dens(i)=LOG10(ync/xkg*1.D+42*938.926)
         cen_dens(i)=ync                       
!     radius in units of 10 km
         radius(i)=star_radius*10.
!     mass in units of solar masses                      
         smass(i)=y(2)                        
!     moment of inertia in units of solar mass * (10km)^2, factor 100
         inertia(i)=(ang_mom/omega(i))*const_1*100.
         WRITE(6,'(E12.3,5E12.4)') cen_dens(i),radius(i),smass(i), omega(i),inertia(i)
      ENDDO          

      END  PROGRAM neutron_star
!
!        4th-Runge-Kutta solution of coupled equations     
!        See any textbook on numerical methods for details
!
      SUBROUTINE runge_kutta_4(x,y,yout,dydx)
      USE constants
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(number_differential_eqs) :: yt, dyt, dym
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN) :: y, dydx
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(OUT) :: yout
      DOUBLE PRECISION :: hh, h6, xh
      DOUBLE PRECISION, INTENT(IN) :: x 


      hh=diff_eq_step*0.5; h6=diff_eq_step/6. ; xh=x+hh
!     first rk-step
      yt=y+hh*dydx
      CALL derivatives(xh,yt,dyt)
!     second rk-step
      yt=y+hh*dyt
      CALL derivatives(xh,yt,dym)
!     third rk-step
      yt=y+diff_eq_step*dym;  dym=dyt+dym
      CALL derivatives(x+diff_eq_step,yt,dyt)
!     fourth rk-step
      yout=y+h6*(dydx+dyt+2.*dym)

      END SUBROUTINE runge_kutta_4
!
!     Here the expressions for the derivatives are set up with monopole terms only
! 
      SUBROUTINE derivatives(r,y,ders)
      USE eos
      USE constants 
      IMPLICIT NONE
      DOUBLE PRECISION :: e_rho, xx
      DOUBLE PRECISION, INTENT(IN)  :: r
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(OUT)  :: ders
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN)  :: y

      IF(y(1) > 0.) THEN
         e_rho=energy_density(y,xx)
!     TOV equation=dp/dr
         ders(1)=tov(r,e_rho,y)
!     dm/dr=const_1*r*r*energydensity, gravitational mass as function of r
         ders(2)=const_1*(r**2)*e_rho
!     d\nu/dr=-dp/dr*2/(e_rho+y(1)), metric function
         ders(3)=-2.*ders(1)/(e_rho+y(1))
!     derivative of the function u used to evaluate omega=y(5)
         ders(4)=(ders(3)*0.5-const_2*(r**2*y(2)-r*ders(2))/ &
                 (r**2-2.*const_2*r*y(2)))*(y(4)+4.*(r**3)*y(5))
!     derivative of angular velocity given by Omega-omega
         ders(5)=y(4)/r**4
!     derivative of the monopole function m_0
!         ders(6)=fourpi*const_1*r*r*((e_rho+y(1))**2)*y(7)/dpdn(xx) + &
!                  (r**4)*(ders(5)**2)*EXP(-y(3))*(1.-2*const_2*y(2)/r)/12./const_2 + &
!                 fourpi*2/3.*const_1*(r**4)*EXP(-y(3))*(e_rho+y(1))*(y(5)**2)
!     derivative of the monopole function p_0
!         ders(7)=1./3.*r*EXP(-y(3))*(y(5)**2)*(2*y(5)-r*y(5)*ders(3)+2*r*ders(5)) + &
!                 (r**3)*EXP(-y(3))*(ders(5)**2)/12. - &
!                 fourpi*const_1*const_2*r*(e_rho+y(1))*y(7)/(1.-2*const_2*y(2)/r) - &
!                 y(6)*const_2*(1.+2*fourpi*const_1*const_2*y(1)*r*r)/(r*r)/  &
!                 ((1.-2*const_2*y(2)/r)**2)
      ENDIF

      END SUBROUTINE derivatives
!     
!     This function uses Brent's method to find the root of a function
!     known to lie between x1 and x2. The root is returned as zbrent
!     and is refined till its accuracy is tol
!
      FUNCTION zbrent(x1,x2,tol)
      USE eos
      IMPLICIT NONE
      DOUBLE PRECISION :: tol, x1, x2, zbrent
      DOUBLE PRECISION, PARAMETER :: EPS= epsilon(x1)
      INTEGER, PARAMETER :: ITMAX=100
      INTEGER :: iter
      DOUBLE PRECISION :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

      a=x1
      b=x2
      fa=press(a)
      fb=press(b)
!      IF((fa > 0.0.AND.fb > 0.0).OR.(fa < 0.0.AND.fb < 0.0)) THEN
!          WRITE(*,*) 'root must be bracketed for zbrent' ; STOP 
!      ENDIF
      c=b
      fc=fb
      DO iter=1,ITMAX
         IF ((fb > 0.0.AND.fc > 0.).OR.(fb < 0.0.AND.fc < 0.0))THEN
            c=a
            fc=fa
            d=b-a
            e=d
         ENDIF
         IF(ABS(fc) < ABS(fb)) THEN
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         ENDIF
         tol1=2.*eps*ABS(b)+0.5*tol
         xm=.5*(c-b)
         IF((ABS(xm) <= tol1) .OR. (fb == 0.))THEN
            zbrent=b
            RETURN
         ENDIF
         IF((ABS(e) >= tol1) .AND.( ABS(fa) > ABS(fb)) ) THEN
            s=fb/fa
            IF(a == c) THEN
               p=2.*xm*s
               q=1.-s
            ELSE
               q=fa/fc
               r=fb/fc
               p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
               q=(q-1.)*(r-1.)*(s-1.)
            ENDIF
            IF(p > 0.) q=-q
            p=ABS(p)
            IF(2.*p < MIN(3.*xm*q-ABS(tol1*q),ABS(e*q))) THEN
               e=d
               d=p/q
            ELSE
               d=xm
               e=d
            ENDIF
         ELSE
            d=xm
            e=d
         ENDIF
         a=b
         fa=fb
         b=b+MERGE(d,SIGN(tol1,xm), ABS(d) > tol1 )
         fb=press(b)
      ENDDO
      WRITE (*,*) 'zbrent exceeding maximum iterations'; STOP
      zbrent=b

      END FUNCTION zbrent

!
!     This subroutine fits the equation of state in terms
!     of a polynomial expansion. Number of terms in the polynomial
!     expansion is given by the variable number_terms
!     The number of data ( density and corresponding energy
!     per particle) are defined by the variable number_of_data
!
      SUBROUTINE  eosfit
      USE constants
      USE eos
      IMPLICIT NONE
      INTEGER :: i, number_of_data
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: w, n, e,sig
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cvm, v, u
      DOUBLE PRECISION :: chisq, xp
      EXTERNAL eos_param

      READ (5,*) number_of_data, number_terms                   
      ALLOCATE (cvm(number_terms,number_terms), &
                v(number_terms,number_terms) )
      ALLOCATE ( u(number_of_data,number_terms) )
      ALLOCATE ( w(number_terms), polynom_terms(number_terms) )
      ALLOCATE (n(number_of_data), e(number_of_data), sig(number_of_data) ) 
      sig=1.
      DO i=1,number_of_data
         READ (5,*)n(i) , e(i)
!              e(i)=e(i)*n(i)+938.927*n(i)
      ENDDO
      CALL svdfit(n,e,sig,number_of_data,polynom_terms, &
                  number_terms, &
                  u,v,w,number_of_data,number_terms,chisq)
      CALL svdvar(v,number_terms,number_terms,w,cvm,number_terms)
      WRITE(6,*) ' Quality of energy fit:'
      WRITE(6,*) ' Number of terms in polynomial expansion:',  number_terms
      DO i=1,number_terms
         WRITE (6,*) i, polynom_terms(i)
      ENDDO
      WRITE (6,'(''CHISQ = '',F12.5)') chisq
      IF(chisq > 0.1) THEN
         WRITE(6,*) 'chi-square tolerance too high, try with more terms'
      ENDIF
      WRITE(6,*) 'Density in [fm^-3]'
      WRITE(6,*) 'Energy and pressure in [MeVfm^-3]'
      WRITE(6,*) 'Density, e-fitted, original e and pressure'
!      DO i=1,number_of_data
!         WRITE(6,'(E10.4,2X,E10.4,2X,E10.4,2X,E10.4)') n(i), rho(n(i)), e(i),press(n(i))
!      ENDDO

      DEALLOCATE (cvm, v)
      DEALLOCATE ( u )
      DEALLOCATE ( w )
      DEALLOCATE (n, e, sig ) 

      END SUBROUTINE eosfit


      SUBROUTINE eos_param(x,afunc,ma)
      IMPLICIT NONE
      INTEGER :: ma, i
      DOUBLE PRECISION :: afunc, x
      DIMENSION afunc(ma)
      afunc(1)=1.D0
      DO i=2,ma
         afunc(i)=afunc(i-1)*x**(1./3.)
      ENDDO

      END SUBROUTINE eos_param



      SUBROUTINE svdvar(v,ma,np,w,cvm,ncvm)
      IMPLICIT NONE
      INTEGER :: ma,ncvm,np,MMAX
      DOUBLE PRECISION ::  cvm(ncvm,ncvm),v(np,np),w(np)
      PARAMETER (MMAX=50)
      INTEGER :: i,j,k
      DOUBLE PRECISION :: sum,wti(MMAX)
      DO i=1,ma
        wti(i)=0.
        IF (w(i) /= 0.) wti(i)=1./(w(i)*w(i))
      ENDDO
      DO i=1,ma
         DO j=1,i
            sum=0.
            DO k=1,ma
               sum=sum+v(i,k)*v(j,k)*wti(k)
            ENDDO
            cvm(i,j)=sum
            cvm(j,i)=sum
         ENDDO
      ENDDO

      END SUBROUTINE svdvar


      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      IMPLICIT NONE
      INTEGER :: m,mp,n,np,NMAX
      DOUBLE PRECISION :: b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      INTEGER :: i,j,jj
      DOUBLE PRECISION :: s,tmp(NMAX)
      DO j=1,n
         s=0.
         IF (w(j) /= 0.) THEN
            DO i=1,m
               s=s+u(i,j)*b(i)
            ENDDO
            s=s/w(j)
         ENDIF
         tmp(j)=s
      ENDDO
      DO j=1,n
         s=0.
         DO jj=1,n
            s=s+v(j,jj)*tmp(jj)
         ENDDO
         x(j)=s
      ENDDO

      END SUBROUTINE svbksb


      SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq)
      IMPLICIT NONE
      INTEGER :: ma,mp,ndata,np,NMAX,MMAX
      DOUBLE PRECISION :: chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np), &
                          x(ndata),  y(ndata),TOL
      PARAMETER (NMAX=1000,MMAX=50,TOL=1.e-5)
      INTEGER :: i,j
      DOUBLE PRECISION :: sum,thresh,tmp,wmax,afunc(MMAX), &
                          b(NMAX)

      DO i=1,ndata
         CALL  eos_param(x(i),afunc,ma)
         tmp=1./sig(i)
         DO j=1,ma
            u(i,j)=afunc(j)*tmp
         ENDDO
         b(i)=y(i)*tmp
      ENDDO
      CALL svdcmp(u,ndata,ma,mp,np,w,v)
      wmax=0.
      DO j=1,ma
         IF (w(j) > wmax)wmax=w(j)
      ENDDO
      thresh=TOL*wmax
      DO j=1,ma
         IF (w(j) < thresh)  w(j)=0.
      ENDDO
      CALL svbksb(u,w,v,ndata,ma,mp,np,b,a)
      chisq=0.
      DO i=1,ndata
         CALL eos_param(x(i),afunc,ma)
         sum=0.
         DO j=1,ma
            sum=sum+a(j)*afunc(j)
         ENDDO
         chisq=chisq+((y(i)-sum)/sig(i))**2
      ENDDO

      END SUBROUTINE svdfit


      FUNCTION pythag(a,b)
      IMPLICIT NONE
      DOUBLE PRECISION ::  a, b, pythag 
      DOUBLE PRECISION :: absa, absb
      absa=ABS(a)
      absb=ABS(b)
      IF (absa > absb) THEN
         pythag=absa*SQRT(1.+(absb/absa)**2)
      ELSE
         IF(absb == 0.) THEN
            pythag=0.
         ELSE
            pythag=absb*SQRT(1.+(absa/absb)**2)
         ENDIF
      ENDIF

      END FUNCTION pythag

      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      IMPLICIT NONE
      INTEGER :: m,mp,n,np,NMAX
      DOUBLE PRECISION :: a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
      INTEGER :: i,its,j,jj,k,l,nm
      DOUBLE PRECISION :: anorm,c,f,g,h,s,scale,x,y,z, &
                          rv1(NMAX),pythag
      g=0.0
      scale=0.0
      anorm=0.0
      DO i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0.0
         s=0.0
         scale=0.0
         IF(i <= m)THEN
            DO k=i,m
               scale=scale+ABS(a(k,i))
            ENDDO
            IF(scale /= 0.0)THEN
               DO k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
               ENDDO
               f=a(i,i)
               g=-sign(SQRT(s),f)
               h=f*g-s
               a(i,i)=f-g
               DO j=l,n
                  s=0.0
                  DO k=i,m
                     s=s+a(k,i)*a(k,j)
                  ENDDO
                  f=s/h
                  DO k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
                  ENDDO
               ENDDO
               DO k=i,m
                  a(k,i)=scale*a(k,i)
               ENDDO
            ENDIF
         ENDIF
         w(i)=scale *g
         g=0.0
         s=0.0
         scale=0.0
         IF((i <= m).and.(i /= n))THEN
            DO k=l,n
               scale=scale+ABS(a(i,k))
            ENDDO
            IF(scale /= 0.0)THEN
               DO k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
               ENDDO
               f=a(i,l)
               g=-sign(SQRT(s),f)
               h=f*g-s
               a(i,l)=f-g
               DO k=l,n
                  rv1(k)=a(i,k)/h
               ENDDO
               DO j=l,m
                  s=0.0
                  DO k=l,n
                     s=s+a(j,k)*a(i,k)
                  ENDDO
                  DO k=l,n
                     a(j,k)=a(j,k)+s*rv1(k)
                  ENDDO
               ENDDO
               DO k=l,n
                  a(i,k)=scale*a(i,k)
               ENDDO
            ENDIF
         ENDIF
         anorm=max(anorm,(ABS(w(i))+ABS(rv1(i))))
      ENDDO
      DO i=n,1,-1
         IF(i < n)THEN
            IF(g /= 0.0)THEN
               DO j=l,n
                  v(j,i)=(a(i,j)/a(i,l))/g
               ENDDO
               DO j=l,n
                  s=0.0
                  DO k=l,n
                     s=s+a(i,k)*v(k,j)
                  ENDDO
                  DO k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
                  ENDDO
               ENDDO
            ENDIF
            DO j=l,n
               v(i,j)=0.0
               v(j,i)=0.0
            ENDDO
         ENDIF
         v(i,i)=1.0
         g=rv1(i)
         l=i
      ENDDO
      DO i=min(m,n),1,-1
         l=i+1
         g=w(i)
         DO j=l,n
            a(i,j)=0.0
         ENDDO
         IF(g /= 0.0)THEN
            g=1.0/g
            DO j=l,n
               s=0.0
               DO k=l,m
                  s=s+a(k,i)*a(k,j)
               ENDDO
               f=(s/a(i,i))*g
               DO k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
               ENDDO
            ENDDO
            DO j=i,m
               a(j,i)=a(j,i)*g
            ENDDO
         ELSE
            DO j= i,m
               a(j,i)=0.0
            ENDDO
         ENDIF
         a(i,i)=a(i,i)+1.0
      ENDDO
      DO k=n,1,-1
         DO its=1,30
            DO l=k,1,-1
               nm=l-1
               IF((ABS(rv1(l))+anorm) == anorm)  goto 2
               IF((ABS(w(nm))+anorm) == anorm)  goto 1
            ENDDO
1           c=0.0
            s=1.0
            DO i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               IF((ABS(f)+anorm) == anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=1.0/h
               c= (g*h)
               s=-(f*h)
               DO j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
               ENDDO
            ENDDO
2           z=w(k)
            IF(l == k)THEN
               IF(z < 0.0)THEN
                  w(k)=-z
                  DO j=1,n
                     v(j,k)=-v(j,k)
                  ENDDO
               ENDIF
               goto 3
            ENDIF
            IF(its == 30) THEN 
               WRITE(*,*) 'no convergence in svdcmp'; STOP
            ENDIF
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.D0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0
            s=1.0
            DO j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               DO jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
               ENDDO
               z=pythag(f,h)
               w(j)=z
               IF(z /= 0.0)THEN
                  z=1.0/z
                  c=f*z
                  s=h*z
               ENDIF
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               DO jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
               ENDDO
            ENDDO
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
         ENDDO
3        CONTINUE    
      ENDDO

      END SUBROUTINE svdcmp

