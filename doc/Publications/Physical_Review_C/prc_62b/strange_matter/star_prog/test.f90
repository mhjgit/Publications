!            Program to calculate the mass (SM), 
!            metric function, angular velocity, moment of inertia and 
!            radius (r) for a slowly rotating neutron star, with a  
!            user supplied EOS. 
!            For more details about the equations used see Hartle,
!            ApJ 150 (1967) p. 1005
!            Only monopole terms included here.
!            This is the Fortran 90 version.
!            Free format is employed.
!
!
!            Coded by : M. Hjorth-Jensen            
!                       Dept. of Physics, Univ. of  
!                       Oslo, Norway, june/1991     
!                       based on white dwarf prog   
!                       whited.c of march 91
!            Language : Fortran 90                              
!            Last Fortran 90 update : 22/june/2001
!
!            You should have the energy per particle as function of 
!            particle density n, where n has units fm-3.
!            Note that the EOS you supply MUST be in units of MeVfm^-3.
!            Although all quantities in the program are dimensionless
!            when solving the differential equations, observe that
!            the units needed for the EoS and other observables are
!            [rho] = MeVfm-3  energy per baryon 
!            user,
!            [n]   = 1/fm^3    particle density                     
!            [P]   = MeV/fm^3  pressure                       
!
!            Recall also that the constants defined 
!            below have the following units
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
!            r=r/r0 with r0=1.D+19 fm= 10km
!            m=m/m_solar

!
!          Observables given by print_output:
!          Prints out central density in fm^-3, radius in km,
!          neutron star mass in terms of solar masses, angular velocity
!          at the surface (in units of s^-1 10km) and Moment of inertia
!          (in units of solar masses * (10km)^2)

!
!     Global constants used to render equations dim-less
!

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
!     constant used to render equations free of constants           
         DOUBLE PRECISION, PUBLIC, PARAMETER :: c1=fourpi*(r0**2)*g
!     constant needed to get mass in terms of solar masses, to be
!     inserted at the end of the calc
         DOUBLE PRECISION, PUBLIC, PARAMETER :: c2=r0/g/sm_solar
!     number of central densites chosen
         INTEGER, PUBLIC, PARAMETER :: number_central_densities=30
!     the differential equations of Hartle, excluding the quadropole terms
         INTEGER, PUBLIC , PARAMETER :: number_differential_eqs= 7
!     number of runge-kutta iterations
         INTEGER, PUBLIC , PARAMETER :: max_rk4_steps=100000
!     Dimensionless step in RK4 procedure
         DOUBLE PRECISION, PUBLIC, PARAMETER :: diff_eq_step = 0.0001
!     Density from which one wants the central density to be calculated 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: first_density = 0.2
!     Last central density 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: last_density = 1.6
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
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: polynom_terms1
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: polynom_terms2
         INTEGER :: number_terms1 , number_terms2
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

            IF ( x <= 0.25 ) THEN
               rho=polynom_terms1(1)
               DO i=2,number_terms1
                  rho=rho+polynom_terms1(i)*(x**(i-1))
               ENDDO
            ELSE
               rho=polynom_terms2(1)
               DO i=2,number_terms2
                  rho=rho+polynom_terms2(i)*(x**(i-1))
               ENDDO
            ENDIF
            END FUNCTION rho
!
!     derivative of rho
!
            FUNCTION drho(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: drho
            DOUBLE PRECISION,  INTENT(IN) :: x
            INTEGER :: i
            IF ( x <= 0.25 ) THEN
               drho=polynom_terms1(2)
               DO i=3,number_terms1
                  drho=drho+FLOAT(i-1)*polynom_terms1(i)*(x**(FLOAT(i-2))) 
               ENDDO
            ELSE
               drho=polynom_terms2(2)
               DO i=3,number_terms2
                  drho=drho+FLOAT(i-1)*polynom_terms2(i)*(x**(FLOAT(i-2))) 
               ENDDO
            ENDIF
            END FUNCTION drho
!
!     second derivative of rho
!
            FUNCTION ddrho(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: ddrho
            DOUBLE PRECISION,  INTENT(IN) :: x
            INTEGER :: i
            IF ( x <= 0.25 ) THEN
               ddrho=2.*polynom_terms1(3)
               DO i=4,number_terms1
                  ddrho=ddrho+FLOAT((i-1)*(i-2))*polynom_terms1(i)*&
                           (x**(FLOAT(i-3))) 
               ENDDO
            ELSE
               ddrho=2.*polynom_terms2(3)
               DO i=4,number_terms2
                  ddrho=ddrho+FLOAT((i-1)*(i-2))*polynom_terms2(i)*&
                           (x**(FLOAT(i-3))) 
               ENDDO
            ENDIF

            END FUNCTION ddrho


!
!     pressure in units of MeV/fm^3
!     pressure = density*d rho/d density - rho
! 
            FUNCTION press(x)
            IMPLICIT NONE
            DOUBLE PRECISION :: press
            DOUBLE PRECISION, INTENT(IN)  :: x
            INTEGER :: i
!  EoS ill-defined for such small densities, yields easily loss of
!  numerical precision

            IF ( x < 0.001) THEN
               press=0.
            ELSE
               press=x*drho(x)-rho(x)
            ENDIF
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
            IF ( x < 0.001) THEN
               dpdn=0.
            ELSE
               dpdn=ddrho(x)*x
            ENDIF

            END FUNCTION dpdn


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
      INTEGER :: j, i, l
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: cen_dens, &
                        radius, smass, inertia, omega
      DOUBLE PRECISION, DIMENSION (number_differential_eqs) :: ynew, y, ders
      DOUBLE PRECISION ::  star_radius, ync, density_step, ang_mom, e_rho, rho_c, xx


!     Read in equation of state, in polynomial form
      READ(5,*)
      READ(5,*) number_terms1 , number_terms2
      ALLOCATE ( polynom_terms1(number_terms1))
      ALLOCATE ( polynom_terms2(number_terms2))
      READ(5,*) 
      DO i =1, number_terms1
         READ(5,*) polynom_terms1(i)
      ENDDO
      READ(5,*)
      DO i =1, number_terms2
         READ(5,*) polynom_terms2(i)
      ENDDO
!     steps in choice of central densities
      density_step=(last_density-first_density)/number_central_densities   

!     Fixed central energy density in units MeV/fm^3, in order to render
!     all equations free of constants and al quantities dim-less        

      central_energy_density=1./c1         

!     loop over central densities
      DO i=1, number_central_densities                         
!     Preparation of starting values
!     central density in fm^-3
         ync=first_density+density_step*i   
!     central pressure in units of MeV/fm^-3
         pressure=press(ync)
!     Dimensionless pressure and energy at the centre (function of central dens)
         y(1)=pressure/central_energy_density;  rho_c=rho(ync)/central_energy_density  
!     Dimensionless start radius
         star_radius=diff_eq_step/10.
!     Dimensionless mass of the star at the centre 
         y(2)=(star_radius**3)/3.              
!     Metric function at the centre of the star
         y(3)=0.
!     Start values of the function w (y(5)) and u (y(4)) used in the calculation of
!     the angular velocity Omega and the derivative dudr of u
         y(4)=0.
         y(5)=0.1
!     Initial value of monopole function m_0
         y(6)=(1./15.)*(rho_c+y(1))*(2.+1/dpdn(ync))* &
              (y(5)*y(5))*(star_radius**5)
!     Initial value of monopole function p_0
         y(7)=0.3333333*(y(5)*star_radius)*(y(5)*star_radius)
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
            y(1)=y(1)+ders
            e_rho=energy_density(y,xx)
!     new pressure, now in units of MeV/fm^-3
            pressure=y(1)*central_energy_density
         ENDDO

!     Preparing output 
!     radius in units of 10 km
         radius(i)=star_radius*10.
!     mass in units of solar masses                      
         smass(i)=c2*y(2)                        
!     angular velocity in units of s^-1 * 10 km
         omega(i)=(y(5)+y(4)/3./(star_radius**3))
!     Calculate the angular momentum J
         ang_mom=(star_radius**4)*ders(5)/6.
!     central density in units of fm^-3
!        cen_dens(i)=LOG10(ync/xkg*1.D+42*938.926)
         cen_dens(i)=ync                       
!     moment of inertia in units of solar mass * (10km)^2, factor 100
         inertia(i)=(ang_mom/omega(i))*100.
         WRITE(6,'(E12.3,7E12.4)') cen_dens(i),radius(i), smass(i), &
         ang_mom,omega(i), inertia(i), (y(6)+(ang_mom**2)/(star_radius**3)+y(2))*c2
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
      DOUBLE PRECISION :: e_rho, xx, tov_value
      DOUBLE PRECISION, INTENT(IN)  :: r
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(OUT)  :: ders
      DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN)  :: y

      IF(y(1) > 0.) THEN
         e_rho=energy_density(y,xx)
!     TOV equation=dp/dr
         tov_value=(y(2)+(r*r*r*y(1)))/(r*r-2.*r*y(2))
         ders(1)=-(e_rho+y(1))*tov_value
!     dm/dr=r*r*energydensity, gravitational mass as function of r
         ders(2)=r*r*e_rho
!     d\nu/dr=-dp/dr*2./(e_rho+y(1)), metric function
         ders(3)=2.*tov_value
!     derivative of the function u used to evaluate omega=y(5)
         ders(4)=(ders(3)*0.5-(y(2)-r*ders(2))/ &
                 (r*r-2.*r*y(2)))*(y(4)+4.*r*r*r*y(5))
!     derivative of angular velocity given by Omega-omega
         ders(5)=y(4)/r**4
!     derivative of the monopole function m_0
         IF ( dpdn(xx) > 0.) THEN
         ders(6)=r*r*((e_rho+y(1))*y(1))*y(7)/dpdn(xx) + &
                  (r*r*r*r)*(ders(5)**2)*EXP(-y(3))*(1.-2*y(2)/r)/12.+ &
                  0.6666667*(r*r*r*r)*EXP(-y(3))*(e_rho+y(1))*(y(5)*y(5))
         ENDIF
!     derivative of the monopole function p_0
         ders(7)=1./3.*r*EXP(-y(3))*(y(5)**2)*(2*y(5)-r*y(5)*ders(3)+2*r*ders(5)) + &
                 (r**3)*EXP(-y(3))*(ders(5)**2)/12. - &
                 (e_rho+y(1))*y(7)/(1.-2.*y(2)/r) - &
                 y(6)*(1.+2.*y(1)*r*r)/(r*r)/  &
                 ((1.-2.*y(2)/r)**2)
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


