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
         INTEGER, PUBLIC, PARAMETER :: number_central_densities=30
!     the differential equations of Hartle, excluding the quadropole terms
         INTEGER, PUBLIC , PARAMETER :: number_differential_eqs= 7
!     number of runge-kutta iterations
         INTEGER, PUBLIC , PARAMETER :: max_rk4_steps=100000
!     Dimensionless step in RK4 procedure
         DOUBLE PRECISION, PUBLIC, PARAMETER :: diff_eq_step = 0.0001
!     Density from which one wants the central density to be calculated 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: first_density = 0.1
!     Last central density 
         DOUBLE PRECISION, PUBLIC, PARAMETER :: last_density = 1.8
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
      INTEGER :: j, i, l
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: cen_dens, &
                        radius, smass, inertia, omega
      DOUBLE PRECISION, DIMENSION (number_differential_eqs) :: ynew, y, ders
      DOUBLE PRECISION ::  star_radius, ync, density_step, ang_mom, e_rho,xx


!     Read in equation of state and make fit of it
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
!      WRITE( 6,*) 'central density, radius, mass, Omega, I'
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
         y(6)=fourpi/15.*const_1*(1.+y(1))*(2.+1/dpdn(ync))* &
              (y(5)**2)*(star_radius**5)
!     Initial value of monople function p_0
         y(7)=1./3.*(y(5)*star_radius)**2
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
            e_rho=energy_density(y,xx)
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
         WRITE(6,'(E12.3,7E12.4)') cen_dens(i),radius(i),smass(i), ang_mom,omega(i),inertia(i), y(6)+(const_1**2)*const_2*  & 
         (ang_mom**2)/(star_radius**3)+y(2)
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
         ders(4)=(ders(3)*0.5-const_2*(y(2)-r*ders(2))/ &
                 (r**2-2.*const_2*r*y(2)))*(y(4)+4.*(r**3)*y(5))
!     derivative of angular velocity given by Omega-omega
         ders(5)=y(4)/r**4
!     derivative of the monopole function m_0
         IF ( dpdn(xx) > 0.) THEN
         ders(6)=const_1*r*r*((e_rho+y(1))**2)*y(7)/dpdn(xx) + &
                  (r**4)*(ders(5)**2)*EXP(-y(3))*(1.-2*const_2*y(2)/r)/12./const_2 + &
                 2/3.*const_1*(r**4)*EXP(-y(3))*(e_rho+y(1))*(y(5)**2)
         ENDIF
!     derivative of the monopole function p_0
         ders(7)=1./3.*r*EXP(-y(3))*(y(5)**2)*(2*y(5)-r*y(5)*ders(3)+2*r*ders(5)) + &
                 (r**3)*EXP(-y(3))*(ders(5)**2)/12. - &
                 const_1*const_2*r*(e_rho+y(1))*y(7)/(1.-2*const_2*y(2)/r) - &
                 y(6)*const_2*(1.+2*const_1*const_2*y(1)*r*r)/(r*r)/  &
                 ((1.-2*const_2*y(2)/r)**2)
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


