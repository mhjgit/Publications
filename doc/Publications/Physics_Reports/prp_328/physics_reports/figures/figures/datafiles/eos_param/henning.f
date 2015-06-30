      program stupid
      implicit real*8(a-h,o-z)
      c = 0.2
      w0 = -16.
      DO i = 0,300
         x = i*0.01/0.16
         e = w0*((2.+c)*x/(1.+c*x)-x*x/(1.+c*x))
         deno = (1+c*x)
         top = 2.+c-2.*x         
         dedx = (w0*top/deno-e*c/deno)
         p = x*x*0.16*dedx
         x1 = x*2.*dedx
         x2 = x*x*w0*(2./deno+c*top/deno/deno)
         x3 = x*x*c*c*e/deno/deno
         x4 = x*x*c*dedx/deno
         dpdn=x1-x2+x3-x4

         WRITE(6,'(E10.4,2X,E10.4,2X,E10.4,2X,E10.4)') 
     &          x*0.16, e,p, dpdn/938.926
      ENDDO
      end
