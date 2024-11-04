       program initialdist
       
       implicit real*8 (a-h,o-z)
       integer, parameter :: Npt = 100000
       real*8, dimension(Npt) :: x,px,y,py,z,dg

       twopi = 4*asin(1.0)
       emass = 0.511e6

       !# energy distribution = a standard deviation (unit: gamma)
       sigdg = 2e3/emass
       
       !# rms size (unit: m)
       sig = 0.2e-3

       !# reference energy (unit: gamma)
       gamma0 = 100e6/emass

       !# geometric emittance (unit: m)
       emitun= 0.5e-6/gamma0
       
       !# mean transverse energy (unit: gamma)
       !!!sigxp = emitun/sig
       sigxp = (emitun/sig)**(2d0)
       
       print*, 'initial conditions ------------------------------'
       print*, 'reference energy', gamma0
       print*, 'energy spread', 2d0*sigdg
       print*, 'emittance', emitun
       print*, 'square root transverse energy', sigxp

       do i = 1, Npt
         call random_number(r0)
         z(i) = (3e-3 * r0) - 1.5e-3
         
         call gaussian(rg)
         dg(i) = sigdg * rg
         
         call gaussian(rg)
         x(i) = sig * rg
         
         call gaussian(rg)
         px(i) = sigxp * rg
         
         call gaussian(rg)
         y(i) = sig * rg
         
         call gaussian(rg)
         py(i) = sigxp * rg
         
         write(10,110)x(i),px(i),y(i),py(i),z(i),dg(i)
       enddo
110    format(6(1x,e15.7))

       end program

       subroutine gaussian(rg)
       real*8 :: r1,r2,rr,rg,phi
         twopi = 4*asin(1.0)
         call random_number(r1)
         call random_number(r2)
         phi = twopi*r1
         rr = sqrt(-2*log(r2))
         rg = rr*cos(phi)
       end subroutine


