      program test
        !implicit none
        implicit real*8 (a-h,o-z)
        integer :: ierr
        integer, parameter :: Npt = 1000
        real*8, parameter :: emass = 0.511d6, gamma0 = 100d6/emass
        real*8 :: beammat(Npt,6), k, freq, len, accgrad
        real*8, parameter :: realnumptl = 200d-12/1.602d-19, q = realnumptl/Npt, m = 9.109d-31, c = 2.997d8
        !!realnumptl is the real number particle
        !!q is the charge of each macroparticle and also the number of e charge per macroparticle, unit is [e], 
        
        !use beam
        
        
        !---------above shows the usage of module in other file
        
        !call data_passing('me')
        
        
        !!! input initial conditions:
        !current = 100e-12
        !energy = 100e6
        !bunchlength = 3e-3
        !energystandeviate = 2e3
        !energyspread = 0
        !rmssize = 0.2e-3
        
        !call initial_distribution(Npt,beammat) !generate distributions and transfer matrix
        
        call read_data(Npt, beammat, ierr)
        !print*, '***', beammat(3,3)
        
        !call linac(Npt, beammat, q, gamma0, 0.61685d0, 1.3d9, 20d0, 16d6)
        !call linac(Npt, beammat, q, gamma0, 0.61685d0, 3.9d9, 5d0, 10d6)
        !call compressor(Npt, beammat, 5d-2, 10d-2, 0.27416d0)   !BC1
        !call linac(Npt, beammat, q, gamma0, 0.61685d0, 1.3d6, 200d0, 16d6)
        !call compressor(Npt, beammat, 5d-2, 20d-2, 0.27416d0)
        !call linac(Npt, beammat, q, gamma0, 0.61685d0, 1.3d6, 400d0, 16d6)
        
        call plothistogram(Npt, q, gamma0,beammat) !plot the figures and print out the parameters
        
        !call felpara()
        
        
      end program

!!! -----------------------------------------------------------------------
      subroutine data_passing(name)
       character(len=*), intent(in) :: name
       print*, 'hi',' ', name  
      end subroutine
!!! ------------------------------------------------------------------------
      subroutine read_data(Npt, beammat, ierr)
          implicit none
          integer, intent(in) :: Npt
          real*8, intent(out) :: beammat
          integer, intent(out) :: ierr
          integer :: i
          character(len=100) :: line
          real*8 :: x, px, y, py, z, dg
          ierr = 0
      
          ! Open the file fort.10 for reading
          open(unit=10, file='fort.10', status='old', action='read', iostat=ierr)
          if (ierr /= 0) then
            print *, "Error opening the file"
            return
          endif
      
          ! Read the data line by line
          do i = 1, Npt
            read(10,'(A)', iostat=ierr) line
            if (ierr /= 0) then
              print *, "Error reading the file or end of file reached"
              exit
            endif
            read(line, *) x, px, y, py, z, dg
            beammat(i,1) = x
            beammat(i,2) = px
            beammat(i,3) = y
            beammat(i,4) = py
            beammat(i,5) = z
            beammat(i,6) = dg
          end do
      
          ! Close the file
          close(10)
      end subroutine read_data
!!! -----------------------------------------------------------------------    
      subroutine linac(Npt, beammat, q, gamma0, k, freq, len, accgrad)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(inout) :: beammat(Npt,6)
        real*8, intent(in) :: q, gamma0, k, freq, len, accgrad
        !real*8, parameter :: m = 9.109d-31, c = 2.997d8
        !for macroparticle
        real*8, parameter :: mc2 = 0.511d6 * q
        
        
        !need to seperate the length into two half
        !the mutval is done by linactrfmat subroutine
        !output the beammat
        
        halflen = len/2d0
        
        call linactrfmat(Npt, beammat, gamma0, k, freq, halflen, accgrad)
        beammat(:,6) = beammat(:,6) + (halflen * q * accgrad) / mc2
        call linactrfmat(Npt, beammat, gamma0, k, freq, halflen, accgrad)
        beammat(:,6) = beammat(:,6) + (halflen * q * accgrad) / mc2
      
      end subroutine
!!! -----------------------------------------------------------------------          !# only do the matrix part, the energy gain is done in the main linac subroutine 
      subroutine linactrfmat(Npt, beammat, gamma0, k, freq, z, accgrad)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(in) :: gamma0, k, freq, z, accgrad
        real*8, intent(inout) :: beammat(Npt,6)
        real*8 :: trfmat(6,6), gamma, beta, thisptl(6), pp
        real*8, parameter :: m = 9.109d-31, c = 2.997d8
        
        do i = 1,Npt
        
          !construction of the linac matrix
          
          !\\\\\ '/q'
          gamma_square = (beammat(i,6)/q + gamma0) * (beammat(i,6)/q + gamma0)
          beta_square = 1 - ( 1 / gamma_square )
          a = sqrt(k)
          
          trfmat = 0d0
          
          trfmat(1,1) = cos(a*z)
          trfmat(1,2) = 1/a * (sin(a*z))
          trfmat(2,1) = -a*sin((a*a)*z)
          trfmat(2,2) = cos(a*z)
          
          trfmat(3,3) = cos(a*z)
          trfmat(3,4) = 1/a * (sin(a*z))
          trfmat(4,3) = -a*sin(a*a*z)
          trfmat(4,4) = cos(a*z)
          
          trfmat(5,5) = 1
          trfmat(5,6) = (1/(gamma_square * beta_square) ) * z
          trfmat(6,5) = 0
          trfmat(6,6) = 1
          
          !*** see the content of matrix---------------------
          !print*, '***trfmat',trfmat(5,6)
          !do i = 1, 6
          !    do j = 1, 6
          !        print *, trfmat(i, j)
          !    end do
          !    print *
          !end do
          !print*, shape(trfmat)
          !----------------------------------------------
          
          !perform the multval
          
          thisptl = beammat(i,:)
          !print*, shape(thisptl)
          beammat(i,:) = matmul(trfmat,thisptl)
          !print*, beammat(i,5)
        end do
        
        
      end subroutine
!!! -----------------------------------------------------------------------      
      subroutine compressor(Npt, beammat, R56, z, k)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(inout) :: beammat(Npt,6)
        real*8, intent(in) ::  R56, z, k
        real*8 :: trfmat(6,6), a, thisptl(6)
        
        !match the transfer matrix
        
        a = sqrt(k)        
        
        trfmat = 0d0
          
        trfmat(1,1) = cos(a*z)
        trfmat(1,2) = 1/a * (sin(a*z))
        trfmat(2,1) = -a*sin((a*a)*z)
        trfmat(2,2) = cos(a*z)
          
        trfmat(3,3) = cos(a*z)
        trfmat(3,4) = 1/a * (sin(a*z))
        trfmat(4,3) = -a*sin(a*a*z)
        trfmat(4,4) = cos(a*z)
          
        trfmat(5,5) = 1
        trfmat(5,6) = R56
        trfmat(6,5) = 0
        trfmat(6,6) = 1        
      
        do i = 1, Npt
          thisptl = beammat(i,:)
          !print*, shape(thisptl)
          beammat(i,:) = matmul(trfmat,thisptl)
        end do
      
      
      end subroutine
!!! -----------------------------------------------------------------------
!!! for transferring the read-in bunch data into matrix 

      subroutine initial_distribution(Npt,beammat)
       implicit real*8 (a-h,o-z) 
       integer, intent(in) :: Npt
       real*8, dimension(Npt) :: x,px,y,py,z,dg
       real*8, intent(out) :: beammat(Npt,6)
       
       twopi = 4*asin(1.0)
       emass = 0.511d6
  
       !energy spread distribution in a standard deviation (unit: gamma)
       sigdg = 2e3/emass
         
       !rms size
       sig = 0.2e-3
  
       !reference energy (unit: gamma)
       gamma0 = 100d6/emass
  
       !emittance 
       emitun= 0.5e-6/gamma0
  
       !square root of mean transverse energy (emitun/sig = \sqrt(MTE)/mc^2) (unit: gamma)
       sigxp = emitun/sig
       
       print*, 'initial conditions ------------------------------'
       print*, 'bunch length', 3e-3
       print*, 'reference energy', gamma0
       print*, 'emittance', emitun
       print*, 'EMS size', sig
       print*, 'square root transverse energy', sigxp

       
       do i = 1, Npt
         call random_number(r0)
         z(i) = 3e-3*r0-1.5e-3
         
         call gaussian(rg)
         dg(i) = sigdg*rg
         call gaussian(rg)
         x(i) = sig*rg
         call gaussian(rg)
         px(i) = sigxp*rg
         call gaussian(rg)
         y(i) = sig*rg
         call gaussian(rg)
         py(i) = sigxp*rg
         write(10,110)x(i),px(i),y(i),py(i),z(i),dg(i)
         
       enddo
110    format(6(1x,e15.7))

       do i = 1, Npt
         beammat(i,1) = x(i)
         beammat(i,2) = px(i)
         beammat(i,3) = y(i)
         beammat(i,4) = py(i)
         beammat(i,5) = z(i)
         beammat(i,6) = dg(i)
       enddo
      !print*, '***', beammat(3,3)
      
      
      end subroutine initial_distribution
      
            
!!! -----------------------------------------------------------------------
      subroutine gaussian(rg)
       real*8 :: r1,r2,rr,rg,phi
         twopi = 4*asin(1.0)
         call random_number(r1)
         call random_number(r2)
         phi = twopi*r1
         rr = sqrt(-2*log(r2))
         rg = rr*cos(phi)

      end subroutine
      
      
!!! -----------------------------------------------------------------------    
      subroutine plothistogram(Npt, q, gamma0,mat)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(in) :: mat(Npt,6), q, gamma0
        integer, parameter :: numbins = 10
        integer :: i, hist(numbins)
        real*8 :: bunchlength, energy, emitun, meanloc, sumdz, mean_pp, mean_xx, mean_px, tlength
        real*8 :: x(Npt), px(Npt), y(Npt), py(Npt), z(Npt), dg(Npt), dz(Npt), pp(Npt), xx(Npt), xp(Npt)
        
        !print*, '***', mat(3,3)
        
        !# define the bin number
        
        !# define the value range
        do i = 1, Npt
          x(i) = mat(i,1)
          px(i) = mat(i,2)
          y(i) = mat(i,3)
          py(i) = mat(i,4)
          z(i) = mat(i,5)
          dg(i) = mat(i,6)
        end do
        
        !# printout the parameters
        
        meanloc = sum(z) / Npt
        dz = z - meanloc
        sumdz = sum(dz)
        !print*, 'meanloc', meanloc
        !print*, 'sumdz', sumdz
        rms_bunchlength = 2d0 * sqrt(sumdz/Npt)
        
          
        energy = maxval(dg) / q + gamma0
        
        do i = 1, Npt
          pp(i) = px(i)*px(i)
          xx(i) = x(i)*x(i)
          xp(i) = px(i)*x(i)
        end do
        mean_pp = sum(pp) / Npt
        mean_xx = sum(xx) / Npt
        mean_px = sum(xp) / Npt
        emitun = sqrt(mean_pp * mean_xx - mean_px * mean_px)
        
        print*, 'output parameters ------------------------------'
        print*, 'rms bunch length', rms_bunchlength
        print*, 'peak energy', energy
        print*, 'emittance', emitun
        
        !!!plot histogram and printout data
        
        !# classificate the particles
        !## bunch length histogram
        !tlength = maxval(z) - minval(z)
        !binwidth = tlength / numbins
        !minz = minval(z)
        !do i = 1, Npt
        !  index = int( (z(i)-minz) / tlength ) + 1
        !  if (index > 0 .and. index <= numbins) then
        !    hist(index) = hist(index) + 1
        !  end if
        !enddo
        
        !# output the file for displaying histograms
        !## bunch length
        !do i = 1, numbins
        !  write(11,*) (minz + binwidth * i), hist(i) 
        !enddo

        do i = 1, Npt
          write(11,*) z(i) 
        enddo
        
        
        !## energy distribution
        do i = 1, Npt
          write(12,110) z(i), (dg(i)/q + gamma0)
        enddo
        
        !## phase space emittance
        do i = 1, Npt
          write(13,110) x(i), px(i), y(i), py(i)
        enddo
        
110    format(6(1x,e15.7))    

      end subroutine

!!! ------------------------------------------------------------------------
