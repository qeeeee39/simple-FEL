      program test
        
        implicit none
        include 'mpif.h'
        
        integer :: ierr, rank, size
        character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
        integer :: name_len
        real*8 :: start_time, end_time, total_time
        
        integer, parameter :: pi = 3.14159d0, emass = 0.511d6, Npt = 100000
        real*8, parameter :: ini_gamma = 100d6/emass, phase=3.0d0
        real*8 :: beammat(Npt,6), gamma0, k, freq, len, accgrad
        
        !# mpi run
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_GET_PROCESSOR_NAME(processor_name, name_len, ierr)
        start_time = MPI_Wtime()
        
        !# main code
        
        call read_data(Npt, beammat, ierr)
        gamma0 = ini_gamma
        
        !0---------------
        call postproc(Npt, gamma0, beammat, 20)
        call linac(Npt, beammat, gamma0, 0.61685d0, 1.3d9, 20d0, 16d6, 5.5d0)
             
        !1---------------
        call postproc(Npt, gamma0, beammat, 21)
        call linac(Npt, beammat, gamma0, 0.61685d0, 3.9d9, 5d0, 10d6, phase)
        !2---------------
        call postproc(Npt, gamma0, beammat, 22)
        call compressor(Npt, beammat, gamma0, 0.05d0, 10d-2, 0.27416d0)
        !3---------------
        call postproc(Npt, gamma0, beammat, 23)
        call linac(Npt, beammat, gamma0, 0.61685d0, 1.3d6, 200d0, 16d6, 0d0/180d0*pi)
        !4---------------
        call postproc(Npt, gamma0, beammat, 24)
        call compressor(Npt, beammat, gamma0, -0.04d0, 20d-2, 0.27416d0)
        !5---------------
        call postproc(Npt, gamma0, beammat, 25)
        call linac(Npt, beammat, gamma0, 0.61685d0, 1.3d6, 400d0, 16d6, 0d0)
        !6---------------
        call postproc(Npt, gamma0, beammat, 26)
        
        !call scanphase(Npt, gamma0, beammat, phase)
        
        call felpara(gamma0, 200d-12, 9.0870689990160753d-005, 1d0, 3d-2, 30d0)
        !# end of mpi run
        end_time = MPI_Wtime()
        total_time = end_time - start_time
        if (rank == 0) then
          print *, 'Total time: ', total_time, ' seconds with ', size, ' processes'
        end if
        call MPI_FINALIZE(ierr)
        
      end program

!!! -----------------------------------------------------------------------      
      subroutine compressor(Npt, beammat, gamma0, R56, steplen, k)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(inout) :: beammat(Npt,6)
        real*8, intent(in) ::  gamma0, R56, steplen, k
        real*8 :: trfmat(6,6), a, thisptl(6), g
        
        !match the transfer matrix
        
        a = sqrt(k)        
        
        trfmat = 0d0
          
        trfmat(1,1) = cos(a*steplen)
        trfmat(1,2) = 1/a * (sin(a*steplen))
        trfmat(2,1) = -a*sin((a*a)*steplen)
        trfmat(2,2) = cos(a*steplen)
        
        trfmat(3,3) = cos(a*z)
        trfmat(3,4) = 1/a * (sin(a*steplen))
        trfmat(4,3) = -a*sin((a*a)*steplen)
        trfmat(4,4) = cos(a*steplen)
          
        trfmat(5,5) = 1
        trfmat(5,6) = 0
        trfmat(6,5) = 0
        trfmat(6,6) = 1        
        
        !compression, energy and spatial modulation
        
        T566 = -3d0/2d0 * R56
        U5666 = 2d0 * R56
        
        do i = 1, Npt
          thisptl = beammat(i,:)
          beammat(i,:) = matmul(trfmat,thisptl)
          g = beammat(i,6)/gamma0
          beammat(i,5) = beammat(i,5) + R56*g + T566*g*g + U5666*g*g*g
        end do
      
      
      end subroutine
!!! -----------------------------------------------------------------------    
      subroutine linac(Npt, beammat, gamma0, k, freq, len, accgrad, phase)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(inout) :: beammat(Npt,6), gamma0
        real*8, intent(in) :: k, freq, len, accgrad, phase
        real*8, parameter :: emass = 0.511d6
        
        halflen = len/2d0
        
        call linac_loc(Npt, beammat, gamma0, k, halflen, accgrad)
        gamma0 = gamma0 + (halflen * 1d0 * accgrad / emass * cos(phase))
        call linac_dg(Npt, beammat, freq, halflen, accgrad, phase)
        call linac_loc(Npt, beammat, gamma0, k, halflen, accgrad)
        gamma0 = gamma0 + (halflen * 1d0 * accgrad / emass * cos(phase))
        
      end subroutine
!!! ----------------------------------------------------------------------- 
      subroutine linac_dg(Npt, beammat, freq, steplen, accgrad, phase)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(in) :: freq, steplen, accgrad, phase
        real*8, intent(inout) :: beammat(Npt,6)
        real*8 :: wnum, edistr
        real*8, parameter :: c = 3.0d8, pi=3.14d0, emass = 0.511d6
        
        wnum = 2d0 * pi * freq / c
        
        do i = 1, Npt
          z_plus = beammat(i,5)
          edistr = cos(phase - wnum * z_plus) - cos(phase)
          beammat(i,6) = beammat(i,6) + (2d0 * steplen * 1d0 * accgrad / emass) * edistr
        enddo
      
      end subroutine
!!! ----------------------------------------------------------------------- 
      subroutine linac_loc(Npt, beammat, gamma0, k, steplen, accgrad)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(in) :: k, steplen, accgrad
        real*8, intent(inout) :: gamma0, beammat(Npt,6)
        real*8 :: trfmat(6,6), thisptl(6)
        real*8 :: gamma_square, beta_square, a
        real*8, parameter :: emass = 0.511d6
        
        gamma_square = gamma0 * gamma0
        beta_square = 1 - ( 1 / gamma_square )
        a = sqrt(k)
        
        trfmat = 0d0
        
        trfmat(1,1) = cos(a*steplen)
        trfmat(1,2) = 1/a * (sin(a*steplen))
        trfmat(2,1) = -a*sin((a*a)*steplen)
        trfmat(2,2) = cos(a*steplen)
        
        trfmat(3,3) = cos(a*steplen)
        trfmat(3,4) = 1/a * (sin(a*steplen))
        trfmat(4,3) = -a*sin(a*a*steplen)
        trfmat(4,4) = cos(a*steplen)
        
        trfmat(5,5) = 1
        trfmat(5,6) = (1/(gamma_square * beta_square)**(3d0/2d0)) * steplen
        trfmat(6,5) = 0
        trfmat(6,6) = 1
          
        do i=1,Npt 
          thisptl = beammat(i,:)
          beammat(i,:) = matmul(trfmat,thisptl)
        end do
        
      end subroutine
!!! -----------------------------------------------------------------------
      subroutine postproc(Npt, gamma0, mat, file_unit)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt, file_unit
        real*8, intent(in) :: mat(Npt,6), gamma0
        real*8, parameter :: emass = 0.511d6
        integer :: i
        real*8 :: bunchlength, energy, emitun, meanloc, sumdz, meane, sumde, mean_pp, mean_xx, mean_px, tlength
        real*8 :: x(Npt), px(Npt), y(Npt), py(Npt), z(Npt), dg(Npt), dz(Npt), de(Npt), pp(Npt), xx(Npt), xp(Npt)
        
        do i = 1, Npt
          x(i) = mat(i,1)
          px(i) = mat(i,2)
          y(i) = mat(i,3)
          py(i) = mat(i,4)
          z(i) = mat(i,5)
          dg(i) = mat(i,6)
        end do

        !# printout data to files
        do i = 1, Npt
          write(file_unit,110) x(i), px(i), y(i), py(i), z(i), dg(i)*0.511 
        enddo
        
        !# parameters calculation
        meanloc = sum(z) / Npt
        dz = z - meanloc
        do i = 1,Npt
          dz(i) = dz(i) * dz(i)
        enddo
        sumdz = sum(dz)
        rms_bunchlength = 2d0 * sqrt(sumdz/Npt)
        
        meane = sum(dg) / Npt
        de = dg - meane
        do i =1, Npt
          de(i) = de(i) * de(i)
        enddo
        sumde = sum(de)
        energyspread = 2d0 * sqrt(sumde/Npt)
        
        energy = meane + gamma0
        
        do i = 1, Npt
          pp(i) = px(i)*px(i)
          xx(i) = x(i)*x(i)
          xp(i) = px(i)*x(i)
        end do
        mean_pp = sum(pp) / Npt
        mean_xx = sum(xx) / Npt
        mean_px = sum(xp) / Npt
        emitun = sqrt(mean_pp * mean_xx - mean_px * mean_px)
        !!emitun = gamma0/2d0 * mean_xx * sqrt(sumdz/Npt) * sqrt(sumdz/Npt) /2d0**(1d0/2d0)
        
        
        !# printout the parameters to shell    
        print*, file_unit, ' ', 'output parameters ------------------'
        print*, 'rms bunch length (unit: m)', rms_bunchlength
        print*, 'beam energy (unit: gamma)', energy
        print*, 'rms energy spread (unit: gamma)', energyspread
        print*, 'emittance (unit: m)', emitun
        
110    format(6(1x,e15.7))    

      end subroutine
!!! ----------------------------------------------------------------------- 
      subroutine felpara(beamenergy, totalcharge, bunchlength, undK, lambdau, undlen)
        implicit real*8 (a-h,o-z)
        real*8, intent(in) :: beamenergy, totalcharge, bunchlength, undK, lambdau, undlen
        real*8 :: lambda, duration, pc, power, undPeriod, bandwidth
        real*8, parameter :: e = 1.602d-19, m = 9.109d-31, c = 2.997d8 
        
        lambda = lambdau * (1+ undK*undK /2d0) / (2d0 * beamenergy * beamenergy)
        
        duration = bunchlength / c
        pc = totalcharge / duration
        power = pc * beamenergy * m * c * c / e
        
        undPeriod = undlen / lambdau
        bandwidth = lambda / undPeriod
        
        print*, 'FEL estimation----------------------'
        print*, 'wavelength', lambda
        print*, 'power', power
        print*, 'bandwidth', bandwidth
        print*, '1/N_und = ', 1/undPeriod
            
      end subroutine felpara
!!! -----------------------------------------------------------------------
      subroutine scanphase(Npt, gamma0, mat, phase)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(in) :: mat(Npt,6), gamma0, phase
        integer, parameter :: numbins = 10
        integer :: i, hist(numbins)
        real*8 :: bunchlength, energy, emitun, meanloc, sumdz, meane, sumde, mean_pp, mean_xx, mean_px, tlength
        real*8 :: x(Npt), px(Npt), y(Npt), py(Npt), z(Npt), dg(Npt), dz(Npt), de(Npt), pp(Npt), xx(Npt), xp(Npt)
        
        do i = 1, Npt
          x(i) = mat(i,1)
          px(i) = mat(i,2)
          y(i) = mat(i,3)
          py(i) = mat(i,4)
          z(i) = mat(i,5)
          dg(i) = mat(i,6)
        end do
        
        !# parameters calculation
        meanloc = sum(z) / Npt
        dz = z - meanloc
        do i = 1,Npt
          dz(i) = dz(i) * dz(i)
        enddo
        sumdz = sum(dz)
        rms_bunchlength = 2d0 * sqrt(sumdz/Npt)
        
        meane = sum(dg) / Npt
        de = dg - meane
        do i =1, Npt
          de(i) = de(i) * de(i)
        enddo
        sumde = sum(de)
        energyspread = 2d0 * sqrt(sumde/Npt)
        
        energy = meane + gamma0
        
        do i = 1, Npt
          pp(i) = px(i)*px(i)
          xx(i) = x(i)*x(i)
          xp(i) = px(i)*x(i)
        end do
        mean_pp = sum(pp) / Npt
        mean_xx = sum(xx) / Npt
        mean_px = sum(xp) / Npt
        emitun = sqrt(mean_pp * mean_xx - mean_px * mean_px)
        
        !# print the result in the file
        open(11, file='scanphase', status='old', position='append')
        write(11,110) phase, rms_bunchlength, energy, energyspread, emitun
        close(11)
110     format(6(1x,e15.7))    

      end subroutine
!!! -----------------------------------------------------------------------
      subroutine read_data(Npt, beammat, ierr)
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: Npt
        real*8, intent(out) :: beammat(Npt,6)
        integer, intent(out) :: ierr
        integer :: i
        character(len=100) :: line
        real*8 :: x, px, y, py, z, dg
        ierr = 0
      
        !# Open the file fort.10 for reading
        open(unit=10, file='fort.10', status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
          print *, "Error opening the file"
          return
        endif
      
        !# Read the data line by line
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
      
        close(10)
          
      end subroutine read_data