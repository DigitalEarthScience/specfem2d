!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!=====================================================================

! for poro solver

  subroutine compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,myrank
  use constants, only: C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: ispec_is_poroelastic,nglob_poroelastic, &
                         NSOURCES,source_time_function,sourcearrays, &
                         islice_selected_source,ispec_selected_source, &
                         ibool,phistore,tortstore,rhoarraystore
  use specfem_par, only: NSTEP,  &
                         time_function_type, name_of_source_file, burst_band_width, f0_source,tshift_src, &
                         factor, t0, DT, SOURCE_IS_MOVING, &
                         time_stepping_scheme, stage_time_scheme, islice_selected_source, &
                         USE_TRICK_FOR_BETTER_PRESSURE, myrank, initialfield
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: accels_poroelastic,accelw_poroelastic

  integer :: it,i_stage
  double precision,external :: get_stf_poro

  !local variables
  integer :: i_source,i,j,iglob,ispec
  double precision :: phi,tort,rho_s,rho_f,rho_bar
  real(kind=CUSTOM_REAL) :: fac_s,fac_w
  !real(kind=CUSTOM_REAL) :: stf_used
  double precision :: stf_used
  double precision :: timeval,t_used

  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is poroelastic
      if (ispec_is_poroelastic(ispec)) then

        ! compute current time
        select case(time_stepping_scheme)
        case (1)
          ! Newmark
          timeval = dble(it-1)*DT
        case (2)
          ! LDDRK: Low-Dissipation and low-dispersion Runge-Kutta
          ! note: the LDDRK scheme updates displacement after the stiffness computations and
          !       after adding boundary/coupling/source terms.
          !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme.
          !       we therefore at an additional -DT to have the corresponding timing for the source.
          timeval = dble(it-1-1)*DT + dble(C_LDDRK(i_stage))*DT
        case (3)
          ! RK: Runge-Kutta
          ! note: similar like LDDRK above, displ(n+1) will be determined after stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          timeval = dble(it-1-1)*DT + dble(C_RK4(i_stage))*DT
        case (4)
          ! symplectic PEFRL
          ! note: similar like LDDRK above, displ(n+1) will be determined after final stage of stiffness/source/.. computations.
          !       thus, adding an additional -DT to have the same timing in seismogram as Newmark
          !
          !       for symplectic schemes, the current stage time step size is the sum of all previous and current coefficients
          !          sum( ALPHA_SYMPLECTIC(1:i_stage) ) * DT
          timeval = dble(it-1-1)*DT + dble(sum(ALPHA_SYMPLECTIC(1:i_stage))) * DT
        case default
          call exit_MPI(myrank,'Error invalid time stepping scheme chosen, please check...')
        end select

        t_used = timeval - t0 - tshift_src(i_source)

        ! source time function
        !stf_used = source_time_function(i_source,it,i_stage)
        stf_used = get_stf_poro(t_used,i_source)

        ! adds source contribution
        ! note: we use sourcearrays for both, collocated force and moment tensor forces
        !       (see setup in setup_souces_interpolation() routine), thus can write for both cases the same loop
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            ! poroelastic material
            phi = phistore(i,j,ispec)
            tort = tortstore(i,j,ispec)

            rho_s = rhoarraystore(1,i,j,ispec)
            rho_f = rhoarraystore(2,i,j,ispec)

            rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

            fac_s = real((1.d0 - phi/tort),kind=CUSTOM_REAL)
            fac_w = real((1.d0 - rho_f/rho_bar),kind=CUSTOM_REAL)

            ! solid contribution
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + &
                        fac_s * sourcearrays(1,i,j,i_source) * stf_used
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
                        fac_s * sourcearrays(2,i,j,i_source) * stf_used

            ! fluid contribution
            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + &
                        fac_w * sourcearrays(1,i,j,i_source) * stf_used
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + &
                        fac_w * sourcearrays(2,i,j,i_source) * stf_used
          enddo
        enddo
      endif
    endif ! if this processor core carries the source and the source element is elastic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_poro

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_poro_adjoint()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nrecloc,NSTEP,it, &
                         ispec_is_poroelastic,ispec_selected_rec_loc, &
                         ibool,initialfield,SIMULATION_TYPE, &
                         phistore,rhoarraystore, &
                         accels_poroelastic,accelw_poroelastic, &
                         source_adjoint,xir_store_loc,gammar_store_loc

  implicit none

  ! local parameters
  integer :: irec_local,i,j,iglob,ispec
  double precision :: phi,rho_s,rho_f,rho_bar

  ! checks if anything to do
  if (initialfield) return

  ! only for adjoint/kernel simulations
  if (.not. SIMULATION_TYPE == 3) return

  do irec_local = 1,nrecloc
    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_poroelastic(ispec)) then

      ! add source array
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          ! solid contribution
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + real(xir_store_loc(irec_local,i)*&
            gammar_store_loc(irec_local,j)*source_adjoint(irec_local,NSTEP-it+1,1),kind=CUSTOM_REAL)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + real(xir_store_loc(irec_local,i)*&
            gammar_store_loc(irec_local,j)*source_adjoint(irec_local,NSTEP-it+1,2),kind=CUSTOM_REAL)

          ! fluid
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - &
                real(rho_f/rho_bar * xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                source_adjoint(irec_local,NSTEP-it+1,1),kind=CUSTOM_REAL)
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - &
                real(rho_f/rho_bar * xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                source_adjoint(irec_local,NSTEP-it+1,2),kind=CUSTOM_REAL)
        enddo
      enddo

    endif ! if element is poroelastic

  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_poro_adjoint

!
!========================================================================
!

  double precision function get_stf_poro(t_used,i_source)

  ! prepares source_time_function array

  use constants, only: IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER,OUTPUT_FILES, &
                       SOURCE_DECAY_MIMIC_TRIANGLE, &
                       C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: NSTEP, NSOURCES, source_time_function, &
                         time_function_type, name_of_source_file, burst_band_width, f0_source,tshift_src, &
                         factor, t0, DT, SOURCE_IS_MOVING, &
                         time_stepping_scheme, stage_time_scheme, islice_selected_source, &
                         myrank, initialfield

  implicit none

  ! local parameters
  double precision :: stf_used, timeval, DecT, Tc, omegat, omega_coa,dummy_t,coeff, t_used, Nc
  double precision :: hdur,hdur_gauss

  integer :: it,i_source,ier,num_file
  integer :: i_stage

  character(len=150) :: error_msg1 = 'Error opening the file that contains the external source: '
  character(len=250) :: error_msg
  logical :: trick_ok

  ! external functions
  double precision, external :: comp_source_time_function_heaviside_hdur
  double precision, external :: comp_source_time_function_Gaussian,comp_source_time_function_dGaussian, &
    comp_source_time_function_d2Gaussian,comp_source_time_function_d3Gaussian,marmousi_ormsby_wavelet,cos_taper
  double precision, external :: comp_source_time_function_Ricker,comp_source_time_function_d2Ricker

  double precision :: stf
  
  ! determines source_time_function value for different source types
  select case (time_function_type(i_source))
  case (1)
      ! Ricker (second derivative of a Gaussian) source time function
      stf = - factor(i_source) * &
              comp_source_time_function_Ricker(t_used,f0_source(i_source))

  case (2)
      ! First derivative of a Gaussian source time function
      stf = - factor(i_source) * &
                comp_source_time_function_dGaussian(t_used,f0_source(i_source))

  case (3,4)
      ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
      stf = - factor(i_source) * &
                  comp_source_time_function_Gaussian(t_used,f0_source(i_source))

  case (5)
    ! Heaviside source time function (we use a very thin error function instead)
    hdur = 1.d0 / f0_source(i_source)
    hdur_gauss = hdur * 5.d0 / 3.d0

    ! convert the half duration for triangle STF to the one for Gaussian STF
    hdur_gauss = hdur_gauss / SOURCE_DECAY_MIMIC_TRIANGLE

    ! quasi-Heaviside
    stf = - factor(i_source) * &
                            comp_source_time_function_heaviside_hdur(t_used,hdur_gauss)

  case (6)
    ! ocean poros type I
    DecT = t0 + tshift_src(i_source)
    Tc = 4.d0 / f0_source(i_source) + DecT
    omega_coa = TWO * PI * f0_source(i_source)

    if (timeval > DecT .and. timeval < Tc) then
      ! source time function from Computational Ocean Acoustics
      omegat =  omega_coa * ( timeval - DecT )
      stf = factor(i_source) * HALF * &
            sin( omegat ) * ( ONE - cos( QUARTER * omegat ) )
      !stf = factor(i_source) * HALF / omega_coa / omega_coa * &
      !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )
    else
      stf = ZERO
    endif

  case (7)
    ! ocean poros type II
    DecT = t0 + tshift_src(i_source)
    Tc = 4.d0 / f0_source(i_source) + DecT
    omega_coa = TWO*PI*f0_source(i_source)
      !Tc = 1.d0 / f0_source(i_source) + DecT ! For source 1 OASES
      !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
      !  stf = factor(i_source) * ( &  ! Source 1 OASES
      !            0.75d0 - cos(omega_coa*t_used) + 0.25d0*cos(TWO*omega_coa*t_used))
      !else
      !  stf = ZERO
      !endif
      if (timeval > DecT .and. timeval < Tc) then
        ! source time function from Computational Ocean Acoustics
        omegat =  omega_coa * ( timeval - DecT )
        !stf = factor(i_source) * HALF / omega_coa / omega_coa * &
        !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
        !     8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( timeval - DecT ) + 1./15.*4./f0_source(i_source))
        stf = factor(i_source) * HALF / omega_coa / omega_coa * &
               ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
                8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )
      else if (timeval > DecT) then
        stf = &
               - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0_source(i_source))
      else
        stf = ZERO
      endif
!      if (timeval > DecT .and. timeval < Tc) then
!        ! source time function from Computational Ocean Acoustics
!        omegat =  omega_coa * ( timeval - DecT )
!        !stf = factor(i_source) * HALF / omega_coa / omega_coa * &
!        !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - &
!        !     8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) -1./15.*( timeval - DecT ) + 1./15.*4./f0_source(i_source))
!        stf = factor(i_source) * HALF / omega_coa / omega_coa * &
!               ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) + &
!                8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )
!      else if (timeval > DecT) then
!        stf = &
!               - factor(i_source) * HALF / omega_coa / 15.d0 * (4.d0 / f0_source(i_source))
!      else
!        stf = ZERO
!      endif

  case (8)
    ! external type
    ! opens external file to read in source time function
    if (it == 1) then
      ! reads in from external source time function file
      open(unit=num_file,file=trim(name_of_source_file(i_source)),status='old',action='read',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening source time function file: ',trim(name_of_source_file(i_source))
        error_msg = trim(error_msg1)//trim(name_of_source_file(i_source))
        call exit_MPI(myrank,error_msg)
      endif
    endif

    ! reads in 2-column file values (time value in first column will be ignored)
    ! format: #time #stf-value
    read(num_file,*,iostat=ier) dummy_t, stf
    if (ier /= 0) then
      print *,'Error reading source time function file: ',trim(name_of_source_file(i_source)),' at line ',it
      print *,'Please make sure the file contains the same number of lines as the number of timesteps NSTEP ',NSTEP
      call exit_MPI(myrank,'Error reading source time function file')
    endif

    ! closes external file
    if (it == NSTEP) close(num_file)

    ! amplifies STF by factor
    ! note: the amplification factor will amplify the external source time function.
    !       in case this is not desired, one just needs to set the amplification factor to 1 in DATA/SOURCE:
    !         factor  = 1.0
    coeff = factor(i_source)
    stf = stf * coeff

  case (9)
    ! burst type
    DecT = t0 + tshift_src(i_source)
    t_used = (timeval-t0-tshift_src(i_source))
    Nc = TWO * f0_source(i_source) / burst_band_width(i_source)
    Tc = Nc / f0_source(i_source) + DecT
    omega_coa = TWO*PI*f0_source(i_source)

      if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
        stf = - factor(i_source) * &
                  0.5d0*(ONE-cos(omega_coa*t_used/Nc))*sin(omega_coa*t_used)
      !else if (timeval > DecT) then
      !  stf = ZERO
      else
        stf = ZERO
      endif

  case (10)
    ! Sinus source time function
    omega_coa = TWO*PI*f0_source(i_source)
      ! First derivative of a Gaussian source time function
      stf = factor(i_source) * sin(omega_coa*t_used)

  case (11)
      ! Marmousi_ormsby_wavelet
      hdur = 1.0 / 35.0
      stf = factor(i_source) * &
                cos_taper(t_used,hdur) * marmousi_ormsby_wavelet(PI*t_used)

  case default
    call exit_MPI(myrank,'unknown source time function')

  end select

  ! return value
  get_stf_poro = stf

  end function get_stf_poro
