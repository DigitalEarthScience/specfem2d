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

! for acoustic solver

  subroutine compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,myrank
  use constants, only: C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic, &
                         NSOURCES,source_type,sourcearrays, &
                         islice_selected_source,ispec_selected_source,ibool,kappastore
  use specfem_par, only: NSOURCES, &
                         tshift_src, &
                         t0, DT,&
                         time_stepping_scheme,islice_selected_source, &
                         myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic
  integer,intent(in) :: it,i_stage
  double precision,external :: get_stf_acoustic
  !local variables
  integer :: i_source,i,j,iglob,ispec
  !real(kind=CUSTOM_REAL) :: stf_used
  double precision :: stf_used
  double precision :: timeval,t_used

  do i_source = 1,NSOURCES
    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is acoustic
      if (ispec_is_acoustic(ispec)) then

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
        stf_used = get_stf_acoustic(t_used,i_source)

        ! collocated force
        ! beware, for an acoustic medium, the source is pressure divided by Kappa of the fluid
        if (source_type(i_source) == 1) then
          ! forward wavefield
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                                  real(sourcearrays(1,i,j,i_source) * stf_used / kappastore(i,j,ispec),kind=CUSTOM_REAL)
            enddo
          enddo

        else if (source_type(i_source) == 2) then
          ! moment tensor
          call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
        endif

      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic

!
!=====================================================================
!

  subroutine compute_add_sources_acoustic_moving_sources(potential_dot_dot_acoustic,it,i_stage)

! This subroutine is the same than the previous one but with a moving source

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,NGLJ,TINYVAL,IMAIN

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic, &
                         NSOURCES,source_type,source_time_function, &
                         islice_selected_source,ispec_selected_source, &
                         hxis_store,hgammas_store,ibool,kappastore,myrank,DT,t0,tshift_src, &
                         coord,nspec,nglob,xigll,zigll,NPROC,xi_source, &
                         gamma_source,coorg,knods,NGNOD,npgeo,iglob_source,x_source,z_source, &
                         vx_source,vz_source, time_stepping_scheme, &
                         SOURCE_IS_MOVING, &
                         hxis,hpxis,hgammas,hpgammas

  use moving_sources_par, only: locate_source_moving

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob,ispec
  double precision :: hlagrange
  double precision :: xsrc,zsrc,timeval,t_used
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  if (time_stepping_scheme == 1) then
    ! Newmark
    timeval = (it-1)*DT
  else
    call exit_MPI(myrank,'Only Newmark time scheme is implemented for moving sources (2)')
  endif

  ! user output
  if ((myrank == 0) .and. (it == 1)) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*) 'Your are using acoustic moving source capabilities. Please cite:'
    write(IMAIN,*) 'Bottero (2018) Full-wave numerical simulation of T-waves and of moving acoustic sources'
    write(IMAIN,*) 'PhD thesis'
    write(IMAIN,*) 'https://tel.archives-ouvertes.fr/tel-01893011'
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Note: subroutine compute_add_sources_acoustic_moving_sources can be greatly'
    write(IMAIN,*) 'optimized. See what is done in init_moving_sources (in moving_sources_par.f90).'
    write(IMAIN,*) 'This is easy to do and would probably greatly improve the computational time'
    write(IMAIN,*)
    ! timing warning
    do i_source = 1,NSOURCES
      if ((abs(tshift_src(i_source)) > 0.0d0) .or. (abs(t0) > 0.0d0)) then
        write(IMAIN,*) 'Source #',i_source
        write(IMAIN,*) ' !! BEWARE !! Parameters tshift and/or t0 are used with moving source !'
        write(IMAIN,*) ' The time step for the moving source is: '
        write(IMAIN,*) '    t_used = (it_l-1)*DT-t0-tshift_src(i_source)'
        write(IMAIN,*) ' And the source position is calculated like:'
        write(IMAIN,*) '  xsrc = x_source + vx_source*t_used'
        write(IMAIN,*)
      endif
    enddo
  endif

  do i_source = 1,NSOURCES
    if (abs(source_time_function(i_source,it,i_stage)) > TINYVAL) then
      t_used = (timeval-t0-tshift_src(i_source))
     ! moves and re-locates sources along x and z axis
      xsrc = x_source(i_source) + vx_source(i_source)*t_used
      zsrc = z_source(i_source) + vz_source(i_source)*t_used

      ! collocated force source
      ! TODO: this would be more efficient compled with first guess as in init_moving_sources_GPU
      !call locate_source_moving(xsrc,zsrc, &
      !                   ispec_selected_source(i_source),islice_selected_source(i_source), &
      !                   NPROC,myrank,xi_source(i_source),gamma_source(i_source),.true.)
      call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                         xsrc,zsrc, &
                         ispec_selected_source(i_source),islice_selected_source(i_source), &
                         NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,NGNOD,npgeo, &
                         iglob_source(i_source),.true.)

      ! print *,ispec_selected_source(i_source) > nspec, "xmin:", &
      !               coord(1,ibool(1,1,ispec_selected_source(i_source))), &
      !               "xmax:", coord(1,ibool(NGLLX,1,ispec_selected_source(i_source)))
      ! define and store Lagrange interpolators (hxis,hpxis,hgammas,hpgammas) at all the sources
      !if (AXISYM) then
      !  if (is_on_the_axis(ispec_selected_source(i_source)) .and. myrank == islice_selected_source(i_source)) then
      !    call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
      !    !do j = 1,NGLJ ! ABAB same result with that loop, this is good
      !    !  hxis(j) = hglj(j-1,xi_source(i),xiglj,NGLJ)
      !    !enddo
      !  else
      !    call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
      !  endif
      !else
        call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
      !endif
      call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

      ! stores Lagrangians for source
      hxis_store(i_source,:) = hxis(:)
      hgammas_store(i_source,:) = hgammas(:)

    endif
  enddo

  ! adds source contributions
  do i_source = 1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    ! .and. acoustic(ispec_selected_source(i_source)) ??
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      if (ispec_is_acoustic(ispec)) then
        ! collocated force
        ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
        ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
        ! to add minus the source to Chi_dot_dot to get plus the source in pressure
        if (source_type(i_source) == 1) then
          ! forward wavefield
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
              sourcearray(1,i,j) = hlagrange

              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                      real(source_time_function(i_source,it,i_stage)*sourcearray(1,i,j) / &
                      kappastore(i,j,ispec),kind=CUSTOM_REAL)
            enddo
          enddo
          ! moment tensor
          else if (source_type(i_source) == 2) then
            call exit_MPI(myrank,'Cannot have moment tensor source in acoustic element')
        endif
      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic_moving_sources

!
!=====================================================================
!

! for acoustic solver for adjoint propagation wave field

  subroutine compute_add_sources_acoustic_adjoint()

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: potential_dot_dot_acoustic,ispec_is_acoustic,NSTEP,it, &
                         nrecloc,ispec_selected_rec_loc, &
                         ibool,source_adjoint,xir_store_loc,gammar_store_loc
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec
  integer :: it_tmp
  real(kind=CUSTOM_REAL) :: stf

  ! time step index
  it_tmp = NSTEP - it + 1

  do irec_local = 1,nrecloc

    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_acoustic(ispec)) then
      ! add source array
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          ! adjoint source of Peter et al. (A8):
          !   f^adj = - sum_i \partial_t^2 (p^syn - p^obs)(T-t) \delta(x - x_i)
          ! note that using the adjoint source derived from the optimization problem, there is no 1/kappa term applied
          ! to the adjoint source. the negative sign also is part of the construction of the adjoint source.
          !
          ! since we don't know which formulation of adjoint source is used for the input, we add the adjoint source as is,
          ! without 1/kappa factor, and with a positive sign.
          stf = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,1)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + stf

        enddo
      enddo
    endif ! if element acoustic
  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_acoustic_adjoint

!
!========================================================================
!

  double precision function get_stf_acoustic(t_used,i_source)

  ! prepares source_time_function array

  use constants, only: IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER, &
                       SOURCE_DECAY_MIMIC_TRIANGLE, &
                       C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: NSTEP,time_function_type,name_of_source_file,         &
                         burst_band_width,f0_source,tshift_src,factor,t0,      &
                         USE_TRICK_FOR_BETTER_PRESSURE, myrank

  implicit none

  ! local parameters
  double precision :: timeval, DecT, Tc, omegat, omega_coa,dummy_t,coeff, t_used, Nc
  double precision :: hdur,hdur_gauss

  integer :: it,i_source,ier,num_file

  character(len=150) :: error_msg1 = 'Error opening the file that contains the external source: '
  character(len=250) :: error_msg

  ! external functions
  double precision, external :: comp_source_time_function_heaviside_hdur
  double precision, external :: comp_source_time_function_Gaussian,comp_source_time_function_dGaussian, &
    comp_source_time_function_d2Gaussian,comp_source_time_function_d3Gaussian,marmousi_ormsby_wavelet,cos_taper
  double precision, external :: comp_source_time_function_Ricker,comp_source_time_function_d2Ricker

  double precision :: stf
  
  ! determines source_time_function value for different source types
  select case (time_function_type(i_source))
  case (1)
    ! Ricker: second derivative of a Gaussian
    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Second derivative of Ricker source time function :
      stf =  - factor(i_source) * &
               comp_source_time_function_d2Ricker(t_used,f0_source(i_source))
    else
      ! Ricker (second derivative of a Gaussian) source time function
      stf = - factor(i_source) * &
              comp_source_time_function_Ricker(t_used,f0_source(i_source))
    endif

  case (2)
    ! first derivative of a Gaussian

    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Third derivative of Gaussian source time function :
      stf = - factor(i_source) * &
                comp_source_time_function_d3Gaussian(t_used,f0_source(i_source))
    else
      ! First derivative of a Gaussian source time function
      stf = - factor(i_source) * &
                comp_source_time_function_dGaussian(t_used,f0_source(i_source))
    endif

  case (3,4)
    ! Gaussian/Dirac type

    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Second derivative of Gaussian :
      stf = - factor(i_source) * &
                 comp_source_time_function_d2Gaussian(t_used,f0_source(i_source))
    else
      ! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
      stf = - factor(i_source) * &
                  comp_source_time_function_Gaussian(t_used,f0_source(i_source))
    endif

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
    ! ocean acoustics type I
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
    ! ocean acoustics type II
    DecT = t0 + tshift_src(i_source)
    Tc = 4.d0 / f0_source(i_source) + DecT
    omega_coa = TWO*PI*f0_source(i_source)
    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Second derivative of source 7 :
      if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
        stf = factor(i_source) * &
                  0.5d0*(ONE-cos(omega_coa*t_used/4.0d0))*sin(omega_coa*t_used)
      else
        stf = ZERO
      endif
    else
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
    endif

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

    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Second derivative of Burst :
      if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
        stf = - factor(i_source) * (0.5d0 * (omega_coa)**2 * &
                  sin(omega_coa*t_used) * cos(omega_coa*t_used/Nc) / Nc**2 - &
                  0.5d0 * (omega_coa)**2 * sin(omega_coa*t_used) * &
                  (ONE-cos(omega_coa*t_used/Nc)) + &
                  (omega_coa)**2 * cos(omega_coa*t_used) * &
                  sin(omega_coa*t_used/Nc) / Nc)
      !else if (timeval > DecT) then
      !  stf = ZERO
      else
        stf = ZERO
      endif
      ! Integral of burst
      !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
      !  stf = - factor(i_source) * ( &
      !  Nc*( (Nc+1.0d0)*cos((omega_coa*(Nc-1.0d0)*t_used)/Nc) + &
      !       (Nc-1.0d0)*cos((omega_coa*(Nc+1.0d0)*t_used)/Nc)) - &
      !  TWO*(Nc**2-1.0d0)*cos(omega_coa*t_used) &
      !  ) / (8.0d0*PI*f0_source(i_source)*(Nc-1)*(Nc+1))
      !else
      !  stf = ZERO
      !endif
      ! Double integral of burst
      !if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
      !  stf = - factor(i_source) * ( &
      !      -sin(TWO*f0_source(i_source)*Pi*t_used)/(8.0d0*f0_source(i_source)**TWO*Pi**2) + &
      !      (Nc**2*sin((TWO*f0_source(i_source)*(Nc-1)*PI*t_used)/Nc))/(16.0d0*f0_source(i_source)**2*(Nc-1)**2*Pi**2) + &
      !      (Nc**2*sin((TWO*f0_source(i_source)*(Nc+1)*PI*t_used)/Nc))/(16.0d0*f0_source(i_source)**2*(Nc+1)**2*Pi**2) )
      !else
      !  stf = ZERO
      !endif
    else
      if (timeval > DecT .and. timeval < Tc) then ! t_used > 0 t_used < Nc/f0_source(i_source)) then
        stf = - factor(i_source) * &
                  0.5d0*(ONE-cos(omega_coa*t_used/Nc))*sin(omega_coa*t_used)
      !else if (timeval > DecT) then
      !  stf = ZERO
      else
        stf = ZERO
      endif
    endif

  case (10)
    ! Sinus source time function
    omega_coa = TWO*PI*f0_source(i_source)
    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! Third derivative of Gaussian source time function :
      stf = -TWO*PI*omega_coa*f0_source(i_source)* &
                                                  factor(i_source) * sin(omega_coa*t_used)
    else
      ! First derivative of a Gaussian source time function
      stf = factor(i_source) * sin(omega_coa*t_used)
    endif

  case (11)
      ! Marmousi_ormsby_wavelet
      hdur = 1.0 / 35.0
      stf = factor(i_source) * &
                cos_taper(t_used,hdur) * marmousi_ormsby_wavelet(PI*t_used)

  case default
    call exit_MPI(myrank,'unknown source time function')

  end select

  ! return value
  get_stf_acoustic = stf

  end function get_stf_acoustic
