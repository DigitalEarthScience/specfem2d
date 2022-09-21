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

! for viscoelastic solver

  subroutine compute_add_sources_viscoelastic(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,myrank
  use constants, only: C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: P_SV,ispec_is_elastic,nglob_elastic,NSOURCES,         &
                         islice_selected_source,ispec_selected_source,         &
                         sourcearrays,ibool,tshift_src,t0,DT,                  &
                         time_stepping_scheme,islice_selected_source,myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage
  double precision,external :: get_stf_viscoelastic

  !local variable
  integer :: i_source,i,j,iglob,ispec
  !real(kind=CUSTOM_REAL) :: stf_used
  double precision :: stf_used
  double precision :: timeval,t_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

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
        stf_used = get_stf_viscoelastic(t_used,i_source)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              ! 2D: x-component uses array(1,..) and z-component (2,..)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(1,i,j,i_source) * stf_used  ! x-direction
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + sourcearrays(2,i,j,i_source) * stf_used  ! z-direction
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              ! 2D: y-component uses array(1,..)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(1,i,j,i_source) * stf_used  ! y-direction
            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic

!
!=====================================================================
!

  subroutine compute_add_sources_viscoelastic_moving_sources(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,TINYVAL,NGLJ,IMAIN
  use constants, only: C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool,coord,nspec,nglob,xigll,zigll,NPROC, &
                         xi_source,gamma_source,coorg,knods,NGNOD,npgeo,iglob_source,x_source,z_source, &
                         vx_source,vz_source,DT,t0,myrank, &
                         time_stepping_scheme,hxis_store,hgammas_store,tshift_src,source_type,ispec_is_acoustic, &
                         hxis,hpxis,hgammas,hpgammas,anglesource,ispec_is_poroelastic,Mxx,Mxz,Mzz,gammax,gammaz,xix,xiz, &
                         AXISYM,xiglj,is_on_the_axis,initialfield,SOURCE_IS_MOVING

  use moving_sources_par, only: locate_source_moving

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage
  double precision,external :: get_stf_viscoelastic

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used
  double precision :: hlagrange
  double precision :: xsrc,zsrc,timeval,t_used
  ! single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  if (time_stepping_scheme == 1) then
    ! Newmark
    timeval = (it-1)*DT
  else
    call exit_MPI(myrank,'Only Newmark time scheme is implemented for moving sources (3)')
  endif

  if ((myrank == 0) .and. (it == 1)) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*) 'Your are using elastic moving source capabilities. Please cite:'
    write(IMAIN,*) 'Bottero (2018) Full-wave numerical simulation of T-waves and of moving acoustic sources'
    write(IMAIN,*) 'PhD thesis'
    write(IMAIN,*) 'https://tel.archives-ouvertes.fr/tel-01893011'
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Note: subroutine compute_add_sources_viscoelastic_moving_sources can be greatly'
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
      ! moves and re-locates sources along x and z-axis
      xsrc = x_source(i_source) + vx_source(i_source)*t_used
      zsrc = z_source(i_source) + vz_source(i_source)*t_used
      ! collocated force source
      if (source_type(i_source) == 1) then
        ! TODO: this would be more efficient compled with first guess as in init_moving_sources_GPU()
        !call locate_source_moving(xsrc,zsrc, &
        !                   ispec_selected_source(i_source),islice_selected_source(i_source), &
        !                   NPROC,myrank,xi_source(i_source),gamma_source(i_source),.true.)
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           xsrc,zsrc, &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,NGNOD,npgeo, &
                           iglob_source(i_source),.true.)

      else if (source_type(i_source) == 2) then
        ! moment-tensor source
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           xsrc,zsrc, &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,NGNOD,npgeo, &
                           iglob_source(i_source),.false.)

      else if (.not. initialfield) then

        call exit_MPI(myrank,'incorrect source type')

      endif

      ispec = ispec_selected_source(i_source)
      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! Lagrange interpolators
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
            !do j = 1,NGLJ ! ABAB same result with that loop, this is good
            !  hxis(j) = hglj(j-1,xi_source(i_source),xiglj,NGLJ)
            !enddo
          else
            call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
          endif
        else
          call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
        endif
        call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

        if (mod(it,10000) == 0) then
            !  write(IMAIN,*) "myrank:",myrank
            ! user output
            if (myrank == islice_selected_source(i_source)) then
              iglob = ibool(2,2,ispec_selected_source(i_source))
              !write(IMAIN,*) 'xcoord: ',coord(1,iglob)
              write(IMAIN,*) 'Problem... it??: ',it,'xcoord: ',coord(1,iglob)," iglob",iglob
              !'source carried by proc',myrank,"  source x:",x_source(i_source)," ispec:",ispec_selected_source(i_source)

              !call flush_IMAIN()
            endif

        endif

        ! stores Lagrangians for source
        hxis_store(i_source,:) = hxis(:)
        hgammas_store(i_source,:) = hgammas(:)

        sourcearray(:,:,:) = 0._CUSTOM_REAL

        ! computes source arrays
        select case (source_type(i_source))
        case (1)
          ! collocated force source
          do j = 1,NGLLZ
            do i = 1,NGLLX
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

              ! source element is acoustic
              if (ispec_is_acoustic(ispec)) then
                sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
              endif

              ! source element is elastic
              if (ispec_is_elastic(ispec)) then
                if (P_SV) then
                  ! P_SV case
                  sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                  sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                else
                  ! SH case (membrane)
                  sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
                endif
              endif

              ! source element is poroelastic
              if (ispec_is_poroelastic(ispec)) then
                sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
              endif

            enddo
          enddo

        case (2)
          ! moment-tensor source
          call compute_arrays_source(ispec,xi_source(i_source),gamma_source(i_source),sourcearray, &
                                     Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)
          ! checks source
          if (ispec_is_acoustic(ispec)) then
            call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
          endif

          ! checks wave type
          if (ispec_is_elastic(ispec)) then
            if (.not. P_SV ) call exit_MPI(myrank,'cannot have moment tensor source in SH (membrane) waves calculation')
          endif

        end select

        ! stores sourcearray for all sources
        sourcearrays(:,:,:,i_source) = sourcearray(:,:,:)

      endif
    endif
  enddo

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

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
        stf_used = get_stf_viscoelastic(t_used,i_source)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(1,i,j,i_source) * stf_used
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + sourcearrays(2,i,j,i_source) * stf_used
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(1,i,j,i_source) * stf_used

              ! debugging source contribution
              !if (iglob == 37905) &
              !write(1234,*) it, dble(sourcearrays(1,i,j,i_source) * source_time_function(i_source,it,i_stage)), &
              !              accel_elastic(1,iglob),source_time_function(i_source,it,i_stage),sourcearrays(1,i,j,i_source)
            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic_moving_sources

!
!=====================================================================
!

! for viscoelastic solver for adjoint propagation wave field

  subroutine compute_add_sources_viscoelastic_adjoint()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: P_SV,accel_elastic,ispec_is_elastic,NSTEP,it, &
                         nrecloc,ispec_selected_rec_loc,ibool, &
                         source_adjoint,xir_store_loc,gammar_store_loc
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec
  integer :: it_tmp
  real(kind=CUSTOM_REAL) :: stfx,stfz

  ! time step index
  it_tmp = NSTEP - it + 1

  do irec_local = 1,nrecloc

    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_elastic(ispec)) then
      ! add source array
      if (P_SV) then
        ! P-SV waves
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            stfx = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,1)
            stfz = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,2)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) + stfx
            accel_elastic(2,iglob) = accel_elastic(2,iglob) + stfz
          enddo
        enddo
      else
        ! SH (membrane) wavescompute_forces_v
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            stfx = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,1)

            accel_elastic(1,iglob) = accel_elastic(1,iglob) + stfx
          enddo
        enddo
      endif
    endif ! if element is elastic

  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_viscoelastic_adjoint

!
!========================================================================
!

  double precision function get_stf_viscoelastic(t_used,i_source)

  ! prepares source_time_function array

  use constants, only: IMAIN,ZERO,ONE,TWO,HALF,PI,QUARTER, &
                       SOURCE_DECAY_MIMIC_TRIANGLE, &
                       C_LDDRK,C_RK4,ALPHA_SYMPLECTIC

  use specfem_par, only: NSTEP, &
                         time_function_type, name_of_source_file, burst_band_width,f0_source,tshift_src, &
                         factor,t0,&
                         myrank

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
    ! ocean viscoelastics type I
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
    ! ocean viscoelastics type II
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
  get_stf_viscoelastic = stf

  end function get_stf_viscoelastic
