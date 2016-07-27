! include elicense.txt'
! this is a main routine for multistage excavation
! REVISION:
!   HNG, Aug 25,2011; HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
subroutine semimage2d(ismpi,gnod,sum_file,format_str)
! import necessary libraries
use global
use string_library, only : parse_file,get_quoted_ivec,get_quoted_rvec
use math_constants
use gll_library
!use mesh_spec
use shape_library
use math_library
!use gauss_library
use global_dof
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
use sparse
use parsolver
use parsolver_petsc
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif
use bc
use blur
use preprocess
use visual
use postprocess
use ieee_arithmetic

implicit none
logical,intent(in) :: ismpi ! .true.=parallel version, .false.=serial version
integer,intent(in) :: gnod(4)!geometrical nodes (corner nodes) per element
character(len=250),intent(in) :: sum_file !Summary file. Not currently used
character(len=20),intent(in) :: format_str !ptail=string with processor
!number. format_str determines format for processor number (eg. 1, 01, 001)

integer :: i,ios,istat,j,k !i,j,k are dummy vars for interation
!ios:input output status indicator (for opening files)
!istat: status indicator for allocation (can be used in other contexts)

integer :: i_elmt,i_nliter,i_node,i_inc,i_srf,i_tstep,i_excav !do-loop indices
integer :: ielmt,igdof,imat,inode,iedof !element ID for gdof, node, etc.

real(kind=kreal),parameter :: r3=3.0_kreal,one_third=one/r3,two_third=two/r3
!constants

real(kind=kreal) :: t,uerr,umax,uxmax,xval !uerr: used to check convergence
!umax: max of displacement magnitude, uxmax: max of displacement components

integer :: cg_iter,cg_tot,nl_iter,nl_tot !cg_iter: conjugate gradient iteration,
!cg_tot: total cg interation, nl_iter: nonlinear iteration, nl_tot: total nl
!iteration
logical :: nl_isconv ! logical variable to check convergence of
! nonlinear (NL) iterations

! dynamic arrays
integer,allocatable::gdof_elmt(:,:),num(:),nvalency(:) !num: g_num for
!particular element. gdof_elmt: matrix of elemental degrees of freedom (per 
!processor). eg. gdof_elmt(:,1) gives all dof in 1st element in the processsor.
!nvalency: number of elements that share each node

real(kind=kreal),allocatable :: gload(:),bload(:),coord(:,:),      &
dprecon(:),ndscale(:),x(:)
real(kind=4),allocatable :: ggll(:,:)
real(kind=kreal),allocatable :: storekmat(:,:,:)
!bload: load computed on all nodes of the element
!coord: coordinates of geometrical nodes
!deriv: derivative of interpolation functions
!dprecon: diagonal preconditioner. Not used for petsc solver
!ndscale: nondimensionalize scale
!nodalnu: nodal displacement for all nodes (not just BC)
!storekmat: stiffness matrix for all elements
!x: solution
!elas_strain: elastic strain for all elements
!stress_local: stress for each element
!stress_global: nodal stress for all elements in processor
!scf: not used
!vmeps: von mises strain

integer,allocatable :: egdof(:) ! placeholder array. holds values of gdof_elmt for a given element.

integer :: map2exodus(4),node_quad4(4)
real(kind=kreal),allocatable :: dshape_quad4(:,:,:)

!double precision
real(kind=kreal),allocatable :: interpfgll(:,:)

character(len=250) :: inp_fname,out_fname,prog
character(len=150) :: path
character(len=20) :: ensight_etype
character(len=80) :: buffer,destag ! this must be 80 characters long
character(len=20) :: ext 
character(len=250) :: case_file,geo_file 
integer :: npart,tinc,tstart,twidth,ts !ets: time set for ensight gold

real(kind=kreal) :: cpu_tstart,cpu_tend,telap,step_telap,max_telap,mean_telap

logical :: gravity,pseudoeq ! gravity load and pseudostatic load
real(kind=kreal),allocatable :: wpressure(:) ! water pressure
logical,allocatable :: submerged_node(:)

real(kind=kreal),allocatable:: nodalg(:,:),nodalcg(:,:),nodalm0(:,:),nodalcm0(:,:),nodalccm0g(:,:),nodalnu0(:,:)
real(kind=kreal),allocatable:: nodalm(:,:),nodalnu(:,:)
real(kind=kreal) :: energy 
logical :: discardm
! blurring
integer :: iup,nup
real(kind=kreal) :: gauss_sigma 
real(kind=kreal) :: fac,halfinvs2,xsq
real(kind=kreal) :: xi(ndim),xj(ndim) 
real(kind=kreal),allocatable :: nodalgauss(:)

! tomo model
integer :: ig,ig1,ig2,ig3,ig4
integer :: ix1,ix2,iy1,iy2
integer :: grid_l1,grid_l2,grid_m1,grid_m2,grid_n1,grid_n2
integer :: grid_n,grid_nx,grid_ny,grid_nz
integer :: grid_wext(6)
real(kind=kreal) :: grid_x0(3),grid_dx(3)
real(kind=kreal),allocatable :: grid_m(:,:)
real(kind=kreal) :: a,b,xp,yp,zp,shape1,shape2,shape3,shape4,xg(2,4)
logical,allocatable :: isnode(:)
character(len=1) :: dumc,dumc0
character(len=60) :: dumstr
character(len=250) :: line
real :: rval
real(kind=kreal),allocatable :: bulkmod(:),shearmod(:)

real(kind=kreal) :: maxx


integer :: ipart
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode
integer :: tot_neq,max_neq,min_neq
integer :: maxngnode
! number of active ghost partitions for a node
integer,allocatable :: ngpart_node(:)
character(len=250) :: errtag ! error message
integer :: errcode
logical :: isopen ! flag to check whether the file is opened

!needed for petsc
character(len=20) :: proc_str
character(len=80) :: ofname


errtag=""; errcode=-1

ipart=myid-1 ! partition id starts from 0 (same as myrank)

ensight_etype='quad4'

! map sequential node numbering to exodus/cubit order for 8-noded hexahedra
map2exodus=(/ 1,2,4,3 /)

! prepare ghost partitions for the communication
if(nproc.gt.1)then
  call prepare_ghost(myid,nproc,ngpart,maxngnode)
endif

!-------------------------------

if(myid==1)then
  write(stdout,'(a)')'setting model properties...'
  if(mattype==0)then
    write(stdout,'(a)')'  model type: BLOCK'
  elseif(mattype==1)then
    write(stdout,'(a)')'  model type: TOMO'
  else
    write(stdout,'(a,i3)')'  ERROR: wrong model type',mattype
    stop
  endif
endif

call allocate_store_arrays

allocate(interpfgll(ngll,ngll))
call precompute_quadrature(gnod,interpfgll)


! compute node valency only once
allocate(num(nenod),nvalency(nnode))

nvalency=0
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,ielmt)
  nvalency(num)=nvalency(num)+1
enddo

! model properties
allocate(isnode(nnode))
allocate(nodalg(nndof,nnode),nodalcg(nndof,nnode),nodalcm0(nndof,nnode),nodalccm0g(nndof,nnode))
nodalg=ZERO

if(model_input.eq.1)then
  ! tomo model
  open(100,file=trim(mat_path)//trim(matfile),action='read',form='unformatted',status='old')
  allocate(ggll(ngll,nelmt))
  read(100)ggll
  do i_elmt=1,nelmt
    nodalg(1,g_num(:,i_elmt))=nodalg(1,g_num(:,i_elmt))+ggll(:,i_elmt)
  enddo
  deallocate(ggll)
  nodalg(1,:)=nodalg(1,:)/nvalency
  close(100)
else
  if(mattype==0)then
    ! block model
    ! Bulk modulus and Shear modulus
    allocate(bulkmod(nmat))
    ! block model
    bulkmod(1)=6000.0_kreal
    bulkmod(2)=4000.0_kreal
    bulkmod(3)=2000.0_kreal
    do i_elmt=1,nelmt
      num=g_num(:,i_elmt)
      nodalg(:,num)=nodalg(:,num)+bulkmod(mat_id(i_elmt)) 
    enddo
    nodalg(1,:)=nodalg(1,:)/nvalency
  else
    ! tomo model
    open(100,file=trim(mat_path)//trim(matfile),action='read',status='old')
    read(100,*)
    read(100,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    call get_quoted_ivec(line,"WholeExtent",grid_wext,6)
    call get_quoted_rvec(line,"Origin",grid_x0,3)
    call get_quoted_rvec(line,"Spacing",grid_dx,3)
    grid_nx=grid_wext(2)-grid_wext(1)+1
    grid_ny=grid_wext(4)-grid_wext(3)+1
    grid_n=grid_nx*grid_ny
    read(100,*)
    read(100,*)
    read(100,*)
    allocate(grid_m(1,grid_n))
    grid_m=zero
  !  ig=1
  !  do i=1,grid_n/2
  !    read(100,'(a)',iostat=ios)line
  !    !print*,'line:',i,line
  !    read(line,*)grid_m(:,ig),grid_m(:,ig+1)
  !    ig=ig+2
  !  enddo
    do i=1,grid_n
      read(100,*)grid_m(:,i)
    enddo
    close(100)
    grid_l1=grid_wext(1)+1
    grid_l2=grid_wext(2)+1
    grid_m1=grid_wext(3)+1
    grid_m2=grid_wext(4)+1
    !print*,'grid_l:',grid_l1,grid_l2,grid_m1,grid_m2
    isnode=.false.
    do i_node=1,nnode
      xp=g_coord(1,i_node)
      yp=g_coord(2,i_node)

      ix1=floor((xp-grid_x0(1))/grid_dx(1))+1
      iy1=floor((yp-grid_x0(2))/grid_dx(2))+1
      !print*,'ix:',ix1,iy1
      if(ix1.le.grid_l1)ix1=grid_l1
      if(iy1.le.grid_m1)iy1=grid_m1
      if(ix1.ge.grid_l2)ix1=grid_l2-1
      if(iy1.ge.grid_m2)iy1=grid_m2-1
      ix2=ix1+1
      iy2=iy1+1

      ig1=(iy1-1)*grid_nx+ix1
      ig2=(iy1-1)*grid_nx+ix2
      ig3=(iy2-1)*grid_nx+ix2
      ig4=(iy2-1)*grid_nx+ix1

      !print*,'node:',i_node,grid_x0,grid_dx
      !print*,'xp:',xp,yp
      !print*,'ix:',ix1,iy1
      !print*,'ig:',ig1,ig2,ig3,ig4
      !print*,grid_m(:,ig1),grid_m(:,ig2),grid_m(:,ig3),grid_m(:,ig4)
      !if(i_node>10)stop
      ! define shape functions for 4 corner points
      ! JN Reddy, Introduction to Finite Elemenet method. P423
      xg(1,1)=grid_x0(1)+grid_dx(1)*(ix1-1)
      xg(2,1)=grid_x0(2)+grid_dx(2)*(iy1-1)
      
      xg(1,2)=xg(1,1)+grid_dx(1)
      xg(2,2)=xg(2,1)
      
      xg(1,3)=xg(1,2)
      xg(2,3)=xg(2,2)+grid_dx(2)
      
      xg(1,4)=xg(1,1)
      xg(2,4)=xg(2,3)

      a=xg(1,2)-xg(1,1)
      b=xg(2,4)-xg(2,1)

      ! shift origin to the cell corner
      xp=xp-xg(1,1)
      yp=yp-xg(2,1)
      xg(1,:)=xg(1,:)-xg(1,1)
      xg(2,:)=xg(2,:)-xg(2,1)

      shape1=(one-xp/a)*(one-yp/b)
      shape2=xp*(one-yp/b)/a
      shape3=xp*yp/(a*b)
      shape4=(one-xp/a)*yp/b

      !print*,'shape:',shape1+shape2+shape3+shape4
     

      ! interpolation
      nodalg(:,i_node)=shape1*grid_m(:,ig1)+shape2*grid_m(:,ig2)+       &
                       shape3*grid_m(:,ig3)+shape4*grid_m(:,ig4)
      !nodalg(:,i_node)=shape1*norm(grid_m(:,ig1))+shape2*norm(grid_m(:,ig2))+       &
      !                shape3*norm(grid_m(:,ig3))+shape4*norm(grid_m(:,ig4))
      

    enddo
  endif
endif
if(nproc.gt.1)then
  call assemble_ghosts_nodal(myid,ngpart,maxngnode,nndof,nodalg,nodalg)
endif
call sync_process() 

call save_nodal_data(ptail,format_str,0,nnode,nodalg,'Original model','g',.false.)
write(stdout,'(a)')'complete!'

!gauss_sigma=4.0_kreal
!write(stdout,'(a)')'computing nodal distances...'
!! number of upper trianglular elements
!nup=nnode*(nnode-1)/2
!call sync_process() 
!allocate(nodalgauss(0:nup),stat=ios)
!if(ios.ne.0)then
!  print*,'ERROR: not enough memory!'
!  stop
!endif
!
!call compute_gaussian_kernel(gauss_sigma,nodalgauss)
!write(stdout,'(a)')'complete!'
!
!! blur original image
!
!write(stdout,'(a)')'blurring model...'
!!call convolve_with_gaussian(ismpi,gnod,gauss_sigma,nodalg,nodalcg,errcode,errtag)
!call convolve_with_gaussian_precomp(ismpi,gnod,nup,nodalgauss,nodalg,nodalcg,errcode,errtag)
!if(nproc.gt.1)then
!  call assemble_ghosts_nodal(myid,ngpart,maxngnode,nndof,nodalcg,nodalcg)
!endif
!!
!! set blurred image as the background image
!nodalg=nodalcg
!call save_nodal_data(ptail,format_str,0,nnode,nodalg,'Blurred model','b',.false.)
!write(stdout,'(a)')'complete!'
!print*,minval(grid_m),maxval(grid_m)
!print*,minval(nodalg),maxval(nodalg)
!!stop
!!-------------------------------------------------------------------------------

! apply displacement boundary conditions
if(myid==1)write(stdout,'(a)',advance='no')'applying BC...'
allocate(gdof(nndof,nnode),stat=istat)
if (istat/=0)then
  write(stdout,*)'ERROR: cannot allocate memory!'
  stop
endif

gdof=1
call modify_gdof(gdof,neq,errcode,errtag)
if(myid==1)write(stdout,*)'complete!'
!-------------------------------------------------------------------------------


! modify ghost GDOFs
if(nproc.gt.1)then
  call prepare_ghost_gdof(myid,ngpart,gdof)
endif

allocate(coord(ngnode,ndim),egdof(nedof),stat=istat)
if (istat/=0)then
  write(stdout,*)'ERROR: cannot allocate memory!'
  stop
endif
!-------------------------------------------------------------------------------

! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than gdof!!!
allocate(gdof_elmt(nedof,nelmt))
gdof_elmt=0
do i_elmt=1,nelmt
  gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) !g=g_g(:,i_elmt)
enddo

!-------------------------------------------------------------------------------
! global indexing
! this process is necessary for petsc implementation and is done in a single
! processor
call sync_process

write(proc_str,'(i10)')myid
ofname='partition/before'//trim(adjustl(proc_str))
open(22,file=ofname,action='write',status='replace')
write(22,*)g_num
close(22)

if(nproc.gt.1.and.solver_type.eq.petsc_solver .and. myid.eq.1)call gindex()
call sync_process

ofname='partition/after'//trim(adjustl(proc_str))
open(22,file=ofname,action='write',status='replace')
write(22,*)g_num
close(22)

!-------------------------------------------------------------------------------    

tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
if(myid==1)then
  write(stdout,*)'degrees of freedoms => total:',tot_neq,' max:',max_neq,      &
  ' min:',min_neq
endif


if(myid==1)write(stdout,'(a)',advance='no')'preprocessing...'

! open summary file
open(unit=10,file=trim(sum_file),status='old',position='append',action='write',&
iostat=ios)
write(10,*)'CG_MAXITER, CG_TOL, NL_MAXITER, NL_TOL'
write(10,*)cg_maxiter,cg_tol,nl_maxiter,nl_tol

if(myid==1)then
  write(stdout,'(a,e12.4,a,i5)')'CG_TOL:',cg_tol,' CG_MAXITER:',cg_maxiter
  write(stdout,'(a,e12.4,a,i5)')'NL_TOL:',nl_tol,' NL_MAXITER:',nl_maxiter
  write(stdout,'(/,a)')'--------------------------------------------'
endif

write(10,*)'Number of time steps'
write(10,*)ntstep
write(10,*)'STEP, CGITER, NLITER, UXMAX, UMAX'
close(10)

allocate(bload(0:neq),gload(0:neq),x(0:neq),dprecon(0:neq),stat=istat)

allocate(storekmat(nedof,nedof,nelmt),stat=istat) 
if (istat/=0)then
  write(stdout,*)'ERROR: cannot allocate memory for storekmat!'
  stop
endif

allocate(ngpart_node(nnode))

if(myid.eq.1)then
  print*,'Image parameter:',eps_img,gam_img,eta_img
endif

if(solver_type.eq.builtin_solver .and.solver_diagscale)allocate(ndscale(0:neq))
allocate(nodalnu0(1,nnode),nodalm0(1,nnode))
allocate(nodalnu(1,nnode),nodalm(1,nnode))

open(111,file=trim(out_path)//'energy,dat',action='write',status='replace')
! initialize
nodalnu0=one
nodalm0=nodalg

!! compute h*g once for all
!!  gauss_sigma=5.0_kreal
!nodalcg=zero
!!call convolve_with_gaussian(ismpi,gnod,gauss_sigma,nodalg,nodalcg,errcode,errtag)
!call convolve_with_gaussian_precomp(ismpi,gnod,nup,nodalgauss,nodalg,nodalcg,errcode,errtag)
!
!! covolution of two gaussians
!call compute_gaussian_kernel(sqrt(two)*gauss_sigma,nodalgauss)

if(solver_type.eq.petsc_solver)then
  ! prepare sparsity of the stiffness matrix
  call prepare_sparse(nup)!,nodalgauss)

  ! petsc solver
  call petsc_initialize()
  call petsc_create_vector()                                                   
  call petsc_matrix_preallocate_size()                                         
  call petsc_create_matrix()                                                   
  call petsc_create_solver()
  if(myid==1)print*,'petsc_initialize: SUCCESS!'
endif

discardm=.true. ! .true. for reconstruction (convolved m term), .false. for segementation
print*,'-----------------------------------------------------------------------'
time_step: do i_tstep=1,ntstep
  write(*,'(a,a,i4,a,i4)')CR,'step:',i_tstep, '/',ntstep  

  !---------------------edge detection------------------------------------------
  if(myid.eq.1)print*,'staring nu solver.....'
  ! initialize variable
  nodalnu=zero
  ! compute stiffness, load, diagonal
  call stiffness_load_edgefield(gdof_elmt,interpfgll,nodalm0,storekmat,dprecon,bload)
  call sync_process()
  print*,'stiffness comutation finished!'
  
  ! diagonal preconditioner
  if(solver_type.eq.builtin_solver)then
    ! assemble diagonal preconditioner
    if(nproc.gt.1)then
      call assemble_ghosts(myid,ngpart,maxngnode,nndof,neq,dprecon,dprecon)
    endif
    call sync_process() 
    dprecon(0)=zero
    
    print*,'nzero in dprecon:',count(dprecon(1:).eq.zero),minval(abs(dprecon(1:))),maxval(abs(dprecon(1:)))
    
    if(solver_diagscale)then
      ! regularize linear equations,
      ndscale=ONE
      do i=1,neq
        ndscale(i)=one/sqrt(abs(dprecon(i)))
      enddo
      ! nondimensionalize stiffnes martix
      ! elastic region
      do i_elmt=1,nelmt
        ielmt=i_elmt
        egdof=gdof_elmt(:,ielmt)
        do i=1,nedof
          do j=1,nedof
            storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*ndscale(egdof(j))
          enddo
        enddo
      enddo
    else
      dprecon(1:)=one/dprecon(1:)
    endif
  endif

 ! solver
  cg_tot=0
  x=zero
  
  bload(0)=zero;
  if(myid.eq.1)print*,'Residual NL:',maxval(abs(bload))
  
  if(solver_type.eq.builtin_solver)then
    ! builtin solver
    call sync_process()
    if(solver_diagscale)then
      ! nondimensionalize load
      bload=ndscale*bload
      call cg_solver(myid,ngpart,maxngnode,neq,nelmt,storekmat,x,bload, &
      gdof_elmt,cg_iter,errcode,errtag)
      x=ndscale*x
    else
      ! pcg solver
      call pcg_solver(myid,ngpart,maxngnode,neq,nelmt,storekmat,x,bload,dprecon, &
      gdof_elmt,cg_iter,errcode,errtag)
    endif
  else
    !petsc solver
    call petsc_set_stiffness_matrix(storekmat)
    if(myid==1)print*,'petsc_set_stiffness_matrix: SUCCESS!'
    call petsc_set_vector(bload)
    if(myid==1)print*,'petsc_set_vector: SUCCESS!','load:',maxval(abs(bload))
    call petsc_set_ksp_operator()
    call petsc_solve(x(1:),cg_iter)
    if(myid==1)print*,'petsc_solve: SUCCESS!'    
  endif
  if(errcode/=0)call error_stop(errtag,stdout,myid)
  cg_tot=cg_tot+cg_iter
  x(0)=zero
  maxx=maxvec(abs(x))
  if(myid.eq.1)print*,'cg iters:',cg_iter,' max x:',maxx

  call sync_process()
  ! update total nodal displacement
  do i=1,nndof
    do j=1,nnode
      if(gdof(i,j)/=0)then
        nodalnu(i,j)=x(gdof(i,j)) ! time steps are not incremental!!
      endif
    enddo
  enddo
  call save_nodal_data(ptail,format_str,i_tstep,nnode,nodalnu,'Edge detector','nu')
  print*,'complete!'
  ! update
  nodalnu0=nodalnu
  !---------------------end edge detection--------------------------------------

!  !---------------------begin model computation---------------------------------
!  if(myid.eq.1)print*,'starting m solver.....'
!  ! initialize variable
!  nodalm=zero;
! 
!  !nodalcg=nodalg
!  ! compute load contributed by the term h*g
!  discardm=.false.
!  call compute_bodyload_model(gdof_elmt,interpfgll,nodalcg,gload)
!!  ! compute stiffness, load, diagonal
!!  discardm=.false.
!!  call stiffness_load_model(gdof_elmt,interpfgll,nodalcg,nodalnu0,gauss_sigma,nup,nodalgauss,storekmat,dprecon,bload,discardm)
!!  call sync_process()
!!  print*,'stiffness computation finished!'
!  
!  ! diagonal preconditioner
!  if(solver_type.eq.builtin_solver)then
!    ! assemble diagonal preconditioner
!    if(nproc.gt.1)then
!      call assemble_ghosts(myid,ngpart,maxngnode,nndof,neq,dprecon,dprecon)
!    endif
!    call sync_process() 
!    dprecon(0)=zero
!    
!    print*,'nzero in dprecon:',count(dprecon(1:).eq.zero),minval(abs(dprecon(1:))),maxval(abs(dprecon(1:)))
!    
!    if(solver_diagscale)then
!      ! regularize linear equations,
!      ndscale=ONE
!      do i=1,neq
!        ndscale(i)=one/sqrt(abs(dprecon(i)))
!      enddo
!      ! nondimensionalize stiffnes martix
!      ! elastic region
!      do i_elmt=1,nelmt
!        ielmt=i_elmt
!        egdof=gdof_elmt(:,ielmt)
!        do i=1,nedof
!          do j=1,nedof
!            storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*ndscale(egdof(j))
!          enddo
!        enddo
!      enddo
!    else
!      dprecon(1:)=one/dprecon(1:)
!    endif
!  endif
!
! ! solver
!  cg_tot=0
!  x=zero
!  
!  gload(0)=zero
!  bload=gload  
!  if(myid.eq.1)print*,'Residual NL:',maxval(abs(bload))
!  
!  if(solver_type.eq.builtin_solver)then
!    ! builtin solver
!    call sync_process()
!    if(solver_diagscale)then
!      ! nondimensionalize load
!      bload=ndscale*bload
!      call cg_solver(myid,ngpart,maxngnode,neq,nelmt,storekmat,x,bload, &
!      gdof_elmt,cg_iter,errcode,errtag)
!      x=ndscale*x
!    else
!      ! pcg solver
!      call pcg_solver(myid,ngpart,maxngnode,neq,nelmt,storekmat,x,bload,dprecon, &
!      gdof_elmt,cg_iter,errcode,errtag)
!    endif
!  else
!    !petsc solver
!    !call petsc_set_stiffness_matrix(storekmat)
!    call petsc_set_stiffness_matrixconv(nup,nvalency,interpfgll,gdof_elmt,nodalnu0,nodalgauss)
!    if(myid==1)print*,'petsc_set_stiffness_matrix: SUCCESS!'
!    call petsc_set_vector(bload)
!    if(myid==1)print*,'petsc_set_vector: SUCCESS!','load:',maxval(abs(bload))
!    call petsc_set_ksp_operatorconv()
!    call petsc_solve(x(1:),cg_iter)
!    if(myid==1)print*,'petsc_solve: SUCCESS!'    
!  endif
!  !if(errcode/=0)call error_stop(errtag,stdout,myid)
!  if(errcode/=0)print*,'WARNING:',trim(errtag)
!  cg_tot=cg_tot+cg_iter
!  x(0)=zero
!  maxx=maxvec(abs(x))
!  if(myid.eq.1)print*,'cg iters:',cg_iter,' max x:',maxx
!
!  call sync_process()
!  ! update total nodal displacement
!  do i=1,nndof
!    do j=1,nnode
!      if(gdof(i,j)/=0)then
!        nodalm(i,j)=x(gdof(i,j)) ! time steps are not incremental!!
!      endif
!    enddo
!  enddo
!  print*,'hey!',minval(nodalm),maxval(nodalm)
!  do i=1,nnode
!    if(ieee_is_nan(nodalm(1,i)).or. .not.ieee_is_finite(nodalm(1,i)))then
!      print*,'whow!:',nodalm(i,j)
!      stop
!    endif
!  enddo
!  call save_nodal_data(ptail,format_str,i_tstep,nnode,nodalm,'Model','m',.true.)
!  print*,'complete!'
!  !---------------------end model computation-----------------------------------
!
!  call compute_image_energy(nodalg,nodalm,nodalnu,energy)
!  write(111,'(i4,1x,e14.6)')i_tstep,energy
!  print*,'Residual in model approximation:',norm(nodalg(1,:)-nodalm(1,:))/maxval(nodalg)
!
!  ! update initial variables
!  nodalm0=nodalm

  print*,'-----------------------------------------------------------------------'
enddo time_step ! i_tstep time stepping loop
close(111)
close(10)

!  write \nabla.(\nu^2\nabla m)
call save_specfem_output(gdof_elmt,interpfgll,nodalnu0,nodalm0)

if(solver_type.eq.builtin_solver .and.solver_diagscale)deallocate(ndscale)
if(mattype==0)deallocate(mat_id,mat_domain,gam,ym,coh,nu,phi,psi)
deallocate(g_coord,g_num,isnode)
deallocate(interpfgll,storekmat)
call free_ghost(ngpart)
!-----------------------------------

return
end subroutine semimage2d
!===========================================


