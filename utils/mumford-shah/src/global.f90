! this module contains global parameters and variables
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
!  precision parameters
module set_precision
implicit none
integer,parameter :: kreal=8 !selected_real_kind(15)
end module set_precision
!============================================

module conversion_constants
use set_precision
implicit none
! CGI to SI units 
real(kind=kreal),parameter :: M2KM=1e-3_kreal
real(kind=kreal),parameter :: CGI2SI_MOMENT=1e-7_kreal
! SI to CGI units 
real(kind=kreal),parameter :: KM2M=1e+3_kreal
real(kind=kreal),parameter :: SI2CGI_MOMENT=1e+7_kreal
end module conversion_constants 
!============================================

module ksp_constants
use set_precision
implicit none

! KSP solver parameters
integer,parameter :: KSP_MAXITER=3000
real(kind=kreal),parameter :: KSP_RTOL=1.0e-8_kreal
real(kind=kreal),parameter :: KSP_ATOL=1.0e-30_kreal
real(kind=kreal),parameter :: KSP_DTOL=1.0e30_kreal
end module ksp_constants
!=============================================

! global parameters/variables
module global
use set_precision
implicit none
character(len=20) :: ptail
logical,parameter :: off=.false., on=.true.
character(len=3) :: method
INTEGER,PARAMETER :: NDIM=2 !NUMBER OF DIMENSIONS. THIS IS A 3D VERSION
INTEGER,PARAMETER :: NNDOF=1 ! NUMBER OF DEGREES OF FREEDOM PER NODE - UX, UY, UZ
INTEGER,PARAMETER,DIMENSION(NNDOF) :: IDOFU=(/ 1 /)!ID FOR EACH DEGREE OF
!freedom (per node)

integer :: model_input !0: DEFAULT full mesh input, 1: specfem2d input
! number of gauss-lobatto-legendre points along the x, y, and z directions.
! we always take ngllx=nglly=ngllz for the simplicity
! ngll=ngllx*nglly
integer :: ngllx,nglly,ngllz,ngll
integer,parameter :: ng=4 ! number of gauss points for FEM

integer :: myrank,myid,nproc !myrank is indexed from 0, myid is indexed from 1
integer :: ngdof !Number of nodal degrees of freedom per processor = nndof*nnode 
integer :: neq !number of equations per processor = ngdof - degrees of freedom
!lost due to constraints
integer,allocatable :: l2gdof(:)!map from local dof (in processor) to global dof
!(in entire system). l2gdof contains the global (system-wide) indices for the dof
!in that processor.

integer :: nsparse
integer,allocatable :: kcol_sparse(:),krow_sparse(:),kgcol_sparse(:),kgrow_sparse(:)

integer :: nenod, nedof ! number of elemental nodes (ie. number of nodes per
!element = ngll), number of elemental degrees of freedom -> nedof=nndof*nenod

integer :: ngnode ! number of geometrical nodes. usually, for FEM ngnode=nenod
integer :: nnode,nelmt! total # of nodes, ie. gll pts, per processor, and number of elements per processor
integer,allocatable :: mat_id(:)
integer :: mattype ! 0: block (DEFAULT), 1: tomo file

integer,allocatable :: g_num(:,:),gdof(:,:),ggdof(:,:)!g_num: global node IDs
!for each element (per processor). gdof: matrix of nodal dof (per processor). 
!eg. gdof(:,1) gives all dof in 1st node of that processor.
!ggdof is the same as gdof, but for all processors (global global dof matrix).

integer :: nedofu ! number of elemental degrees of freedoms for displacement
integer :: nedofphi ! number of elemental degrees of freedoms for gravity

! acceleration due to gravity
real(kind=kreal),parameter :: agrav=9.81_kreal
real(kind=kreal),allocatable :: g_coord(:,:) ! global coordinates
integer :: nmat !number of material domains
integer,allocatable :: mat_domain(:)
real(kind=kreal),allocatable :: gam(:),rho(:),ym(:),nu(:),coh(:),phi(:),psi(:)
logical,allocatable :: water(:)
integer,parameter :: ELASTIC=1,VISCOELASTIC=11
integer,parameter :: nmaxwell=1
integer :: nmat_viscoelas
real(kind=kreal),allocatable :: muratio_blk(:,:),viscosity_blk(:,:)

logical :: allelastic,iseqload,iswater,istraction,phinu,isimage
! pseudostatic coefficients for earthquake loading eqkh=ah/g, eqkv=av/g
real(kind=kreal) :: eqkx,eqky,eqkz
! where ah and av are horizontal and vertical pseduostatic accelerations

integer,parameter :: nst=6 ! number of unique stress components
! sx,sy,sz,tauxy,tauyz,tauzx
character(len=250) :: file_head,inp_path,out_path,part_path,mat_path
! displacement BC, ghost, traction, and water surface files
character(len=250) :: confile,idfile                                             
character(len=250),dimension(3) :: coordfile
character(len=250) :: matfile,uxfile,uyfile,uzfile,gfile,trfile,wsfile
integer :: cg_maxiter,nl_maxiter,nexcav,ninc,nsrf,ntstep
real(kind=kreal) :: cg_tol,nl_tol,dtstep

real(kind=kreal) :: eps_img,gam_img,eta_img

real(kind=kreal) :: max_elmtsize

! order of viscoelastic algorithm
! 1: First order, 2: Second order Simo and Hughes, 21: Second order Zienckiewicz
integer,parameter :: VISCO_ORDER=21 

! store array
real(kind=kreal),allocatable :: storederiv(:,:,:,:),storejw(:,:)

! convolution
integer,allocatable :: nnzero_diagconv(:)
! solver types
! DEVELOPER OPTIONS
! diagonally scale equations and use CG solver without preconditioning
logical,parameter ::  solver_diagscale=.true.
!for builtin solver
integer,parameter :: smart_solver=0  ! select appropriate solver automatically
integer,parameter :: builtin_solver=1! select builtin conjugate gradient solver
integer,parameter :: petsc_solver=2  ! select PETSC solver
integer :: solver_type=2!petsc_solver !smart_solver
! save options
type savedata_options
  logical :: disp,stress,porep,psigma,maxtau,nsigma,scf,vmeps
end type savedata_options
type(savedata_options) :: savedata

! others
character(len=1),parameter :: CR=achar(13) ! carriage return to overwrite
!previous line

integer :: stdout=6

end module global
!============================================

