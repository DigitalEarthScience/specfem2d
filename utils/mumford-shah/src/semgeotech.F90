! this is a main program SPECFEM3D_GEOTECH
! REVISION:
!   HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
program semgeotech
! import necessary libraries
use global
use string_library, only : parse_file
!use math_constants
use input
use mesh
use mesh_spec
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
use visual

implicit none
integer :: funit,i,ios,istat,j,k
integer :: i_elmt,i_node,i_inc,i_srf,ielmt,igdof,imat,inode

integer :: gnod(4),gnum_quad(4),map2exodus(4),ngllxy,node_quad4(4)

character(len=250) :: arg1,inp_fname,out_fname,prog
character(len=150) :: path
character(len=20), parameter :: wild_char='********************'
character(len=20) :: ensight_etype
character(len=80) :: buffer,destag ! this must be 80 characters long
character(len=20) :: ext,format_str
character(len=250) :: case_file,geo_file,sum_file
integer :: npart,nt,tinc,tstart,twidth,ts ! ts: time set for ensight gold

real(kind=kreal) :: cpu_tstart,cpu_tend,telap,step_telap,max_telap,mean_telap

logical :: ismpi !.true. : MPI, .false. : serial
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode

character(len=250) :: errtag ! error message
integer :: errcode
logical :: isopen ! flag to check whether the file is opened

myid=1; nproc=1;
errtag=""; errcode=-1

call start_process(ismpi,myid,nproc,stdout)
!ipart=myid-1 ! partition id starts from 0

call get_command_argument(0, prog)
!----input and initialisation----
if (command_argument_count() <= 0) then
  call error_stop('ERROR: no input file!',stdout,myid)
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  if(myid==1)then
    write(stdout,'(a)')'Usage: '//trim(prog)//' [Options] [input_file]'
    write(stdout,'(a)')'Options:'
    write(stdout,'(a)')'    --help        : Display this information.'
    write(stdout,'(a)')'    --version     : Display version information.'
  endif
  !call sync_process
  call close_process()
elseif(trim(arg1)==('--version'))then
  if(myid==1)then
    write(stdout,'(a)')'SPECFEM3D_GEOTECH 1.1 Beta'
    write(stdout,'(a)')'This is free software; see the source for copying '
    write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
    write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  endif
  !call sync_process
  call close_process()
endif

! starting timer
call cpu_time(cpu_tstart)

! get processor tag
ptail=proc_tag(myid,nproc)

! get input file name
call get_command_argument(1, inp_fname)

! read input data
call read_input(ismpi,inp_fname,errcode,errtag)
if(errcode/=0)call error_stop(errtag,stdout,myid)
!call sync_process()
if(model_input.eq.1)call read_mesh(ismpi,inp_fname,errcode,errtag)
if(errcode/=0)call error_stop(errtag,stdout,myid)
tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myid==1)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',max_nelmt,' min:',min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',max_nnode,' min:',min_nnode
endif

if (trim(method)/='sem')then
  write(errtag,'(a)')'ERROR: wrong input for sem3d!'
  call error_stop(errtag,stdout,myid)
endif

call parse_file(inp_fname,path,file_head,ext)

ensight_etype='quad4'
ts=1 ! time set
tstart=1; tinc=1
!if(nexcav==0)then
!  nt=nsrf
!else
!  nt=nexcav+1 ! include 0 excavation stage (i.e., initial)
!  tstart=0
!endif
nt=ntstep
twidth=ceiling(log10(real(nt)+1.))

call compute_max_elementsize()

! create spectral elements
if(model_input.ne.1)then
  if(myid==1)write(stdout,'(a)',advance='no')'creating spectral elements...'
  call quad2spec(ndim,ngnode,nelmt,nnode,ngllx,nglly,errcode,errtag)
  if(errcode/=0)call error_stop(errtag,stdout,myid)
  if(myid==1)write(stdout,*)'complete!'
endif

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myid==1)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',max_nelmt,' min:',min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',max_nnode,' min:',min_nnode
endif

!call sync_process

!stop
nenod=ngll !(ngllx*nglly*ngllz) ! number of elemental nodes (nodes per element)
! number of elemental degrees of freedom
nedof=nndof*nenod

ngllxy=ngllx*nglly

! geometrical nodes (corner nodes) in EXODUS/CUBIT order
! bottom nodes
gnod(1)=1;
gnod(2)=ngllx
gnod(3)=ngllxy;
gnod(4)=gnod(3)-ngllx+1

! map sequential node numbering to exodus/cubit order for 4-noded quadrilateral
map2exodus=(/ 1,2,4,3 /)

case_file=trim(out_path)//trim(file_head)//trim(ptail)//'.case'
open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  call error_stop(errtag,stdout,myid)
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
write(11,'(a,a/)')'model:    ',trim(file_head)//trim(ptail)//'.geo'

write(11,'(a)')'VARIABLE'
write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','background',' ',  &
trim(file_head)//trim(ptail)//'.g'
!write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','blurred',' ',  &
!trim(file_head)//trim(ptail)//'.b'
!write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','model',' ',  &
!trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.m'

if(savedata%disp)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','nu',' ',  &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.nu'
endif
write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','dnu',' ',  &
trim(file_head)//trim(ptail)//'.dnu'
write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','dm',' ',  &
trim(file_head)//trim(ptail)//'.dm'
write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','wdm',' ',  &
trim(file_head)//trim(ptail)//'.wdm'
write(11,'(a)')'TIME'
write(11,'(a,i10)')'time set:',ts
write(11,'(a,i10)')'number of steps:',nt
write(11,'(a,i10)')'filename start number:',tstart
write(11,'(a,i10)')'filename increment:',tinc
write(11,'(a)',advance='no')'time values: '

do i=1,nt
  write(11,'(e12.5)',advance='yes')real(i)*dtstep
enddo
close(11)

! Format string
write(format_str,*)twidth
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//',a)'

! write geo file for inital stage (original)
! open Ensight Gold geo file to store mesh data
write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//trim(ptail)//'.geo'
npart=1
destag='unstructured meshes'
call write_ensight_geocoord(geo_file,destag,npart,nnode,real(g_coord),funit)

! writes element information
buffer=ensight_etype
write(funit)buffer
write(funit)nelmt*(ngllx-1)*(nglly-1)

! do not substract 1 for ensight file
do i_elmt=1,nelmt
    do j=1,nglly-1
      do i=1,ngllx-1
        ! corner nodes in a sequential numbering
        node_quad4(1)=(j-1)*ngllx+i
        node_quad4(2)=node_quad4(1)+1

        node_quad4(3)=node_quad4(1)+ngllx
        node_quad4(4)=node_quad4(3)+1

        ! map to exodus/cubit numbering and write
        !write(funit)g_num(node_quad4(map2exodus),i_elmt)
        gnum_quad=g_num(node_quad4(map2exodus),i_elmt)
        write(funit)gnum_quad !g_num(node_quad4(map2exodus),i_elmt)
      enddo
    enddo
enddo
close(funit)

! open summary file
sum_file = trim(out_path)//trim(file_head)//'_summary'//trim(ptail)
open(unit=10,file=trim(sum_file),status='replace',action='write',iostat=ios)
write(10,*)'--------------------------------------------'
write(10,*)'Result summary produced by SPECFEM3D_GEOTECH'
write(10,*)'--------------------------------------------'
close(10)

if(myid==1)write(stdout,'(a)')'--------------------------------------------'

! solver type
if(solver_type.eq.smart_solver)then
  if(ismpi)then
    solver_type=petsc_solver
  else
    solver_type=builtin_solver
  endif
endif
if(ismpi)then
  if(solver_type.eq.builtin_solver)then
    if(myid.eq.1)write(stdout,'(a)')'solver type: builtin parallel solver'
  elseif(solver_type.eq.petsc_solver)then
    if(nproc.gt.1)then
      if(myid.eq.1)write(stdout,'(a)')'solver type: PETSc parallel solver'
    else
      if(myid.eq.1)write(stdout,'(a)')'solver type: PETSc uniprocess solver'
    endif
  else
    write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
    call error_stop(errtag,stdout,myid)
  endif
else
  if(solver_type.eq.builtin_solver)then
    write(stdout,'(a)')'solver type: builtin serial solver'
  elseif(solver_type.eq.petsc_solver)then
    write(stdout,'(a)')'solver type: PETSc parallel solver'
    !write(errtag,'(a,i4)')'ERROR: currently PETSC solver cannot be used for serial version!'
    !call error_stop(errtag,stdout,myid)
  else
    write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
    call error_stop(errtag,stdout,myid)
  endif
endif   
! call main routines
call semimage2d(ismpi,gnod,sum_file,format_str)
!-----------------------------------

! compute elapsed time
call cpu_time(cpu_tend)
telap=cpu_tend-cpu_tstart
max_telap=maxscal(telap)
mean_telap=sumscal(telap)/real(nproc,kreal)

write(format_str,*)ceiling(log10(real(max_telap)+1.))+5 ! 1 . and 4 decimals
format_str='(3(f'//trim(adjustl(format_str))//'.4,1X))'
open(10,file=trim(sum_file),status='old',position='append',action='write')
write(10,*)'ELAPSED TIME, MAX ELAPSED TIME, MEAN ELAPSED TIME'
write(10,fmt=format_str)telap,max_telap,mean_telap
close(10)
!-----------------------------------

if(myid==1)then
  write(stdout,*) ! write new line
  write(stdout,'(a)')'--------------------------------------------'
  inquire(stdout,opened=isopen)
  if(isopen)close(stdout)
endif

call sync_process
call close_process()
contains

! routine must be called before converting hex8 to spectral elements
subroutine compute_max_elementsize 
use math_library,only:distance
implicit none
integer :: i_elmt,num(4)
real(kind=kreal),dimension(ndim) :: x1,x2,x3,x4
real(kind=kreal) :: d1,d2
real(kind=kreal) :: maxall_elmtsize,maxdiag
! find largest size of the element
max_elmtsize=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  x1=g_coord(:,num(1))
  x2=g_coord(:,num(2))
  x3=g_coord(:,num(3))
  x4=g_coord(:,num(4))
  d1=distance(x1,x3,2)
  d2=distance(x2,x4,2)
  maxdiag=max(d1,d2)
  if(maxdiag.gt.max_elmtsize)max_elmtsize=maxdiag
enddo
maxall_elmtsize=maxscal(max_elmtsize)
if(myid.eq.1)print*,'maximum element size across the diagonal:',maxall_elmtsize
end subroutine compute_max_elementsize 

end program semgeotech
!===========================================

