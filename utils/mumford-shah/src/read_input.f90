module input
contains
! this subroutine reads the input information from a structured ASCII text file
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO:
!   - prompt warning or error for unknown argument/s
subroutine read_input(ismpi,inp_fname,errcode,errtag,ispartmesh)
use global
use math_constants,only:zero,zerotol
use string_library
implicit none

integer :: i
character(len=*),intent(in) :: inp_fname
logical,intent(in) :: ismpi
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
logical,optional,intent(in) :: ispartmesh
character(len=250) :: line
character(len=800) ::tag
character(len=80) :: strval,token
character(len=1) :: tmp_char
character(len=80),dimension(50) :: args
character(len=20) :: ptail_inp
integer :: id,ind,ios,narg,slen

integer :: bc_stat,preinfo_stat,mesh_stat,material_stat,control_stat, &
save_stat,image_stat
integer :: mat_count,nwmat
integer :: i_elmt,ielmt,i_node,inode,i_mat,imat,tmp_nelmt,tmp_nnode !,mat_domain

character(len=20) :: format_str
character(len=250) :: fname
character(len=150) :: data_path

integer :: ipart,nproc_inp ! partition ID
integer :: iselastic,ismatpart,istat,ival,issave,nexcavid_all
integer,allocatable :: ivect(:)
real(kind=kreal) :: rval
real(kind=kreal),allocatable :: rvect(:)

errtag="ERROR: unknown!"
errcode=-1

ipart=myid-1 ! partition ID starts from 0

! reading main input information
if(myid==1)write(*,'(a)',advance='no')'reading main input file...'
preinfo_stat=-1
bc_stat=-1
mesh_stat=-1
material_stat=-1
control_stat=-1
image_stat=0
save_stat=0

iselastic=0
ismatpart=1
isimage=.false.

savedata%disp=.false.
savedata%stress=.false.
savedata%porep=.false.
savedata%psigma=.false.
savedata%maxtau=.false.
savedata%nsigma=.false.
savedata%scf=.false.

mat_count=0
mattype=0

! default value
method='sem'
if(ismpi.and.nproc.gt.1)then
  inp_path='./partition/'
  part_path='./partition/'
  mat_path='./partition/'
else
  inp_path='./input/'
  mat_path='./input/'
endif
out_path='./output/'

model_input=0

! default CG and NL parameters
cg_tol=zerotol; cg_maxiter=100
nl_tol=zerotol; nl_maxiter=100

ninc=1 ! number of load increments
ntstep=1 ! number of time steps
dtstep=zero ! time step interval

if(ismpi.and.nproc.gt.1)then
  ptail_inp=trim(ptail)
else
  ptail_inp=""
endif
! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "'//trim(inp_fname)//'" cannot be opened!'
  return
endif
do
  read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
  if (ios/=0)exit
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle

  ! check for line continuation                                                 
  tag=''                                                                        
  do                                                                            
    call last_char(line,tmp_char,ind)                                           
    if(tmp_char.eq.'&')line(ind:ind)=''                                         
    tag=trim(tag)//trim(line)                                                   
    if(tmp_char.ne.'&')exit                                                     
    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    if(ios /= 0)then                                                            
      write(errtag,'(a)')'ERROR: line continuation incomplete!'                 
      return                                                                    
    endif                                                                       
  enddo   
  call first_token(tag,token)
  ! read pre information
  if (trim(token)=='preinfo:')then
    !print*,preinfo_stat
    !stop
    if(preinfo_stat==1)then
      write(errtag,*)'ERROR: copy of line type preinfo: not permitted!'
      return
    endif
    preinfo_stat=-1;
    call split_string(tag,',',args,narg)
    if(ismpi.or. (present(ispartmesh).and.ispartmesh))then
      nproc_inp=get_integer('nproc',args,narg);
      if(present(ispartmesh).and.ispartmesh)nproc=nproc_inp
      if(nproc_inp/=nproc)then
        write(errtag,*)'ERROR: number of processors and images must be equal!'
        return
      endif
      !phead=get_string('phead',args,narg)
    endif
    call seek_integer('type',ival,args,narg,istat)
    if(istat==0)model_input=ival
    ngllx=get_integer('ngllx',args,narg);
    nglly=get_integer('nglly',args,narg);
    ngllz=get_integer('ngllz',args,narg);
    ngll=ngllx*nglly*ngllz ! total GLL points
    if(model_input.eq.1)then
      nenod=ngll
    else 
      nenod=get_integer('nenod',args,narg);
    endif
    ! number of elemental degrees of freedom
    ! nedof=nndof*nenod
    ! number of geometrical nodes
    ngnode=get_integer('ngnod',args,narg);
    call seek_string('inp_path',strval,args,narg)
    if (.not. isblank(strval))inp_path=trim(strval)
    slen=len_trim(inp_path)
    if(inp_path(slen:slen)/='/')inp_path=trim(inp_path)//'/'
    if(ismpi)then
      call seek_string('part_path',strval,args,narg)
      if (.not. isblank(strval))part_path=trim(strval)
      slen=len_trim(part_path)
      if(part_path(slen:slen)/='/')part_path=trim(part_path)//'/'
    endif
    call seek_string('out_path',strval,args,narg)
    if (.not. isblank(strval))out_path=trim(strval)
    slen=len_trim(out_path)
    if(out_path(slen:slen)/='/')out_path=trim(out_path)//'/'

    preinfo_stat=1
    cycle
  endif
  ! read mesh information
  if (trim(token)=='mesh:')then
    if(mesh_stat==1)then
      write(errtag,*)'ERROR: copy of line type mesh: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    coordfile(1)=get_string('xfile',args,narg)
    coordfile(2)=get_string('yfile',args,narg)
    if(ndim.gt.2)coordfile(3)=get_string('zfile',args,narg)
    confile=get_string('confile',args,narg)
    if(model_input.ne.1)then
      idfile=get_string('idfile',args,narg)
      if(ismpi.and.nproc.gt.1)gfile=get_string('gfile',args,narg)
    endif

    mesh_stat=1
    cycle
  endif

  ! read material list
  if (trim(token)=='material:')then
    if(material_stat==1)then
      write(errtag,*)'ERROR: copy of line type material: not permitted!'
      return
    endif
    material_stat=-1
    call split_string(tag,',',args,narg)
    call seek_integer('ispart',ival,args,narg,istat)
    if(istat==0)ismatpart=ival
    ! if not partitioned default location is ../input
    if(ismatpart==0)mat_path='./input/'
    if(model_input.eq.1)mat_path=trim(inp_path)
    call seek_string('matpath',strval,args,narg)
    if (.not. isblank(strval))mat_path=trim(strval)
    slen=len_trim(mat_path)
    if(mat_path(slen:slen)/='/')mat_path=trim(mat_path)//'/'
    matfile=get_string('matfile',args,narg)
    call seek_integer('type',ival,args,narg,istat)
    if(istat==0)mattype=ival

    material_stat=1
    cycle
  endif

  ! read image loading information
  if (trim(token)=='image:')then
    if(image_stat==1)then
      write(errtag,*)'ERROR: copy of line type image: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    eps_img=get_real('eps',args,narg)
    gam_img=get_real('gam',args,narg)
    eta_img=get_real('eta',args,narg)
    image_stat=1
    isimage=.true.
    cycle
  endif

  ! read control information
  if (trim(token)=='control:')then
    if(control_stat==1)then
      write(errtag,*)'ERROR: copy of line type control: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    !cg_tol=get_real('cg_tol',args,narg)
    !cg_maxiter=get_integer('cg_maxiter',args,narg)
    !nl_tol=get_real('nl_tol',args,narg)
    !nl_maxiter=get_integer('nl_maxiter',args,narg)
    call seek_real('cg_tol',rval,args,narg,istat)
    if(istat==0)cg_tol=rval
    call seek_integer('cg_maxiter',ival,args,narg,istat)
    if(istat==0)cg_maxiter=ival
    call seek_real('nl_tol',rval,args,narg,istat)
    if(istat==0)nl_tol=rval
    call seek_integer('nl_maxiter',ival,args,narg,istat)
    if(istat==0)nl_maxiter=ival
    ! time step variables
    call seek_real('dt',rval,args,narg,istat)
    if(istat==0)dtstep=rval
    call seek_integer('ntstep',ival,args,narg,istat)
    if(istat==0)ntstep=ival
    !---------------------------

    call seek_integer('ninc',ival,args,narg,istat)
    if(istat==0)ninc=ival
    control_stat=1
    cycle
  endif

  ! read save options
  if (trim(token)=='save:')then
    if(save_stat==1)then
      write(errtag,*)'ERROR: copy of line type save: not permitted!'
      return
    endif
    save_stat=-1
    call split_string(tag,',',args,narg)
    call seek_integer('disp',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%disp=.true.
    call seek_integer('stress',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%stress=.true.
    call seek_integer('porep',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%porep=.true.
    call seek_integer('psigma',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%psigma=.true.
    call seek_integer('nsigma',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%nsigma=.true.
    call seek_integer('scf',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%scf=.true.
    call seek_integer('vmeps',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%vmeps=.true.

    save_stat=1
    cycle
  endif

  write(errtag,'(a)')'ERROR: invalid line type: "'//trim(token)//'"!'
  return


enddo ! do

! check input status
if (preinfo_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read pre information! make sure the line with "preinfo:" token is correct.'
  return
endif

! check input status
if (mesh_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read mesh information! make sure the line with "mesh:" token is correct.'
  return
endif

! check material status
if (material_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read material information! make sure the line with "material:" token is correct.'
  return
endif

! check control status
if (control_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read control information! make sure the line with "control:" token is correct.'
  return
endif
if(myid==1)write(*,*)'complete!'
!--------------------------------------------------------

if(model_input.eq.1)then
  errcode=0
  return
endif

! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

if(myid==1)write(*,'(a)',advance='no')'reading mesh & material properties...'
!write(format_str,*)ceiling(log10(real(nproc)+1))
!format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
! read coordinates information
do i=1,ndim
  !open(unit=11,file=trim(inp_path)//trim(coordfile(i)),action='read',status='old')
  ! open output file
  !write(fname, fmt=format_str)trim(inp_path)//trim(coordfile(i))//'_proc',ipart
  fname=trim(data_path)//trim(coordfile(i))//trim(ptail_inp)
  !print*,fname
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  read(11,*)tmp_nnode
  if(i==1)then
    nnode=tmp_nnode
    allocate(g_coord(ndim,nnode))
  endif
  if(tmp_nnode/=nnode)then
    write(errtag,'(a)')'ERROR: total number of nodes mismatch!'
    return
  endif
  if(ismpi.and.nproc.gt.1)then
    do i_node=1,nnode
      read(11,*)inode,g_coord(i,inode)
    enddo
  else
    do i_node=1,nnode
      read(11,*)g_coord(i,i_node)
    enddo
  endif
enddo
close(11)
!sync all
! read connectivity
!open(unit=11,file=trim(inp_path)//trim(confile),action='read',status='old')
! open output file
!write(fname, fmt=format_str)trim(inp_path)//trim(confile)//'_proc',ipart
fname=trim(data_path)//trim(confile)//trim(ptail_inp)
!print*,fname
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)nelmt
allocate(g_num(nenod,nelmt))
!print*,nenod,nelmt,ismpi,fname
if(ismpi.and.nproc.gt.1)then
  do i=1,nelmt
    read(11,*)ielmt,g_num(:,ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)g_num(:,i)
  enddo
endif
close(11)


! read material id
!open(unit=11,file=trim(inp_path)//trim(idfile),action='read',status='old')
! open output file
!write(fname, fmt=format_str)trim(inp_path)//trim(idfile)//'_proc',ipart
fname=trim(data_path)//trim(idfile)//trim(ptail_inp)
!print*,fname
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)tmp_nelmt
if(tmp_nelmt/=nelmt)then
  write(errtag,'(a)')'ERROR: total number of elements mismatch!'
  !,fname,tmp_nelmt,nelmt
  return
endif
allocate(mat_id(nelmt))
if(ismpi.and.nproc.gt.1)then
  do i=1,nelmt
    read(11,*)ielmt,mat_id(ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)mat_id(i)
  enddo
endif
close(11)

! for partmesh library following information is read in the library itself      
if(mattype.eq.0)then
  if(.not.present(ispartmesh).or. .not.ispartmesh)then 
    ! read material lists
    ! open output file
    !write(fname, fmt=format_str)trim(inp_path)//trim(matfile)//'_proc',ipart
    if(ismatpart==0)then ! material file not partitioned
      fname=trim(mat_path)//trim(matfile)
    elseif(ismatpart==1)then ! material file partitioned
      fname=trim(mat_path)//trim(matfile)//trim(ptail_inp)
    else
      write(errtag,'(a)')'ERROR: illegal ismatpart value!'
      return
    endif
    !print*,fname
    open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
    if( ios /= 0 ) then
      write(errtag,'(a)')'ERROR1: file "'//trim(fname)//'" cannot be opened!'
      return
    endif
    read(11,*)
    read(11,*)nmat
    allocate(mat_domain(nmat),gam(nmat),rho(nmat),ym(nmat),coh(nmat),nu(nmat),phi(nmat),psi(nmat),  &
    water(nmat))
    do i=1,nmat
      read(11,*)imat,mat_domain(i),gam(i),ym(i),nu(i),phi(i),coh(i),psi(i)
    enddo
    nmat_viscoelas=count(mat_domain.eq.VISCOELASTIC)
    allocate(muratio_blk(nmaxwell,nmat_viscoelas),viscosity_blk(nmaxwell,nmat_viscoelas))
    ! read visoelastic properties
    do i=1,nmat_viscoelas
      read(11,*)imat,muratio_blk(:,i),viscosity_blk(:,i)
    enddo
    if(minval(mat_id)<1 .or. maxval(mat_id)>nmat)then
      write(errtag,'(a)')'ERROR: material IDs must be consistent with the defined material regions!'
      return
    endif
    rho=gam/9.81_kreal ! Kg/m3
  endif
endif
! read bc
errcode=0
if(myid==1)write(*,*)'complete!'

end subroutine read_input
end module input
