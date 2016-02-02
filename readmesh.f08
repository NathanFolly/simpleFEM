!This subroutine is part of the simpleFEM project
! Written by Maximilian Hartig

! It reads a mesh file created by the gmsh program and translates it 
!into a format usable by simpleFEM

! since we create each element here anyway, we will simply add the boundaryconditions as well



subroutine readmesh(element,node, meshfilename, quadstart, quadend, quadcounter, nnodes)
use mytypes
implicit none
type(elementtype), allocatable, intent(inout):: element(:)
type(nodetype), allocatable, intent(inout), target :: node(:)
integer, intent(inout) :: nnodes
character*50, intent(in) :: meshfilename
integer, intent(inout) :: quadstart, quadend, quadcounter! this variable is used to count the number of quad type elements

! for the boundary conditions
type(bcobject), allocatable :: bc(:)
character*20, parameter :: bcfile='bcfile'
integer :: nbc

type nametagtype
	integer :: objectdim
	integer :: associatedobject
	character*40 :: name
end type nametagtype

type(nametagtype), allocatable :: nametag(:)
character*40, allocatable :: taglist(:) 



integer:: filestatus, nelements, readstatus, eltype, ntags, nphysnames
integer, allocatable:: eltags(:)
integer:: i,j
integer, dimension(3) :: linehead
integer, allocatable :: dummy(:)
character*50 :: line, dchar

!!!!!! need an interface to pass an allocatable array:
interface
	subroutine readBC(bcfile, bc, nbc)
		use mytypes
		implicit none
		character*20, intent(in) :: bcfile
		type(bcobject), allocatable, intent(inout):: bc(:)
		integer, intent(inout) :: nbc
	end subroutine readBC
end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! actual code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(5,file=meshfilename, status='old',iostat=filestatus)
if (filestatus .ne. 0) then
	print *,'Error reading the meshfile. Fortran error Code:', filestatus
	call exit()
end if

rewind(5)


! Reading the nodes
10 read(5,'(A)', end =11) line
	if (index(line, '$PhysicalNames').ne. 0 ) then !check i fthere are an nametags specified
		read (5,*,iostat=readstatus) nphysnames
		allocate(nametag(nphysnames))
		allocate(taglist(nphysnames))
		do i=1,nphysnames
			read(5,*,iostat=readstatus) nametag(i)%objectdim, nametag(i)%associatedobject, nametag(i)%name
			taglist(nametag(i)%associatedobject)=nametag(i)%name
		end do
		goto 10
 	else if (index(line, '$Nodes').ne. 0) then !here the nodes start
		read(5,*,iostat=readstatus)  nnodes	! the next line states how many nodes there are
		allocate(node(nnodes))
		do i=1,nnodes	! Then the node data starts and can be read
			read(5,*)  node(i)%num, node(i)%x, node(i)%y, node(i)%z
			! todo: strain and stressvectors are currently for 2-D only. adjust this if coe is made 3-D
			allocate(node(i)%state%strainvector(3))
			allocate(node(i)%state%stressvector(3))
		end do
		goto 44
	else
		goto 10
	endif

11 print *,'Reached the end of the file without being able to find the nodes'

! before reading the elements, we need to know about the boundary conditions:
44 call readBC(bcfile, bc, nbc)



12 read(5,'(A)', end =13) line
 	if (index(line, '$Elements').ne. 0) then !here the Elements start
		read(5,*,iostat=readstatus)  nelements	! the next line states how many elements there are
		allocate(element(nelements))
		do i=1,nelements	! Then the element data starts and can be read
			read(5,'(A)') line
			read(line,*)  linehead
			eltype=linehead(2)
			ntags=linehead(3)
			element(i)%num=linehead(1)
			element(i)%kind=linehead(2)
			if (element(i)%kind==15) then !kind 15 is a one-node element
				allocate(dummy(3+ntags+1))
				read(line,*) dummy
				element(i)%name=taglist(dummy(4))
				allocate(element(i)%node(1))
				element(i)%node(1)%p=>node(dummy(4+ntags))
				element(i)%ndof_local=2*1
				! todo: insert name support for 1-node-elements
			else if (element(i)%kind==1) then !kind 1 is a 2-node line element
				allocate(dummy(3+ntags+2))
				read(line,*) dummy
				element(i)%name=taglist(dummy(4))
				allocate(element(i)%node(2))
				do j=1,2
					element(i)%node(j)%p=>node(dummy(3+ntags+j))
				end do
				element(i)%ndof_local=2*2
			else if (element(i)%kind==3) then 	! kind 3 is a quadrilateral 4-node element
				if (quadstart==0) then
					quadstart = i
				endif
				quadcounter=quadcounter+1
				allocate(dummy(3+ntags+4))
				read(line,*)  dummy
				element(i)%name=taglist(dummy(4))
				allocate(element(i)%node(4))
				do j=1,4
					element(i)%node(j)%p=>node(dummy(3+ntags+j))
				end do
				element(i)%ndof_local=2*size(element(i)%node)
			else 
				print *,'This element type is not accepted by simpleFEM for now'
			end if
			call appendbc(element(i))
			quadend=quadstart+quadcounter-1
			deallocate(dummy)
		end do
		goto 14
	else
		goto 12		
	endif

13	print*,'reached the end of the file without encountering the $Elements flag'

14 close(5)

print *, 'successfully read the mesh file'


contains 
subroutine appendbc(el)
	implicit none
	type(elementtype), intent(inout) :: el
	integer :: i
	do i = 1,nbc
		if (bc(i)%boundaryname==el%name) then
			el%hasbc = .true.
			el%bcnature = el%bcnature+bc(i)%bcnature
			el%bc(bc(i)%bcnature,:)=bc(i)%conditions(:)
		end if
	end do

end subroutine appendbc

end subroutine readmesh




