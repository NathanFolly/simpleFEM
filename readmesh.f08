!This subroutine is part of the simpleFEM project
! Written by Maximilian Hartig

! It reads a mesh file created by the gmsh program and translates it 
!into a format usable by simpleFEM




subroutine readmesh(mesh, meshfilename, quadstart, quadend, quadcounter, nnodes)
use mytypes
implicit none
type(elementtype), allocatable, intent(inout):: mesh(:)
integer, intent(inout) :: nnodes
character*50, intent(in) :: meshfilename
integer, intent(inout) :: quadstart, quadend, quadcounter! this variable is used to count the number of quad type elements

type(nodetype), allocatable :: node(:)
integer:: filestatus, nelements, readstatus, eltype, ntags
integer, allocatable:: eltags(:)
integer:: i,j
integer, dimension(3) :: linehead
integer, allocatable :: dummy(:)
character*50 :: line, dchar


open(5,file=meshfilename, status='old',iostat=filestatus)
if (filestatus .ne. 0) then
	print *,'Error reading the meshfile. Fortran error Code:', filestatus
	call exit()
end if

rewind(5)

10 read(5,'(A)', end =11) line
 	if (index(line, '$Nodes').ne. 0) then !here the nodes start
		read(5,*,iostat=readstatus)  nnodes	! the next line states how many nodes there are
		allocate(node(nnodes))
		do i=1,nnodes	! Then the node data starts and can be read
			read(5,*)  node(i)%num, node(i)%x, node(i)%y, node(i)%z
		end do
		goto 12
	else
		goto 10
	endif

11 print *,'Reached the end of the file without being able to find the nodes'


12 read(5,'(A)', end =13) line
 	if (index(line, '$Elements').ne. 0) then !here the Elements start
		read(5,*,iostat=readstatus)  nelements	! the next line states how many elements there are
		allocate(mesh(nelements))
		do i=1,nelements	! Then the element data starts and can be read
			read(5,'(A)') line
			read(line,*)  linehead
			eltype=linehead(2)
			ntags=linehead(3)
			mesh(i)%num=linehead(1)
			mesh(i)%kind=linehead(2)
			!for now, we do not use the tags provided by gmash
			if (mesh(i)%kind==15) then !kind 15 is a one-node element
				allocate(dummy(3+ntags+1))
				read(line,*) dummy
				allocate(mesh(i)%node(1))
				mesh(i)%node(1)=node(dummy(4+ntags))
			else if (mesh(i)%kind==1) then !kind 1 is a 2-node line element
				allocate(dummy(3+ntags+2))
				read(line,*) dummy
				allocate(mesh(i)%node(2))
				do j=1,2
					mesh(i)%node(j)=node(dummy(3+ntags+j))
				end do
			else if (mesh(i)%kind==3) then 	! kind 3 is a quadrilateral 4-node element
				if (quadstart==0) then
					quadstart = i
				endif
				quadcounter=quadcounter+1
				allocate(dummy(3+ntags+4))
				read(line,*)  dummy
				allocate(mesh(i)%node(4))
				do j=1,4
					mesh(i)%node(j)=node(dummy(3+ntags+j))
				end do
			else 
				print *,'This element type is not accepted by simpleFEM for now'
			end if
			quadend=quadstart+quadcounter-1
			deallocate(dummy)
		end do
		goto 14
	else
		goto 12		
	endif

13	print*,'reached the end of the file without encountering the $Elements flag'

14 close(5)



end subroutine readmesh




