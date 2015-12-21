
! This is the main Program of simpleFEM written by Maximilian Hartig

module mytypes
	! we want some properties defined. Since we do assume constant global properties, we will later on create only one instance of this type.
	! IN later programs, we could also attach one property tape to each nodetype and make the properties variables.

	type propertytype
		real :: E, nu
	end type propertytype

	! We need to define a format for our mash data
	! What we want to do is to iterate through all elements and get the nodes that are involved
	! So and element centered data structure seems to make sense

	type statetype
		real :: sigx, sigy, tau, epsx, epsy, epsxy
	end type statetype


	type nodetype
		integer :: num
		real :: x, y, z
		type(statetype) :: state
	 end type nodetype

	 type elementtype
	 	! Right now, the node instances are stored within each element
	 	! Pointers might be the better idea!!
	 	integer :: num, kind
	 	character*20 :: name
	 	type(nodetype), allocatable :: node(:)
	end type elementtype


end module mytypes


program main

.#include '/opt/local/lib/petsc/include/petsc.h90'

use mytypes
implicit none

type(elementtype), allocatable, target :: mesh(:) ! essentially the mesh is just an array of elements
type(elementtype), pointer :: meshpointer(:)
real, dimension(4,8) :: kele ! This is the dummy element stiffness matrix that
!gets updated by the subroutine generateesm
real, allocatable 	::	Kglobal(:,:) ! The global stiffness matrix
real, allocatable	::	u(:), b(:)
character*50, parameter :: meshfile='testmesh_2D_box_quad.msh'

integer :: quadstart=0, quadend=0, quadcounter=0, nnodes=0
integer :: i,j,k, iglobal, jglobal


interface	! need this interface so that we can pass an allocatable array to 
			!subroutine readmesh and allocate it there according to the number 
			!of elements which we find in the meshfile
	subroutine readmesh(mesh, meshfilename, quadstart, quadend, quadcounter, nnodes)
		use mytypes
		implicit none
		type(elementtype), allocatable, intent(inout):: mesh(:)
		integer, intent(inout) :: nnodes
		character*50, intent(in) :: meshfilename
		integer, intent(inout) :: quadstart, quadend, quadcounter
	end subroutine readmesh
end interface


!!!!!!!!!!!! Can be deleted !!!!!!!!!!!!!
!in order to test the program let's create some nodes that make sense
! type(elementtype)::testelement
 type(propertytype)::properties

! allocate(testelement%node(4))
! testelement%node(1)%x=0
! testelement%node(1)%y=0
! testelement%node(2)%x=0
! testelement%node(2)%y=1
! testelement%node(3)%x=1
! testelement%node(3)%y=1
! testelement%node(4)%x=1
! testelement%node(4)%y=0

 properties%E = 2E11
 properties%nu = 0.3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read in the mesh created by the external program gmsh
meshpointer=>mesh(:)



call readmesh(mesh, meshfile, quadstart, quadend, quadcounter, nnodes)
! we assemble the Matrix
! not using any assembly strategy. just adding up everything
print *,'done reading mesh'

allocate(Kglobal(4*quadcounter,2*nnodes))
allocate(u(2*nnodes))
allocate(b(4*quadcounter))
u=0
b=0


do i=quadstart, quadend 
	call generateesm(kele, mesh(i), properties)
	do j=1,4 ! ideally vary leftmost index of the bigger Matrix more quickly because fortran is column-major -> that's quicker
		iglobal = 4*(i-quadstart)+j
		do k=1,4
			jglobal=2*mesh(i)%node(k)%num-1 ! for the x-displacement of node j
			Kglobal(iglobal,jglobal)= Kglobal(iglobal,jglobal)+kele(j,2*k-1)
			jglobal=2*mesh(i)%node(k)%num ! for the y- displacement of node j
			Kglobal(iglobal,jglobal)=Kglobal(iglobal,jglobal)+ kele(j,2*k)
		end do
	end do
end do

open(11, file='Kprint')
write(11,*)  Kglobal(50,:)

! now, in theory, all we have to do is to solve the system:






!call generateesm(kele,testelement,properties)




end program main