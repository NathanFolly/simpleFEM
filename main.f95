
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



use mytypes
implicit none

!Tryout for petsc:
#include </home/hartig/petsc/include/petsc/finclude/petscsys.h>
#include </home/hartig/petsc/include/petsc/finclude/petscvec.h>
#include </home/hartig/petsc/include/petsc/finclude/petscmat.h>
#include </home/hartig/petsc/include/petsc/finclude/petscksp.h>
#include </home/hartig/petsc/include/petsc/finclude/petscpc.h>
!!!!!!!!!

!!!!!!!!\
! These Variables are created for use with PETSc

Vec 			xpet,bpet,upet
Mat 			Apet
KSP 			ksppet
PC 				pcpet
PetscReal 		normpet, tolpet
PetscErrorCode 	ierrpet
PetscInt 		ipet, npet, colpet(3), itspet, i1pet, i2pet, i3pet
PetscBool		flgpet
PetscMPIInt		sizepet, rankpet
PetscScalar		nonepet, onepet, valuepet(3)

!!!!!!!


type(elementtype), allocatable, target :: mesh(:) ! essentially the mesh is just an array of elements
type(elementtype), pointer :: meshpointer(:)
real, dimension(8,8) :: kele ! This is the dummy element stiffness matrix that
!gets updated by the subroutine generateesm
real, allocatable 	::	Kglobal(:,:) ! The global stiffness matrix
real, allocatable	::	u(:), b(:)
character*50, parameter :: meshfile='testmesh_2D_box_quad.msh'

integer :: ndof_nodal, ndof_global !The respective degrees of freedom 

integer :: quadstart=0, quadend=0, quadcounter=0, nnodes=0
integer :: i,j,k, ii, jj


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
 !type(elementtype)::testelement
 type(propertytype)::properties

 ! allocate(testelement%node(4))
 ! testelement%node(1)%x=0.
 ! testelement%node(1)%y=0.
 ! testelement%node(2)%x=0.
 ! testelement%node(2)%y=1.
 ! testelement%node(3)%x=1.
 ! testelement%node(3)%y=1.
 ! testelement%node(4)%x=1.
 ! testelement%node(4)%y=0.

 properties%E = 2E11
 properties%nu = 0.3

!call generateesm(kele,testelement,properties)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read in the mesh created by the external program gmsh
meshpointer=>mesh(:)



call readmesh(mesh, meshfile, quadstart, quadend, quadcounter, nnodes)

! We know what problem we want to solve (planar strain) so we know the number of degrees of freedom per node:
ndof_nodal= 2
!determinin ghte globval degree of freedom is easy:
ndof_global= ndof_nodal*nnodes
! we assemble the Matrix
! not using any assembly strategy. just adding up everything
print *,'done reading mesh'

allocate(Kglobal(ndof_global,ndof_global))
allocate(u(ndof_global))
allocate(b(ndof_global))
u=0
b=0


do k=quadstart, quadend !for all the quad elements do
	call generateesm(kele, mesh(k), properties) 
	do i=1,8
		ii = globaldof(mesh(k),i)
		do j=1,8
		jj = globaldof(mesh(k),j)
		Kglobal(ii,jj)=Kglobal(ii,jj)+kele(i,j)
		end do
	end do 
end do

open(11, file='Kprint')
write(11,*)  Kglobal(50,:)

! now, in theory, all we have to do is to solve the system:






!call generateesm(kele,testelement,properties)


!!!!!!!!!!!!!!!!!!!!!!  PETSC-PART !!!!!!!!!!!!!!!!!!!!!!!!


! Initialization

call PetscInitialize(PETSC_NULL_CHARACTER,ierrpet)
call MPI_Comm_size(PETSC_COMM_WORLD, sizepet, ierrpet)
if (sizepet .ne. 1) then
	call MPI_Comm_rank(PETSC_COMM_WORLD,rankpet, ierrpet)
	if (rankpet .eq. 0) then 
		write(6,*) 'This is a uniprocessor example only!'
	endif
	SETERRQ(PETSC_COMM_WORLD, 1, ' ', ierrpet)
endif
nonepet = -1.0
onepet	= 1.0
npet 	= 10
i1pet 	= 1
i2pet 	= 2
i3pet	= 3
call PetScOptionsGetInt(PETSC_NULL_CHARACTER, '-n',npet, flgpet, ierrpet)

! Creating the matrix

call MatCreate(PETSC_COMM_WORLD, Apet, ierrpet)
call MatSetSizes(Apet, PETSC_DECIDE, PETSC_DECIDE, npet, npet, ierrpet)
call MatSetFromOptions(Apet, ierrpet)
call MatSetUp(Apet,ierrpet)

! Assembling the matrix

valuepet(1) = -1.0
valuepet(2) = 2.0
valuepet(3) = -1.0

do 50 ipet=1,npet-2
	colpet(1) = ipet-1
	colpet(2) = ipet
	colpet(3) = ipet+1
	call MatSetValues(Apet,i1pet,ipet,i3pet,colpet,valuepet,INSERT_VALUES,ierrpet)
50	continue
	i= npet -1
	colpet(1) = npet -2
	colpet(2) = npet -1	
	call MatSetValues(Apet,i1pet,ipet,i3pet,colpet,valuepet,INSERT_VALUES,ierrpet)
	ipet = 0
	colpet(1) = 0
	colpet(2) = 1
	valuepet(1)	= 2.0
	valuepet(2)	= -1.0
	call MatSetValues(Apet,i1pet,ipet,i3pet,colpet,valuepet,INSERT_VALUES,ierrpet)
	call MatAssemblyBegin(Apet, MAT_FINAL_ASSEMBLY, ierrpet)
	call MatAssemblyEnd(Apet, MAT_FINAL_ASSEMBLY, ierrpet)


! creating the vectors : one is created from scratch and then duplicated

call VecCreate(PETSC_COMM_WORLD, xpet, ierrpet)
call VecSetSizes(xpet, PETSC_DECIDE, npet, ierrpet)
call VecSetFromOptions(xpet, ierrpet)
call VecDuplicate(xpet, bpet, ierrpet)
call VecDuplicate(xpet, upet, ierrpet)

! Setting the exact solution. then compute the right hand side vector

call VecSet(upet, onepet, ierrpet)
call MatMult(Apet, upet, bpet, ierrpet)

!!!!!!! creating the linear solver and settin various options
call KSPCreate(PETSC_COMM_WORLD,ksppet, ierrpet)

! Setting the operators.
!  HERE, the Matrix that defines the linear system also serves as preconditioning matrix 
call KSPSetOperators(ksppet,Apet,Apet,ierrpet)

!setting linear solver defaults for this problem (optional)
! Those options could also be configured at runtime

call KSPGetPC(ksppet, pcpet, ierrpet)
call PCSetType(pcpet, PCJACOBI, ierrpet)
tolpet = 1.d-7
call KSPSetTolerances(ksppet,tolpet,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierrpet)

! setting the runtime options

call KSPSetFromOptions(ksppet, ierrpet)

!!!!!!!!!!!!! solving the linear system
call KSPSolve(ksppet, bpet, xpet, ierrpet)

! viewing the solver info

call KSPView(ksppet, PETSC_VIEWER_STDOUT_WORLD, ierrpet)

!!!! check solution and cleanup

call VecAXPY(xpet,nonepet, upet, ierrpet)
call VecNorm(xpet,NORM_2, normpet, ierrpet)
call KSPGetIterationNumber(ksppet, itspet, ierrpet)
if (normpet .gt. 1.e-12) then
	write(6,100) normpet, itspet
else
	write(6,200) itspet
endif

100 format('Norm of error = ', e11.4,'Iterations = ', i5)
200 format('Norm of error < 1.e-12, Iterations 0 ', i5)

! Freeing workspace:

call VecDestroy(xpet, ierrpet)
call VecDestroy(upet, ierrpet)
call VecDestroy(bpet, ierrpet)
call MatDestroy(Apet, ierrpet)
call KSPDestroy(ksppet, ierrpet)
call PetscFinalize(ierrpet)

!!!!!!!!!!!!!!!!!!! End of the petsc part





contains

function globaldof(element, localdof) ! returns the global degree of freedom in depenande of the local degree of freedom
	integer :: globaldof
	type(elementtype), intent(in) :: element
	integer, intent(in) :: localdof
	integer:: localnodenum, nodedof, globalnodenum

	if (mod(localdof,2)==0) then 
		localnodenum= localdof/2
		nodedof=2
	else
		localnodenum=(localdof+1)/2
		nodedof=1
	endif

	globalnodenum=element%node(localnodenum)%num

	globaldof=2*globalnodenum-(2-nodedof)

end function globaldof


! 
end program main