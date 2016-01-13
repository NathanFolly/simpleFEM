
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

!Including the relevant petsc header files:
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscviewer.h>
!!!!!!!!!

!!!!!!!!\
! These Variables are created for use with PETSc

Vec 			u, b
Mat 			Apet
KSP 			ksppet
PC 				pcpet
PetscReal 		normpet, tolpet
PetscErrorCode 	ierrpet
PetscInt 		ipet, npet, colpet(3), itspet, i1pet, i2pet, i3pet
PetscBool		flgpet
PetscMPIInt		sizepet, rankpet
PetscScalar		nonepet, onepet, valuepet(3)

PetscViewer 	viewer
PetscReal 		petsckele(8,8), petscdummy ! We need these petsc type variables because 
				!the MatSetValues subroutine refuses to produce the correct output receiving 
				!anything else than petsc-variables

!!!!!!!


type(elementtype), allocatable, target :: mesh(:) ! essentially the mesh is just an array of elements
type(elementtype), pointer :: meshpointer(:)
real, dimension(8,8) :: kele=0 ! This is the dummy element stiffness matrix that
!gets updated by the subroutine generateesm
character*50, parameter :: meshfile='testmesh_2D_box_quad.msh'

integer :: ndof_nodal, ndof_global !The respective degrees of freedom 

integer :: quadstart=0, quadend=0, quadcounter=0, nnodes=0
integer :: i,j,k, ii, jj

integer, dimension(8) :: rowmap, columnmap !vectors mapping the local dof to the global dof

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


! Assembly of the global matrix is one directlyin PetSC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PETSC-PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Initialization

! PETSc initialize checks for an existing petsc database. since we are not using one in this case,
! set the first argument in the initialization subroutine to be PETSC_NULL_CHARACTER

call PetscInitialize(PETSC_NULL_CHARACTER,ierrpet)



! MPI_Comm_size: number of processes running in the specified communicator.
! in this case:
! PETSC_COMM_WORLD : all petsc MPI processes currently runnign that petsc is aware of
! hence sizepet: all the currently running processes petsc is aware of 
call MPI_Comm_size(PETSC_COMM_WORLD, sizepet, ierrpet)
! if (sizepet .ne. 1) then
! 	! MPI_Comm_rank determines the rank of the calling proces (this one) inside the specified
! 	! group of processes (here : all of them)
! 	call MPI_Comm_rank(PETSC_COMM_WORLD,rankpet, ierrpet)
! 	if (rankpet .eq. 0) then 
! 		write(*,*) 'This is a uniprocessor calculation only!'
! 	endif

! 	SETERRQ(PETSC_COMM_WORLD, 1, ' ', ierrpet)
! endif


! I believe those variables are not needed

nonepet = -1.0
onepet	= 1.0
npet 	= 10
i1pet 	= 1
i2pet 	= 2
i3pet	= 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! This subroutine checks whether we have input an option when calling the program
! e.g. : ./simpleFEM -n 100
! then the program sets the integer variable npet to 100 


!  call PetScOptionsGetInt(PETSC_NULL_CHARACTER, '-n',npet, flgpet, ierrpet)


! however, we know the dimension of the problem from reading the mesh and it is 
! the global degree of freedom





! Creating the matrix

call MatCreate(PETSC_COMM_WORLD, Apet, ierrpet)
call MatSetSizes(Apet, PETSC_DECIDE, PETSC_DECIDE, ndof_global, ndof_global, ierrpet)
call MatSetFromOptions(Apet, ierrpet)
call MatSetType(Apet,MATAIJ,ierrpet)
call MatSetUp(Apet,ierrpet)


! Assembling the matrix
!! CAREFUL! petsc matrices use indices from 0 N-1 while fortran uses indices from 1 to N

do k = quadstart, quadend
	call generateesm(kele,mesh(k),properties)
	petsckele=kele
	do i=1,8
		do j=1,8
		call MatsetValue(Apet, globaldof(mesh(k),i)-1, globaldof(mesh(k),j)-1,petsckele(i,j),ADD_VALUES, ierrpet)
		end do
	end do
	! do i=1,8
	! 	rowmap(i) = globaldof(mesh(k),i)
	! 	columnmap(i) = globaldof(mesh(k),i)
	! end do
	! call MatSetValues(Apet,8,rowmap-1,8,columnmap-1,petsckele,ADD_VALUES,ierrpet)	
end do

call MatAssemblyBegin(Apet, MAT_FINAL_ASSEMBLY, ierrpet)
call MatAssemblyEnd(Apet, MAT_FINAL_ASSEMBLY, ierrpet)

! Let's look at the matrix we assembled :
call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierrpet)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierrpet)
call PetscViewerFileSetName(viewer, 'Kay', ierrpet)
call MatView(Apet,viewer,ierrpet)


! creating the vectors : one is created from scratch and then duplicated

call VecCreate(PETSC_COMM_WORLD, u, ierrpet)
call VecSetSizes(u, PETSC_DECIDE, ndof_global, ierrpet)
call VecSetFromOptions(u, ierrpet)
call VecDuplicate(u, b, ierrpet)
! do we need this? : call VecDuplicate(xpet, upet, ierrpet)


! setting the force vector
call VecSet(b,0,ierrpet)

! creating the linear solver
call KSPCreate(PETSC_COMM_WORLD,ksppet,ierrpet)
! setting the operators, using the same matrix that defines the linear system as preconditioning matrix
call KSPSetOperators(ksppet,Apet,Apet,ierrpet)
! creating a pointer to the preconditioning context of our Krylov space solver: 
call KSPGetPC(ksppet, pcpet, ierrpet)
! ... so that we can easily manipulate it here:
call PCSetType(pcpet,PCJACOBI,ierrpet)
!setting the tolerances for the solver:
tolpet=1.d-7
call KSPSetTolerances(ksppet,tolpet, PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER, ierrpet)
! we can also manipulate the solver options when later calling the program with flags.
! For this to be possible, we need to call the following subroutine:
call KSPSetFromOptions(ksppet,ierrpet)

!!!! now we solve:!!!!!!
call KSPSolve(ksppet,b,u,ierrpet)



! viewing the solver info

call KSPView(ksppet, PETSC_VIEWER_STDOUT_WORLD, ierrpet)

!!!! check solution and cleanup

! wont work becaus we have no exact solution:  call VecAXPY(xpet,nonepet, upet, ierrpet)
! we could check for the residual with that :   call VecNorm(xpet,NORM_2, normpet, ierrpet)
call KSPGetIterationNumber(ksppet, itspet, ierrpet)
! if (normpet .gt. 1.e-12) then
! 	write(6,100) normpet, itspet
! else
! 	write(6,200) itspet
! endif

100 format('Norm of error = ', e11.4,'Iterations = ', i5)
200 format('Norm of error < 1.e-12, Iterations 0 ', i5)





! Freeing workspace:

!call VecDestroy(xpet, ierrpet)
call VecDestroy(u, ierrpet)
call VecDestroy(b, ierrpet)
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