
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
		real ::  ux, uy, uz ! the displacements
		real, allocatable :: stressvector(:) !stresstensor as vector, in the case of 2D: (sigx, sigy, tauxy)
		real, allocatable :: strainvector(:) !Straintensor as vector, in the case of 2D: (epsx, epsy, gamm=2*epsxy)
		real :: vmises
		!todo: make all these pointers 	
	end type statetype


	type nodetype
		integer :: num
		real :: x, y, z
		type(statetype) :: state
	 end type nodetype

 	type bcobject
		integer :: bcnature
		character*40:: boundaryname
		real, dimension(3) ::conditions
	end type bcobject


	type ndptrarr ! the only way to create an array of pointers
		type(nodetype), pointer :: p
	end type ndptrarr

	 type elementtype
	 	! Right now, the node instances are stored within each element
	 	! Pointers might be the better idea!!
	 	integer :: ndof_local
	 	integer :: num, kind
	 	character*40 :: name
	 	logical :: hasbc=.false.
	 	integer :: bcnature=0
	 	real, dimension(2,3) :: bc=0
	 	type(ndptrarr), allocatable :: node(:)
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
PetscScalar, allocatable:: disppet(:)

PetscViewer 	viewer
PetscReal 		petsckele(8,8), petscdummy ! We need these petsc type variables because 
				!the MatSetValues subroutine refuses to produce the correct output receiving 
				!anything else than petsc-variables
PetscReal 		petforce(8)

!!!!!!!


type(elementtype), allocatable, target :: element(:) ! essentially the mesh is just an array of elements
type(nodetype), allocatable, target :: node(:)
real, dimension(8,8) :: kele=0 ! This is the dummy element stiffness matrix that
!gets updated by the subroutine generateesm
real, allocatable :: fele(:) !element force vector
real, allocatable :: displacement(:) ! Vector to store displacements in for postprocessing after the petsc-vector is distroyed
real, allocatable :: stressvector(:)
integer :: vecsize !needed for allocating the displacement(:) array

character*50, parameter :: meshfile='testmesh_2D_box_quad.msh'
character*50, parameter :: ppfile='postprocessingfile.msh'
character*100 :: cmdmessage='' ! used when we pass a command to the system

integer :: ndof_local=8, ndof_nodal, ndof_global !The respective degrees of freedom  

integer :: quadstart=0, quadend=0, quadcounter=0, nnodes=0
integer :: i,j,k, ii, jj

integer, dimension(8) :: rowmap, columnmap !vectors mapping the local dof to the global dof


integer :: filestatus

interface	! need this interface so that we can pass an allocatable array to 
			!subroutine readmesh and allocate it there according to the number 
			!of elements which we find in the meshfile
	subroutine readmesh(element, node, meshfilename, quadstart, quadend, quadcounter, nnodes)
		use mytypes
		implicit none
		type(elementtype), allocatable, intent(inout):: element(:)
		type(nodetype), allocatable, intent(inout), target:: node(:)
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

! We make a copy of the original meshfile to write the data into for external postprocessing
! unfortunately, if we make system calls in the postprocessing section of 
! this file, they do not seem to work. The problem might be PETSc

call execute_command_line('cp '//meshfile//' '//ppfile, wait=.false.)


! We read in the mesh created by the external program gmsh

call readmesh(element,node, meshfile, quadstart, quadend, quadcounter, nnodes)

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

! creating the vectors : one is created from scratch and then duplicated

call VecCreate(PETSC_COMM_WORLD, u, ierrpet)
call VecSetSizes(u, PETSC_DECIDE, ndof_global, ierrpet)
call VecSetFromOptions(u, ierrpet)
call VecDuplicate(u, b, ierrpet)
! do we need this? : call VecDuplicate(xpet, upet, ierrpet)


! setting the force vector to equal zero
call VecSet(b,0,ierrpet)

! Assembling the matrix
!! CAREFUL! petsc matrices use indices from 0 N-1 while fortran uses indices from 1 to N

do k = 1, size(element)
	! natural boundary conditions (nodal forces) are written to the force vector
	if ((element(k)%bcnature==1).or.(element(k)%bcnature==3)) then
		call lumpstress(element(k), fele)
		petforce=fele
		do i=1,element(k)%ndof_local
			rowmap(i)=globaldof(element(k),i)-1
		end do
		call VecSetValues(b,size(fele),rowmap,petforce,ADD_VALUES,ierrpet)
		
	else if (element(k)%kind == 3) then
		call generateesm(kele,element(k),properties)
		petsckele=kele ! because the Matsetvalue subroutine only accepts PETSC-scalars / PETSC-vectors
		do i=1,element(k)%ndof_local
			do j=1,element(k)%ndof_local
			call MatsetValue(Apet, globaldof(element(k),i)-1, globaldof(element(k),j)-1,petsckele(i,j),ADD_VALUES, ierrpet)
			end do
		end do
	end if
	! do i=1,8
	! 	rowmap(i) = globaldof(element(k),i)
	! 	columnmap(i) = globaldof(element(k),i)
	! end do
	! call MatSetValues(Apet,8,rowmap-1,8,columnmap-1,petsckele,ADD_VALUES,ierrpet)	
end do

! The necessary calls to complete the assembley of Vectors an matrices
call VecAssemblyBegin(b, ierrpet)
call VecAssemblyEnd(b, ierrpet)
call MatAssemblyBegin(Apet, MAT_FINAL_ASSEMBLY, ierrpet)
call MatAssemblyEnd(Apet, MAT_FINAL_ASSEMBLY, ierrpet)

! Let's look at the matrix we assembled :
call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierrpet)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierrpet)
call PetscViewerFileSetName(viewer, 'Kay', ierrpet)
call MatView(Apet,viewer,ierrpet)




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
!call KSPGetIterationNumber(ksppet, itspet, ierrpet)
! if (normpet .gt. 1.e-12) then
! 	write(6,100) normpet, itspet
! else
! 	write(6,200) itspet
! endif

100 format('Norm of error = ', e11.4,'Iterations = ', i5)
200 format('Norm of error < 1.e-12, Iterations 0 ', i5)




! call VecView(u, PETSC_VIEWER_STDOUT_WORLD, ierrpet)
! copying the displacements to a vector that lives outside petsc:
call VecGetSize(u, vecsize,ierrpet)

allocate(disppet(vecsize))
!idxpet = (/(i,i=0,vecsize-1)/)
sizepet = vecsize
call VecGetValues(u,sizepet,(/(i,i=0,vecsize-1)/),disppet, ierrpet)

!then assigning the displacements to the relevant nodes
call assigndisplacements(disppet,element)



! Freeing workspace:

!call VecDestroy(xpet, ierrpet)
call VecDestroy(u, ierrpet)
call VecDestroy(b, ierrpet)
call MatDestroy(Apet, ierrpet)
call KSPDestroy(ksppet, ierrpet)
call PetscFinalize(ierrpet)

!!!!!!!!!!!!!!!!!!! End of the petsc part


!!!!!Begin the postprocessing /  calculating the stresses

!calling the subroutine for stress recovery:

do i=1,size(element)
	if (element(i)%kind==3) then
		call recoverstress(element(i),properties)
	end if
end do


!writing the  data to the meshfile

open(5,file=ppfile, position='append', status='old',iostat=filestatus)

if (filestatus .ne. 0) then
	print *,'Error reading the post processing file . Fortran error Code:', filestatus
	call exit()
end if

write(5,'(A)') '$NodeData'
write(5,*) 1 !write the number of string tags
write(5,'(A)') 'Von Mises Stress' !the string tag
write(5,*) 1 ! Wehave one tag of type real (the point in time)
write(5,*) 0.0
write(5,*) 3 ! Three integer tags
write(5,*) 0 !Time step

write(5,*) 1 ! 1-component scalar field
write(5,*) nnodes ! number of associated nodal values

do i =1,nnodes
	write(5,*) node(i)%num, node(i)%state%vmises
end do	
write(5,'(A)') '$EndNodeData'

close(5)



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

	globalnodenum=element%node(localnodenum)%p%num

	globaldof=2*globalnodenum-(2-nodedof)

end function globaldof

subroutine lumpstress(element, forcevector)
	! reads the natural BC on an element(stress) and lumps it to the nodes 
	real, allocatable, intent(inout) :: forcevector(:)
	type(elementtype), intent(in) :: element
	integer :: i
!	real :: nodalforcex, nodalforcey, nodalforcez
	if (allocated(forcevector).and.(size(forcevector).ne.element%ndof_local)) then 
		deallocate(forcevector)
		allocate(forcevector(element%ndof_local))
	else if (.not. allocated(forcevector)) then
		allocate(forcevector(element%ndof_local))

	end if
	forcevector =0
	! stress is uniform over an element -> take the stress, multiply by the appropriate node distance, divide by number of nodes
	if (ndof_nodal==2) then
		do i = 1, size(element%node)-1
			forcevector(2*i-1) = forcevector(2*i-1)+element%bc(1,1)*abs(element%node(i)%p%y-element%node(i+1)%p%y)/2.
			forcevector(2*i) = forcevector(2*i)+element%bc(1,2)*abs(element%node(i)%p%x-element%node(i+1)%p%x)/2.
			forcevector(2*(i+1)-1) = forcevector(2*(i+1)-1)+element%bc(1,1)*abs(element%node(i)%p%y-element%node(i+1)%p%y)/2.
			forcevector(2*(i+1)) = forcevector(2*(i+1))+element%bc(1,2)*abs(element%node(i)%p%x-element%node(i+1)%p%x)/2.
		end do
	end if
end subroutine lumpstress

subroutine assigndisplacements(displacements, element)
	real(kind=8), intent(in) :: displacements(:)
	type(elementtype) :: element(:)

	do i=1,size(element)
		do j=1,size(element(i)%node)
		element(i)%node(j)%p%state%ux=displacements(globaldof(element(i),2*j-1))
		element(i)%node(j)%p%state%uy=displacements(globaldof(element(i),2*j))
		end do
	end do
end subroutine assigndisplacements

! 
end program main