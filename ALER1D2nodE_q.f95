program bla
! Solving a slightly more complex FE static case: one-dimensional, vertically loaded elastic rod. composed of 2-noded elements
!
!          q       
!   | | | | | | | |
!   v v v v v v v v
! |================ 
! |-------L=1m-----|
! round rod with D=2.5mm -> A=19.63 mm^2

! Disclaimer: this program uses LAPACK subroutines. during compilation, LAPACK should be included: gfortran -framework accelerate ALER1D2nodE.f95


implicit none


integer :: nelem, nnodes, i,j, info
real, parameter :: E=2.1E11, h=10E-3, w=10E-3, L=1., F=7E3, qq=10 !in SI units, qq in [N/m]
real :: initelemL
real, allocatable :: A(:,:), D(:,:), AA(:,:)
real, allocatable :: b(:), u(:), x(:), tfunctvec(:), derivedtfuncvec(:),q(:),y(:)
integer, allocatable :: ipiv(:)
real :: FTHM, M

type nodetype
	real :: u,w, momentum, Fy, Fx, x
	real :: belongstoelement(2)
end type nodetype

type elementtype
	type(nodetype) :: node(3)
	real :: length
end type elementtype

type(nodetype), allocatable :: node(:)
type(elementtype), allocatable :: element(:)


! Calculating the FTHM
FTHM=(w*h**3)/12
print *,FTHM


print *,'Number of elements?'
read *,nelem

nnodes=nelem+1



allocate(A(nnodes,nnodes), AA(nnodes,nnodes), D(nnodes,nnodes), u(nnodes), b(nnodes), ipiv(nnodes), tfunctvec(nnodes))
allocate(x(nnodes), node(nnodes), element(nnodes-2), q(nnodes),y(nnodes), derivedtfuncvec(nnodes))




! We#ll just assume  a constant force q:

do i=1,nnodes
	q(i)= qq/nelem
end do

! Defining the vector of discrete x:
initelemL=L/nelem

do i=1,nnodes
	x(i)=initelemL*(i-1)
	node(i)%x=x(i)
end do

! Assigning nodes to elements:

do i=1,size(element)
	element(i)%node(1)=node(i)
	element(i)%node(2)=node(i+1)
	element(i)%node(3)=node(i+2)
end do

! creating the matrix A consisting of the integral over the derivatives of the test functions:
A=0
A(1,1)=1 ! for the BC @x=0
do i=1,size(element)
	A(i+1,i+1) = E*FTHM*(1./(element(i)%node(2)%x-element(i)%node(1)%x)+1./(element(i)%node(3)%x-element(i)%node(2)%x))
	A(i+1,i) = -E*FTHM*(1./(element(i)%node(2)%x-element(i)%node(1)%x))
	A(i+1,i+2) = -E*FTHM*(1./(element(i)%node(3)%x-element(i)%node(2)%x))
end do
A(nnodes,nnodes)=1 ! for the BC @ x=l

AA=A ! We copy the matrix A because we'll need it a second time to compute the second case and the LAPACK solving subroutine will change A during the solving process



! maybe I'll need to set A(1,1) nad A(nnode,nnode) let's see

!creating matrix L and subsequently the vector b

D=0
do i=1,size(element)
	D(i+1,:)=(element(i)%node(3)%x-element(i)%node(1)%x)*1./2.
end do
D(:,1)=0
!oD(:,nnodes)=0

b=matmul(D,q)
b(1)=-qq/2./(E*FTHM)
!calling the LAPACK subroutine:
print *,b
!
print *,'Solving the LE'
call SGESV(nnodes,1, A, nnodes, ipiv, b, nnodes, info) !k1,k2,nrhs,A,nnodes,ipiv,b,nnodes, info )


y=integrate(integrate(b))


! Writing the data to files
open(unit=1, file='bendingline.dat')
do i=1,nnodes
	write(1,*) x(i), y(i)
end do
close(unit=1)

open(unit=2, file='bendingangle.dat')
do i=1,nnodes
	write(2,*) x(i), b(i)
end do 
close(unit=2)



!now let's change the boundary conditions:
b=matmul(D,q)
b(1)=0

A=AA !restore Matrix A from the copy we made since the original one has been altered during the solving procedure by the SGESV subroutine
call SGESV(nnodes,1, A, nnodes, ipiv, b, nnodes, info) !k1,k2,nrhs,A,nnodes,ipiv,b,nnodes, info )

! we need to integrate again:
y=integrate(integrate(b))

! aaaand let's write the output in a file
open(unit=1, file='bendinglineB.dat')
do i=1,nnodes
	write(1,*) x(i), y(i)
end do
close(unit=1)

open(unit=2, file='bendingangleB.dat')
do i=1,nnodes
	write(2,*) x(i), b(i)
end do 
close(unit=2)

!A text file that contains the commands for gnuplot
open(unit=10, file='gp.txt')
write(10,*) 'set term png'
write(10,*) 'f(x) = -1./',FTHM*E,'*',qq,'*(x**4/24.-1./6*x**3+1./4*x**2)'
write(10,*) 'g(x) = -1./',FTHM*E,'*',qq,'*(x**2/2.-x+1./2.)'
write(10,*) 'set output "bendinglines.png"'
write(10,*) 'set title "Comparing bending lines from the FEM solver and the analytical solution."'
write(10,*) 'plot [x=0:',L,'] f(x) title "analytical", "bendingline.dat" title "FEM"'
write(10,*) 'set output "bendingangles.png"'
write(10,*) 'set key right bottom'
write(10,*)	'set title "Comparing bending angles from the FEM solver and the analytical solution."'
write(10,*) 'plot [x=0:',L,'] g(x) title "analytical", "bendingangle.dat" title "FEM"'
write(10,*) 'g(x) = -1./',FTHM*E,'*',qq,'*(x**2/2.-x/2.)'
write(10,*) 'f(x) = -1./',FTHM*E,'*',qq,'*(x**4/24.-x**3/12+x/24)'
write(10,*) 'set output "bendinglines_fixedends.png"'
write(10,*) 'set title "Comparing bending lines from the FEM solver and the analytical solution."'
write(10,*) 'plot [x=0:',L,'] f(x) title "analytical", "bendinglineB.dat" title "FEM"'
write(10,*) 'set output "bendingangles_fixedends.png"'
write(10,*) 'set key right bottom'
write(10,*)	'set title "Comparing bending angles from the FEM solver and the analytical solution."'
write(10,*) 'plot [x=0:',L,'] g(x) title "analytical", "bendingangleB.dat" title "FEM"'
close(unit=10)

call SYSTEM('gnuplot gp.txt')

contains

function integrate(u)
	real :: u(:)
	real, allocatable :: integrate(:)
	real :: integrationlenth

	allocate(integrate(size(u)))
	integrationlenth = L/nelem
	integrate(1)=0
	do i=2,size(integrate)
		integrate(i)= integrate(i-1)+u(i-1)*integrationlenth
	end do
end function integrate
! function elementfunction(x,element)
! 	real :: x, xkm1, xk, xkp1
! 	real ::elementfunction
! 	type(element) :: element

! 	xk = element%nodes(2)%x
! 	xkm1= element%nodes(1)%x
! 	xkp1 = element%nodes(3)%x

! 	if (xkm1<x<xk) do
! 		elementfunction = (x-xkm1)/(xk-xkm1)
! 	else if (xk<x<xkp1) do
! 		elementfunction = (xkp1-x)/(xkp1-xk)
! 	else do
! 		elementfunction = 0
! 	end if

! end function elementfunction

! function derivedelementfunction(x,element)
! 	real :: x, xkm1, xk, xkp1
! 	real ::elementfunction
! 	type(element) :: element

! 	xk = element%nodes(2)%x
! 	xkm1= element%nodes(1)%x
! 	xkp1 = element%nodes(3)%x

! 	if (xkm1<x<xk) do
! 		elementfunction = (x)/(xk-xkm1)
! 	else if (xk<x<xkp1) do
! 		elementfunction = (-x)/(xkp1-xk)
! 	else do
! 		elementfunction = 0
! 	end if
! end function derivedelementfunction
end program bla

