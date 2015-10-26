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
real, parameter :: E=2.1E9, h=3E-3, w=5E-3, L=1., F=7E3
real :: initelemL
real, allocatable :: A(:,:), D(:,:)
real, allocatable :: b(:), u(:), x(:), tfunctvec(:), derivedtfuncvec(:),q(:)
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
FTHM=(w*h**4)/12



print *,'Number of elements?'
read *,nelem

nnodes=nelem+1



allocate(A(nnodes,nnodes), D(nnodes,nnodes), u(nnodes), b(nnodes), ipiv(nnodes), tfunctvec(nnodes), derivedtfuncvec(nnodes))
allocate(x(nnodes), node(nnodes), element(nnodes-2), q(nnodes))




! We#ll just assume  a constant force q:

do i=1,nnodes
	q(i)= 10
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
do i=1,size(element)
	A(i,i) = E*FTHM*(1/(element(i)%node(2)%x-element(i)%node(1)%x)-1/(element(i)%node(3)%x-element(i)%node(2)%x))
end do

! maybe I'll need to set A(1,1) nad A(nnode,nnode) let's see

!creating matrix L and subsequently the vector b

D=0
do i=2,size(element)
	D(i,i)=(element(i)%node(3)%x-element(i)%node(1)%x)/element(i)%node(2)%x*1/3
	print *, D(i,:)
end do

b=matmul(D,q)
print *,'x=',x
print *,'q=', q
print *, 'D=', D
print *,'vector b before solving',b

!calling the LAPACK subroutine:

!
print *,'Solving the LE'
call SGESV(nnodes,1, A, nnodes, ipiv, b, nnodes, info) !k1,k2,nrhs,A,nnodes,ipiv,b,nnodes, info )
print *,info
print *,b(:)

open(unit=1, file='data.dat')
do i=1,nnodes
	write(1,*) x(i), b(i)
end do
close(unit=1)

contains

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

