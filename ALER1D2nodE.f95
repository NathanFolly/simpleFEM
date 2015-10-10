program bla
! Solving the easiest of all FE static cases: one-dimensional, axially loaded elastic rod. composed of 2-noded elements
!
! |---------------- -> F=10kN
! |-------L=1m-----|
! round rod with D=2.5mm -> A=19.63 mm^2

! Disclaimer: this program uses LAPACK subroutines. during compilation, LAPACK should be included: gfortran -framework accelerate ALER1D2nodE.f95


implicit none


integer :: nelem, nnodes, i,j, info
real, parameter :: E=2.1E9, crossec=19.634375E-6, L=1., F=7E3
real :: initelemL
real, allocatable :: A(:,:)
real, allocatable :: b(:), u(:)
integer, allocatable :: ipiv(:)
real :: blubb(2,2)


! blubb(1,1)=11
! blubb(1,2)=12
! blubb(2,1)=21
! blubb(2,2)=22



! print *,blubb
print *,'Number of elements?'
read *,nelem

nnodes=nelem+1


allocate(A(nnodes,nnodes),u(nnodes),b(nnodes), ipiv(nnodes))
initelemL=L/nelem


print *, 'creating the matrix A'
A=0
A(1,1)=1

do i=2,nnodes
A(i,i)=1
A(i,i-1)=-1
end do

b(1)=0   !That's the BC of zero-displacement of the first node
do i=2,nnodes
b(i)=((F/crossec)/E +1)*initelemL
end do


print *, 'This is the matrix A'
do i=1,nnodes
	print *, A(i,:)
end do

!calling the LAPACK subroutine:

!
print *,'Solving the LE'
call SGESV(nnodes,1, A, nnodes, ipiv, b, nnodes, info) !k1,k2,nrhs,A,nnodes,ipiv,b,nnodes, info )
print *,info
print *,b(:)

end program bla

