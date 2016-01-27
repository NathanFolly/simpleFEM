! Subroutine that permits to recover the stress from the displacements that have 
! been calculated.
! Uses the more rudimentary approach (1) according to:
! http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf 

subroutine recoverstress(element, property)


use mytypes
implicit none

type(elementtype),intent(inout) :: element
type(propertytype),intent(in) :: property

real, allocatable:: ux(:), uy(:), x(:), y(:)
real, dimension(3):: stressvector=0, strainvector=0 ! this is for 2-D cases only

real :: eta=0, xi=0

real, dimension(2,2) :: Dhat=0
real, dimension(4,2) :: localcoordvec=0 ! the vector containing the local coordinates of the nodes (from 1 to 4)
real, dimension(2,4) :: u=0
real, dimension(3,3) :: A=0

real:: E=0, nu=0

integer :: i

!get the properties
E= property%E
nu = property%nu

!get the displacements
allocate(ux(size(element%node)))
allocate(uy(size(element%node)))


do i = 1, size(ux)
     ux(i)= element%node(i)%state%ux
     uy(i)= element%node(i)%state%uy
end do

u(:,1)=ux(:)
u(:,2)=uy(:)

!get the nodal coordinates:
allocate(x(size(element%node)))
allocate(y(size(element%node)))
do i=1,size(x)
     x(i)=element%node(i)%x
     y(i)=element%node(i)%y
end do


! creatin the vector of nodal points in local coords (xi, eta):
localcoordvec(1,:)= (/-1.,-1./) 
localcoordvec(2,:)= (/-1.,1./)
localcoordvec(3,:)= (/1.,1./)
localcoordvec(4,:)= (/1.,-1./)

! The strain-stress relation matrix:
      A(1,1) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      A(1,2) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      A(2,1) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      A(2,2) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      A(3,3) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu)



! for each node, evaluate the equation for the strain, then for the stress:
do i=1, size(x) !that is the number of nodes
	xi= localcoordvec(i,1)
	eta= localcoordvec(i,2)
	Dhat= matmul(u,matmul(C(xi,eta),Jinv(xi,eta)))
	
	strainvector(1)= Dhat(1,1)
	strainvector(2)= Dhat(2,2)
	strainvector(3)= Dhat(1,2)+Dhat(2,1)
	element%node(i)%state%strainvector= strainvector

	print *,strainvector

	stressvector= matmul(A,strainvector)
	element%node(i)%state%stressvector=stressvector
end do







! The derivatives of the shape functions
contains

function C(xi, eta)
	real, dimension (4,2) :: C
	real, intent(in):: xi, eta

 	C(1,1) = -0.1D1 / 0.4D1 + eta / 0.4D1 
    C(1,2) = -0.1D1 / 0.4D1 + xi / 0.4D1
    C(2,1) = -0.1D1 / 0.4D1 - eta / 0.4D1
    C(2,2) = 0.1D1 / 0.4D1 - xi / 0.4D1
    C(3,1) = 0.1D1 / 0.4D1 + eta / 0.4D1
    C(3,2) = 0.1D1 / 0.4D1 + xi / 0.4D1
    C(4,1) = 0.1D1 / 0.4D1 - eta / 0.4D1
    C(4,2) = -0.1D1 / 0.4D1 - xi / 0.4D1
end function C



! The inverse of the Jacobian:


function Jinv(xi,eta)
	real, dimension(2,2) :: Jinv
	real, intent(in):: xi, eta

	Jinv(1,1) = 2 * (y(1) * xi - y(2) * xi + y(3) * xi - y(4) * xi - y&
     &(1) + y(2) + y(3) - y(4)) / (x(1) * eta * y(3) - x(1) * eta * y(4)&
     & - x(2) * eta * y(3) + x(2) * eta * y(4) - x(3) * eta * y(1) + x(3&
     &) * eta * y(2) + x(4) * eta * y(1) - x(4) * eta * y(2) + x(1) * y(&
     &2) * xi - x(1) * y(3) * xi - x(2) * y(1) * xi + x(2) * y(4) * xi +&
     & x(3) * y(1) * xi - x(3) * y(4) * xi - x(4) * y(2) * xi + x(4) * y&
     &(3) * xi - x(1) * y(2) + x(1) * y(4) + x(2) * y(1) - x(2) * y(3) +&
     & x(3) * y(2) - x(3) * y(4) - x(4) * y(1) + x(4) * y(3))
      Jinv(1,2) = -2 * (x(1) * xi - x(2) * xi + x(3) * xi - x(4) * xi - &
     &x(1) + x(2) + x(3) - x(4)) / (x(1) * eta * y(3) - x(1) * eta * y(4&
     &) - x(2) * eta * y(3) + x(2) * eta * y(4) - x(3) * eta * y(1) + x(&
     &3) * eta * y(2) + x(4) * eta * y(1) - x(4) * eta * y(2) + x(1) * y&
     &(2) * xi - x(1) * y(3) * xi - x(2) * y(1) * xi + x(2) * y(4) * xi &
     &+ x(3) * y(1) * xi - x(3) * y(4) * xi - x(4) * y(2) * xi + x(4) * &
     &y(3) * xi - x(1) * y(2) + x(1) * y(4) + x(2) * y(1) - x(2) * y(3) &
     &+ x(3) * y(2) - x(3) * y(4) - x(4) * y(1) + x(4) * y(3))
      Jinv(2,1) = -2 * (y(1) * eta - y(2) * eta + y(3) * eta - y(4) * et&
     &a - y(1) - y(2) + y(3) + y(4)) / (x(1) * eta * y(3) - x(1) * eta *&
     & y(4) - x(2) * eta * y(3) + x(2) * eta * y(4) - x(3) * eta * y(1) &
     &+ x(3) * eta * y(2) + x(4) * eta * y(1) - x(4) * eta * y(2) + x(1)&
     & * y(2) * xi - x(1) * y(3) * xi - x(2) * y(1) * xi + x(2) * y(4) *&
     & xi + x(3) * y(1) * xi - x(3) * y(4) * xi - x(4) * y(2) * xi + x(4&
     &) * y(3) * xi - x(1) * y(2) + x(1) * y(4) + x(2) * y(1) - x(2) * y&
     &(3) + x(3) * y(2) - x(3) * y(4) - x(4) * y(1) + x(4) * y(3))
      Jinv(2,2) = 2 * (x(1) * eta - x(2) * eta + x(3) * eta - x(4) * eta&
     & - x(1) - x(2) + x(3) + x(4)) / (x(1) * eta * y(3) - x(1) * eta * &
     &y(4) - x(2) * eta * y(3) + x(2) * eta * y(4) - x(3) * eta * y(1) +&
     & x(3) * eta * y(2) + x(4) * eta * y(1) - x(4) * eta * y(2) + x(1) &
     &* y(2) * xi - x(1) * y(3) * xi - x(2) * y(1) * xi + x(2) * y(4) * &
     &xi + x(3) * y(1) * xi - x(3) * y(4) * xi - x(4) * y(2) * xi + x(4)&
     & * y(3) * xi - x(1) * y(2) + x(1) * y(4) + x(2) * y(1) - x(2) * y(&
     &3) + x(3) * y(2) - x(3) * y(4) - x(4) * y(1) + x(4) * y(3))
end function Jinv




end subroutine recoverstress


