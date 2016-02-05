subroutine genESM(element, kele)

use mytypes
implicit none

type(elementtype) :: element
real, dimension(8,8) :: kele
real, dimension(4):: x=0, y=0 !the nodal x and y values
real :: E, nu ! material properties
!real, dimension(4,2) :: C=0 ! The matrix that holds the derivatives of the shape functions
real :: Jacobian=0 !the jacobain

real, dimension(3,3), target :: planestrain=0, planestress=0 ! the stress-strain relations for the respective cases
real, dimension(3,3), pointer :: ss(:,:)
integer :: i

!writing the nodal coordinates in the appropriate arrays for easier handling
do i=1, size(element%node)
	x(i)= element%node(i)%p%x
	y(i)= element%node(i)%p%y
end do

E=element%properties%E
nu= element%properties%nu


! defining the stress-strain relation matrices
planestrain(1,1) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
planestrain(1,2) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
planestrain(2,1) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
planestrain(2,2) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
planestrain(3,3) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu)

!todo: introduce the plane stress formulation


! todo: introduce a function that selects the right stress-strain relation matrix

ss=>planestrain
	

! all we need to do now is to call the assembling subroutine
kele=0

call gaussint(kele)





contains

function C(xi, eta) !The derivatives of the shape functions: C(:,1)=(d N_i)/(d xi) C(:,2)=(dN_i)/(d eta)
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


! The derivatives relative to the global coordinates are 
!
!				| dN_1/dx 	dN_1/dy |
!				| dN_2/dx 	dN_2/dy |
! [C][Jinv] = 	| dN_3/dx 	dN_3/dy |
! 				| dN_4/dx 	dN_4/dy |

function AA(xi, eta)
	real, dimension(3,8) :: AA
	real, dimension(4,2) :: globaldiff=0
	real :: xi, eta
	AA=0
	globaldiff= matmul(C(xi,eta),Jinv(xi,eta))

	do i=1,4
		AA(1,2*i-1)= globaldiff(i,1)
		AA(2, 2*i) = globaldiff(i,2)
		AA(3,2*i-1)= globaldiff(i,2)
		AA(2, 2*i) = globaldiff(i,1)
	end do
end function AA



! Inverse of the Jacoby matrix
! this is:
!
!		 |  dxi/dx   dxi/dy  |
!  Jinv= | deta/dx  deta/dy  |


function Jinv(xi,eta) !inverse of the jacobian matrix
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


! and the determinant of the Jacobian:

function djac(xi,eta)
	real, intent(in) :: xi,eta
	real :: djac

	djac = -xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2) * yy&
     &(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3&
     &) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1 &
     &+ xx(1) * eta * yy(3) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1&
     &) * yy(2) * xi / 0.8D1 - xx(1) * yy(3) * xi / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 - xx(2) * yy(1) * xi&
     & / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D&
     &1 + xx(3) * eta * yy(2) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 - xx(&
     &3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(4) * et&
     &a * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * x&
     &i / 0.8D1
end function djac




! we have the differential matrix operator:
!	 | d/dx   0   |
! A= |  0    d/dy |
!	 | d/dy  d/dx |
!
! to produce the differential equation in displacement formulation:
! [A]^T [B] [A] {e} = {f}
!
! if we substitute the continuous displacement {e} ={e_x, e_y}^T with the discrete one:
! {ê} = sum(N_i*{e}_i) or rather {ê}= [N]{e_i} then we have:
!
!		| dN_1/dx 		0 		dN_2/dx 		0 		dN_3/dx  		0 		dN_4/dx 		0	|
!  AA= 	| 			dN_1/dy 		0 		dN_2/dy			0 		dN_3/dy 	 	0 		dN_4/dy	|
!		| dN_1/dy 	dN_1/dx 	dN_2/dy  	dN_2/dx  	dN_3/dy 	dN_3/dx 	dN_4/dy 	dN_4/dx |
!
! and in local coordinates the equivalent of A is:
! 
!		| dxi/dx*d/dxi+ deta/dx*d/deta          0                      |
! A= 	|      0                         dxi/dy*d/dxi + deta/dy*d/deta |
!       | dxi/dy*d/dxi + deta/dy*d/deta  dxi/dx*d/dxi+ deta/dx*d/deta  |
!

! so, the entire element stiffness matrix prior to integration is:
! [k] = [AA]^T [SS] [AA]

function k(xi, eta)
 	real, dimension(8,8) ::k
 	real :: xi, eta

 	k=matmul(transpose(AA(xi,eta)),matmul(planestrain,AA(xi,eta)))*djac(xi,eta)
 	
 end function k


subroutine gaussint(dummymatrix)

	implicit none
	real, dimension(8,8), intent(inout) :: dummymatrix
	real, dimension(3) :: xi, eta
	real:: w1, w2, w
	integer::i,j
	real, dimension(8,8):: intermedres=0

	xi=(/sqrt(3./5.),-sqrt(3./5.),0./)
	eta=(/sqrt(3./5.),-sqrt(3./5.),0./)
	intermedres=0
! the three point gauss legendre quadrature
	do i =1,3
		if (i==3) then
			w1=8./9.
		else
			w1=5./9.
		end if

		do j=1,3
			if (j==3) then
				w2=8./9.
			else
				w2=5./9.
			end if
			
			w=w1*w2
			intermedres=intermedres+w*k(xi(i),eta(j))
		end do
	end do
	dummymatrix = intermedres
end subroutine gaussint


end subroutine genESM