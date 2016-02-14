subroutine genESM(element, kele, bele)

use mytypes
implicit none

type(elementtype) :: element
real(kind=dp), dimension(8,8) :: kele
real(kind=dp), dimension(8) :: bele ! we have to integrate the nodal forces too
real(kind=dp), dimension(8) :: forcevec=0 ! The vector of the nodal forces
real(kind=dp), dimension(4):: x=0, y=0 !the nodal x and y values
real(kind=dp) :: E, nu ! material properties
!real, dimension(4,2) :: C=0 ! The matrix that holds the derivatives of the shape functions
!real(kind=dp) :: Jacobian=0 !the jacobain

real(kind=dp), dimension(3,3), target :: planestrain=0, planestress=0 ! the stress-strain relations for the respective cases
real(kind=dp), dimension(3,3), pointer :: ss(:,:)
integer :: i

!writing the nodal coordinates in the appropriate arrays for easier handling
do i=1, size(element%node)
	x(i)= element%node(i)%p%x
	y(i)= element%node(i)%p%y
end do

! filling up the forcevector with the nodal forces
do i=1, size(element%node)
	forcevec(2*i-1:2*i) = element%node(i)%p%force
end do


! allocating the 
! if (allocated(bele).and.(size(bele).ne.element%ndof_local)) then 
! 	deallocate(bele)
! 	allocate(bele(element%ndof_local))
! else if (.not. allocated(bele)) then
! 	allocate(bele(element%ndof_local))
! end if



E=element%properties%E
nu= element%properties%nu


! defining the stress-strain relation matrices
! todo: create stresss-strain relation matrix outside the souroutine and only reference it here
! 		defining them new everytime the subroutine is called will influence performance nagitively
!       
planestrain(1,1) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
planestrain(1,2) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
planestrain(2,1) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
planestrain(2,2) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
planestrain(3,3) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu)

!todo: introduce the plane stress formulation

planestress(1,1) = 1.D0
planestress(1,2) = nu
planestress(2,1) = nu
planestress(2,2) = 1.D0
planestress(3,3) = (1.D0-nu)/2.D0

planestress=planestress*(E)/(1.D0 - nu**2)
!planestress=planestress*1.D0/(1.D0 - nu**2)

! todo: introduce a function that selects the right stress-strain relation matrix

ss=>planestress	
	

! all we need to do now is to call the assembling subroutine
kele=0
bele=0

call gaussintstiffness(kele)
bele=forcevec
!call gaussintforce(bele)




contains

function shapefunc(xi,eta) ! the shape functions (N1,N1,N2,N2,N3,N3,N4,N4) needed to integrate the nodal forces 
	real(kind=dp), dimension (8) :: shapefunc
	real(kind=dp), intent(in) :: xi, eta

	shapefunc(1) = 1.D0/4.D0*(1.D0-xi)*(1.D0-eta)
	shapefunc(2) = 1.D0/4.D0*(1.D0-xi)*(1.D0-eta)
	shapefunc(3) = 1.D0/4.D0*(1.D0-xi)*(1.D0+eta)
	shapefunc(4) = 1.D0/4.D0*(1.D0-xi)*(1.D0+eta)
	shapefunc(5) = 1.D0/4.D0*(1.D0+xi)*(1.D0+eta)
	shapefunc(6) = 1.D0/4.D0*(1.D0+xi)*(1.D0+eta)
	shapefunc(7) = 1.D0/4.D0*(1.D0+xi)*(1.D0-eta)
	shapefunc(8) = 1.D0/4.D0*(1.D0+xi)*(1.D0-eta)

end function shapefunc

function C(xi, eta) !The derivatives of the shape functions: C(:,1)=(d N_i)/(d xi) C(:,2)=(dN_i)/(d eta)
	real(kind=dp), dimension (4,2) :: C
	real(kind=dp), intent(in):: xi, eta

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
	real(kind=dp), dimension(3,8) :: AA
	real(kind=dp), dimension(4,2) :: globaldiff=0
	real(kind=dp) :: xi, eta
	AA=0
	globaldiff= matmul(C(xi,eta),Jinv(xi,eta))

	do i=1,4
		AA(1,2*i-1)= globaldiff(i,1)
		AA(2, 2*i) = globaldiff(i,2)
		AA(3,2*i-1)= globaldiff(i,2)
		AA(3, 2*i) = globaldiff(i,1)
	end do
end function AA
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


! Inverse of the Jacoby matrix
! this is:
!
!		 |  dxi/dx   dxi/dy  |
!  Jinv= | deta/dx  deta/dy  |


function Jinv(xi,eta) !inverse of the jacobian matrix
	real(kind=dp), dimension(2,2) :: Jinv
	real(kind=dp), intent(in):: xi, eta

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
	real(kind=dp), intent(in) :: xi,eta
	real(kind=dp) :: djac

	djac = -x(1) * y(2) / 0.8D1 + x(1) * y(4) / 0.8D1 + x(2) * y&
     &(1) / 0.8D1 - x(2) * y(3) / 0.8D1 + x(3) * y(2) / 0.8D1 - x(3&
     &) * y(4) / 0.8D1 - x(4) * y(1) / 0.8D1 + x(4) * y(3) / 0.8D1 &
     &+ x(1) * eta * y(3) / 0.8D1 - x(1) * eta * y(4) / 0.8D1 + x(1&
     &) * y(2) * xi / 0.8D1 - x(1) * y(3) * xi / 0.8D1 - x(2) * eta &
     &* y(3) / 0.8D1 + x(2) * eta * y(4) / 0.8D1 - x(2) * y(1) * xi&
     & / 0.8D1 + x(2) * y(4) * xi / 0.8D1 - x(3) * eta * y(1) / 0.8D&
     &1 + x(3) * eta * y(2) / 0.8D1 + x(3) * y(1) * xi / 0.8D1 - x(&
     &3) * y(4) * xi / 0.8D1 + x(4) * eta * y(1) / 0.8D1 - x(4) * et&
     &a * y(2) / 0.8D1 - x(4) * y(2) * xi / 0.8D1 + x(4) * y(3) * x&
     &i / 0.8D1
end function djac





! and in local coordinates the equivalent of A is:
! 
!		| dxi/dx*d/dxi+ deta/dx*d/deta          0                      |
! A= 	|      0                         dxi/dy*d/dxi + deta/dy*d/deta |
!       | dxi/dy*d/dxi + deta/dy*d/deta  dxi/dx*d/dxi+ deta/dx*d/deta  |
!

! so, the entire element stiffness matrix prior to integration is:
! [k] = [AA]^T [SS] [AA]

function k(xi, eta)
 	real(kind=dp), dimension(8,8) ::k
 	real(kind=dp) :: xi, eta

 	k=matmul(transpose(AA(xi,eta)),matmul(ss,AA(xi,eta)))

 	
 end function k


subroutine gaussintstiffness(dummymatrix)

	implicit none
	real(kind=dp), dimension(8,8), intent(inout) :: dummymatrix
	real(kind=dp), dimension(3) :: xi, eta
	real(kind=dp):: w1, w2, w
	integer::i,j,l
	real(kind=dp), dimension(8,8):: intermedres=0.D1
	intermedres=0.
	xi=(/sqrt(3./5.),-sqrt(3./5.),0./)
	eta=(/sqrt(3./5.),-sqrt(3./5.),0./)
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
			intermedres=intermedres+w*djac(xi(i),eta(j))*k(xi(i),eta(j))
!			print *, djac(xi(i),eta(j))
		end do
	end do
	dummymatrix = intermedres
end subroutine gaussintstiffness

subroutine gaussintforce(dummyvector)
	real(kind=dp), dimension(8), intent(inout) :: dummyvector
	real(kind=dp), dimension(3) :: xi, eta
	real(kind=dp):: w1, w2, w
	integer::i,j,l
	real(kind=dp), dimension(8):: intermedres=0.D1
	intermedres=0.
	xi=(/sqrt(3./5.),-sqrt(3./5.),0./)
	eta=(/sqrt(3./5.),-sqrt(3./5.),0./)
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
			intermedres=intermedres+w*djac(xi(i),eta(j))*forcevec*shapefunc(xi(i),eta(j))
		end do
	end do
	dummyvector = intermedres
end subroutine gaussintforce


end subroutine genESM