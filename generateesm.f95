! part of the simpleFEM Program by Maximilian Hartig
! This is essentially a subroutine that generates the element stiffness matrix for a quadrilateral element
! Input: list  of Nodes involved

! We create the stiffness matrix for a four-node quadrilateral element
! The matrix is a 8x4 matrix. After matrix multiplication with the displacement vector, it should equal the force (or zero)
! k x u = f
! Where k=8x4, u=8x1={ux1,uy1,ux2,uy2,ux3,uy3,ux4,uy4} and f=8x1={fx1,fy1,...,fy4}
!

subroutine generateesm(kk, element, property)

use mytypes


!This interface helps us to create a function pointer so as to select which entry of the element stiffness matrix to integrate in the subroutine gaussint 
abstract interface
	function kxx(xx,yy,xi,eta,E,nu)
		real:: kxx
		real, dimension(4) :: xx,yy
		real :: xi, eta, E, nu
	end function kxx
end interface


type(elementtype), intent(in) :: element
real, dimension(8,8), intent(inout)::kk
real, dimension(8,8) :: kkdummy ! I had to create this dummy because I made an error with the orientation of the element stiffness matrix and ended up calculating the transpose of k in Maple
type(propertytype), intent(in):: property

type(nodetype), dimension(4) :: node
real, dimension(8) :: localcoords
real, dimension(2,2) :: Jacoby
real, dimension(4) :: xx, yy
real::placeholder=0


! This pointer helps us select which function to integrate:

procedure(kxx), pointer :: kpntr => null()

! Testing
! let's see what k looks like:
print *,k(xx,yy,1.,1.,property%E,property%nu)




! !We create the vectors of nodal coordinates that are important for the evaluation of our integrals
! do i=1,4
! 	xx(i)=element%node(i)%x
! 	yy(i)=element%node(i)%y
! end do


! do i=1,8
! 	do j=1,4
! 		placeholder=0
! 		call gaussint(placeholder,i,j)
! 		kdummy(i,j) = placeholder
! 	end do
! end do

! ! We need to transpose the dummy element stiffness matrix to orient it the right way.
! ! This is necessary due to my mistake in creating the lengthy terms below
! k=transpose(kdummy)


contains




! subroutine gaussint(placeholder,ipos,jpos)

! 	implicit none
! 	real:: placeholder
! 	integer, intent(in):: ipos,jpos
! 	real, dimension(3) :: xi, eta
! 	real:: w1, w2, w
! 	integer::i,j
! 	real:: intermedres=0

! 	call funcselector(ipos,jpos)

! 	xi=(/sqrt(3./5.),-sqrt(3./5.),0./)
! 	eta=(/sqrt(3./5.),-sqrt(3./5.),0./)
! 	intermedres=0

! 	do i =1,3
! 		if (i==3) then
! 			w1=8./9.
! 		else
! 			w1=5./9.
! 		end if

! 		do j=1,3
! 			if (j==3) then
! 				w2=8./9.
! 			else
! 				w2=5./9.
! 			end if
			
! 			w=w1*w2
! 			intermedres=intermedres+w*kpntr(xx,yy,xi(i),eta(j),property%E,property%nu)
! 		end do
! 	end do
! 	placeholder = intermedres
! end subroutine gaussint








! subroutine funcselector(ipos,jpos) !Selects the correct function from the k-entries below
! 	integer, intent(in) :: ipos, jpos

! 	if (ipos==1) then 
! 		if (jpos==1) then
! 			kpntr=k11
! 		else if (jpos==2) then 
! 			kpntr=k12
! 		else if (jpos==3) then 
! 			kpntr=k13
! 		else if (jpos==4) then 
! 			kpntr=k14
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==2) then
! 		if (jpos==1) then
! 			kpntr=k21
! 		else if (jpos==2) then 
! 			kpntr=k22
! 		else if (jpos==3) then 
! 			kpntr=k23
! 		else if (jpos==4) then 
! 			kpntr=k24
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==3) then
! 		if (jpos==1) then
! 			kpntr=k31
! 		else if (jpos==2) then 
! 			kpntr=k32
! 		else if (jpos==3) then 
! 			kpntr=k33
! 		else if (jpos==4) then 
! 			kpntr=k34
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==4) then
! 		if (jpos==1) then
! 			kpntr=k41
! 		else if (jpos==2) then 
! 			kpntr=k42
! 		else if (jpos==3) then 
! 			kpntr=k43
! 		else if (jpos==4) then 
! 			kpntr=k44
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==5) then
! 		if (jpos==1) then
! 			kpntr=k51
! 		else if (jpos==2) then 
! 			kpntr=k52
! 		else if (jpos==3) then 
! 			kpntr=k53
! 		else if (jpos==4) then 
! 			kpntr=k54
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==6) then
! 		if (jpos==1) then
! 			kpntr=k61
! 		else if (jpos==2) then 
! 			kpntr=k62
! 		else if (jpos==3) then 
! 			kpntr=k63
! 		else if (jpos==4) then 
! 			kpntr=k64
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==7) then
! 		if (jpos==1) then
! 			kpntr=k71
! 		else if (jpos==2) then 
! 			kpntr=k72
! 		else if (jpos==3) then 
! 			kpntr=k73
! 		else if (jpos==4) then 
! 			kpntr=k74
! 		else
! 			call exit()
! 		endif
! 	else if (ipos==8) then
! 		if (jpos==1) then
! 			kpntr=k81
! 		else if (jpos==2) then 
! 			kpntr=k82
! 		else if (jpos==3) then 
! 			kpntr=k83
! 		else if (jpos==4) then 
! 			kpntr=k84
! 		else
! 			call exit()
! 		endif
! 	else
! 		call exit()
! 	endif
! end subroutine funcselector




! Those functions are copy-pasted from the Maple program and represent the entries of the element stiffness matrix before integration
function k(xx,yy,xi,eta,E,nu)
	implicit none
     real, dimension(8,8) :: k
	real, dimension(4), intent(in) :: xx, yy ! Those are the coordinates of the points involved
	real, intent(in)::xi, eta, E, nu

     k(:,:)=0

	 k(1,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1&
     & - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + n&
     &u) / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) +&
     & eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - et&
     &a * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx&
     &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     &+ xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(1,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     ! & / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     ! & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi *&
     ! & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     ! &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     ! &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     ! &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(1,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     ! &0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + &
     ! &xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) &
     ! &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     ! &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     ! &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     ! & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     ! &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     ! & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     ! &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     ! &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     ! & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     ! & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0&
     ! &.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx&
     ! &(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     ! &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     ! &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     ! &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     ! &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     ! &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     ! &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     ! &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(1,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     ! & / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     ! & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     ! &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     ! &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     ! & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     ! &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(1,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta &
     ! &* yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) &
     ! &+ yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2&
     ! &D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1&
     ! &) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     ! &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     ! &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     ! &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     ! &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     ! &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     ! &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     ! &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     ! &D1)))
     !  k(1,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.&
     ! &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * y&
     ! &y(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3&
     ! &) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (&
     ! &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     ! &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(1,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     ! & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     ! & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.&
     ! &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     ! &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     ! & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     ! & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     ! & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     ! & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     ! & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     ! & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     ! & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(1,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     ! &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     ! &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * &
     ! &(xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2)&
     ! & + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(2,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     ! & / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     ! & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi *&
     ! & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     ! &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     ! &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     ! &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(2,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (&
     ! &0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * y&
     ! &y(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) -&
     ! & xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx&
     ! &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     ! &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     ! &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     ! &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     ! &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     ! & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     ! &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     ! &+ xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D&
     ! &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(2,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     ! &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     ! &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     ! & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (et&
     ! &a * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2&
     ! &) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(2,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     ! &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     ! &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
     ! &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     ! & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2&
     ! &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     ! &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     ! &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     ! &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     ! &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     ! & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     ! &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     ! &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     ! &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)&
     ! &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     ! &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     ! &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1)&
     ! & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     ! &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     ! &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     ! &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     ! &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     ! &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     ! &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     ! &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     ! &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(2,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     ! &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     ! &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu)&
     ! & / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - &
     ! &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     ! &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta &
     ! &* xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) &
     ! &+ xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(2,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     ! & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     ! & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi&
     ! & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     ! &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D&
     ! &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     ! &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1))&
     ! & * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1&
     ! & * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy&
     ! &(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     ! &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     ! &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     ! &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     ! &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     ! &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     ! &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     ! &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) -&
     ! & 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1&
     ! &) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     ! &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     ! &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     ! &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     ! &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     ! &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     ! &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     ! &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     ! &D1)))
     !  k(2,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     ! &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     ! &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     ! & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta&
     ! & * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2)&
     ! & + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(2,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     ! & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     ! & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
     ! &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     ! & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2&
     ! &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     ! &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     ! &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     ! &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     ! &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     ! & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     ! &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     ! &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     ! &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)&
     ! &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     ! &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     ! &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) &
     ! &- 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(&
     ! &1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     ! & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     ! & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     ! & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     ! & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     ! & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     ! & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     ! & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(3,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     ! &0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + &
     ! &xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) &
     ! &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     ! &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     ! &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     ! & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     ! &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     ! & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     ! &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     ! &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     ! & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     ! & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0&
     ! &.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx&
     ! &(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     ! &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     ! &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     ! &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     ! &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     ! &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     ! &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     ! &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(3,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     ! &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     ! &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     ! & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (et&
     ! &a * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2&
     ! &) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(3,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu&
     ! &) / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + &
     ! &eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta&
     ! & * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(3,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (&
     ! &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     ! &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(3,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     ! &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * &
     ! &yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + &
     ! &yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1&
     ! & * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) &
     ! &- xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     ! &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     ! &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     ! &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     ! & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     ! &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     ! &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     ! &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     ! &)))
     !  k(3,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     ! &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     ! &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     ! & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     ! & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     ! &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     ! &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     ! &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     ! &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     ! &4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D&
     ! &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(&
     ! &1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) &
     ! &- yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi&
     ! & * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + &
     ! &xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(3,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     ! &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.&
     ! &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta *&
     ! & yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) +&
     ! & yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D&
     ! &1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1)&
     ! & - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     ! &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     ! &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     ! &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     ! &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     ! & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     ! & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     ! &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     ! &D1)))
     !  k(3,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     ! &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     ! & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(4,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     ! & / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     ! & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     ! &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     ! &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     ! & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     ! &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(4,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     ! &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     ! &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
     ! &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     ! & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2&
     ! &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     ! &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     ! &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     ! &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     ! &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     ! & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     ! &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     ! &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     ! &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)&
     ! &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     ! &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     ! &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1)&
     ! & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     ! &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     ! &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     ! &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     ! &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     ! &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     ! &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     ! &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     ! &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(4,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (&
     ! &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     ! &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(4,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0&
     ! &.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     ! &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(4,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     ! &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     ! &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     ! &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     ! &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     ! &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     ! &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     ! &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     ! &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     ! &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     ! &yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi&
     ! & * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) /&
     ! & (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - et&
     ! &a * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy&
     ! &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     ! &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     ! &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     ! &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     ! &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     ! & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     ! &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     ! &+ xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * &
     ! &xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + &
     ! &xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(4,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     ! &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     ! &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     ! &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     ! &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi &
     ! &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     ! &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1&
     ! & * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx&
     ! &(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     ! &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     ! &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     ! &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     ! &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     ! &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     ! &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     ! &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) *&
     ! & E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 *&
     ! & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     ! &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0&
     ! &.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) &
     ! &+ xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     ! &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     ! &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     ! &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     ! & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     ! &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     ! &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     ! &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     ! &)))
     !  k(4,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     ! &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     ! &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     ! &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) &
     ! &/ (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - e&
     ! &ta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + y&
     ! &y(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     ! & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     ! & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     ! & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     ! &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     ! &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     ! &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     ! & + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta *&
     ! & xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) +&
     ! & xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(4,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     ! &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     ! &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     ! &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     ! &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi&
     ! & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     ! &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D&
     ! &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     ! &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) &
     ! &* E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 &
     ! &* (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(&
     ! &2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - &
     ! &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     ! & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     ! &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     ! &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     ! &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     ! &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     ! & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     ! & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     ! &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     ! &D1)))
     !  k(5,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta &
     ! &* yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) &
     ! &+ yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2&
     ! &D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1&
     ! &) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     ! &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     ! &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     ! &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     ! &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     ! &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     ! &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     ! &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     ! &D1)))
     !  k(5,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     ! &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     ! &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu)&
     ! & / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - &
     ! &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     ! &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta &
     ! &* xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) &
     ! &+ xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(5,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     ! &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * &
     ! &yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + &
     ! &yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1&
     ! & * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) &
     ! &- xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     ! &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     ! &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     ! &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     ! & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     ! &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     ! &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     ! &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     ! &)))
     !  k(5,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     ! &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     ! &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     ! &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     ! &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     ! &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     ! &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     ! &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     ! &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     ! &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     ! &yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi&
     ! & * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) /&
     ! & (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - et&
     ! &a * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy&
     ! &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     ! &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     ! &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     ! &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     ! &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     ! & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     ! &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     ! &+ xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * &
     ! &xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + &
     ! &xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(5,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu)&
     ! & / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     ! &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     ! & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     ! & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     ! & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     ! &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     ! &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     ! &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     ! &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     ! &3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta *&
     ! & xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - &
     ! &nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(5,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     ! &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     ! &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     ! &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     ! &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     ! &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     ! &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     ! &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     ! &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     ! &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     ! &0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     ! &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * yy&
     ! &(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) /&
     ! & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     ! &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     ! &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     ! & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     ! &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     ! & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     ! &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     ! &) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     ! &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     ! &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1)&
     ! & - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - &
     ! &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi *&
     ! & xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx&
     ! &(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(5,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0&
     ! &.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D&
     ! &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * y&
     ! &y(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + y&
     ! &y(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 &
     ! &* (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) -&
     ! & xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     ! &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     ! &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     ! &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     ! &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     ! &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     ! &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     ! &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     ! &)))
     !  k(5,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     ! &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     ! &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     ! &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     ! &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     ! &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     ! &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     ! &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     ! &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     ! &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     ! &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     ! &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     ! &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     ! &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     ! & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     ! & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     ! &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     ! &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     ! &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     ! &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     ! &4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1)&
     ! & - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3)&
     ! & + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1&
     ! & - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1&
     ! &) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) -&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi &
     ! &* xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + x&
     ! &x(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(6,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.&
     ! &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * y&
     ! &y(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3&
     ! &) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (&
     ! &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     ! &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     ! &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     ! &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     ! &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     ! & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     ! &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     ! &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     ! &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(6,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     ! & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     ! & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi&
     ! & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     ! &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D&
     ! &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     ! &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1))&
     ! & * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1&
     ! & * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy&
     ! &(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     ! &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     ! &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     ! &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     ! &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     ! &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     ! &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     ! &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) -&
     ! & 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1&
     ! &) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     ! &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     ! &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     ! &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     ! &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     ! &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     ! &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     ! &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     ! &D1)))
     !  k(6,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     ! &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     ! &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     ! & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     ! & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     ! &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     ! &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     ! &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     ! &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     ! &4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D&
     ! &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(&
     ! &1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) &
     ! &- yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi&
     ! & * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + &
     ! &xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(6,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     ! &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     ! &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     ! &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     ! &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi &
     ! &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     ! &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1&
     ! & * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx&
     ! &(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     ! &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     ! &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     ! &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     ! &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     ! &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     ! &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     ! &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) *&
     ! & E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 *&
     ! & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     ! &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     ! &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     ! &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     ! &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     ! &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     ! &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     ! &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     ! &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0&
     ! &.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) &
     ! &+ xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     ! &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     ! &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     ! &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     ! & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     ! &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     ! &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     ! &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     ! &)))
     !  k(6,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     ! &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     ! &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     ! &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     ! &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     ! &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     ! &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     ! &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     ! &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     ! &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     ! &0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     ! &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * yy&
     ! &(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) /&
     ! & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     ! &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     ! &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     ! & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     ! &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     ! & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     ! &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     ! &) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     ! &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     ! &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     ! &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     ! &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     ! &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     ! & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     ! &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     ! &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     ! &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1)&
     ! & - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - &
     ! &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi *&
     ! & xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx&
     ! &(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
     !  k(6,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.&
     ! &1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy(&
     ! &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi&
     ! & * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - &
     ! &nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(6,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 -&
     ! & 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3)&
     ! & - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy&
     ! &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     ! &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     ! &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     ! &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     ! & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     ! &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     ! &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     ! &1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi *&
     ! & xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     ! &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     ! &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     ! &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     ! &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     ! &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     ! &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     ! &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     ! &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     ! &yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     ! &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / &
     ! &(0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta&
     ! & * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(6,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 +&
     ! & nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2)&
     ! & + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - e&
     ! &ta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + x&
     ! &x(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     ! & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     ! & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     ! & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     ! &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     ! &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     ! &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     ! & + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi &
     ! &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     ! &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     ! &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     ! &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * &
     ! &E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * &
     ! &(xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2)&
     ! & + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.&
     ! &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     ! & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     ! &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     ! &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     ! &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     ! &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     ! &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     ! &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     ! &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     ! &)))
     !  k(7,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     ! &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     ! &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     ! &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     ! & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     ! & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     ! & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     ! &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     ! &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     ! &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     ! &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     ! &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     ! &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     ! &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     ! & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     ! & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.&
     ! &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     ! &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     ! & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     ! & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     ! & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     ! & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     ! & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     ! & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     ! & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(7,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     ! & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     ! &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     ! &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     ! & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta&
     ! & * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2)&
     ! & + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(7,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     ! &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     ! &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     ! & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.&
     ! &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta *&
     ! & yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) +&
     ! & yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D&
     ! &1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1)&
     ! & - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     ! &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     ! &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     ! &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     ! &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     ! & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     ! & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     ! &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     ! &D1)))
     !  k(7,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     ! &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     ! &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     ! &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) &
     ! &/ (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - e&
     ! &ta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + y&
     ! &y(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     ! & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     ! & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     ! & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     ! &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     ! &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     ! &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     ! & + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta *&
     ! & xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) +&
     ! & xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(7,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0&
     ! &.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     ! &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     ! &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     ! &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     ! &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     ! &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     ! &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     ! &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     ! &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     ! &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     ! &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D&
     ! &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * y&
     ! &y(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + y&
     ! &y(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 &
     ! &* (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) -&
     ! & xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     ! &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     ! &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     ! &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     ! &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     ! &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     ! &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     ! &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     ! &)))
     !  k(7,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 -&
     ! & 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3)&
     ! & - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy&
     ! &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     ! &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     ! &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     ! &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     ! & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     ! &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     ! &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     ! &1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi *&
     ! & xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     ! &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     ! &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     ! &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     ! &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     ! &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     ! &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     ! &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     ! &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     ! &yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     ! &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / &
     ! &(0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta&
     ! & * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(&
     ! &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     ! &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     ! &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     ! & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     ! & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     ! &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     ! & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     ! & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * x&
     ! &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     ! &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(7,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu&
     ! &) / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + &
     ! &eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta &
     ! &* xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
     !  k(7,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     ! &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     ! & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(8,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     ! &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     ! &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     ! &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     ! &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     ! &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     ! &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     ! &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     ! &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     ! &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     ! &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     ! &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     ! &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     ! &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     ! &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     ! &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     ! & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     ! &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     ! & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     ! & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     ! &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     ! &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     ! &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     ! &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     ! & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     ! &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     ! & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     ! &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     ! &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     ! &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     ! &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     ! &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     ! &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     ! &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     ! &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     ! &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     ! &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     ! &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     ! &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     ! &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * &
     ! &(xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2)&
     ! & + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(8,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     ! & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     ! &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     ! &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     ! & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     ! & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     ! &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     ! &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     ! &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     ! &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     ! & * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     ! & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     ! & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
     ! &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     ! & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2&
     ! &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     ! &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     ! &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     ! &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     ! &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     ! & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     ! &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     ! &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     ! &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)&
     ! &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     ! &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     ! &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) &
     ! &- 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(&
     ! &1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     ! & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     ! & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     ! & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     ! & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     ! & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     ! & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     ! & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     ! &.4D1)))
     !  k(8,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     ! &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     ! &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     ! &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     ! &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     ! &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     ! & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     ! &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     ! &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     ! &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     ! &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     ! &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     ! & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     ! & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     ! & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     ! &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     ! & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     ! & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     ! & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     ! &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     ! & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(8,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     ! &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     ! & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     ! &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     ! &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     ! &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     ! &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     ! & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     ! &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     ! & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     ! &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     ! &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     ! &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     ! &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     ! &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     ! &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     ! &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     ! &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     ! & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     ! & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     ! &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     ! &) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi&
     ! & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     ! &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     ! &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     ! &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     ! &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     ! & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     ! &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     ! &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     ! &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D&
     ! &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     ! &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     ! & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     ! & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     ! & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     ! &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     ! &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     ! &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     ! & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) &
     ! &* E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 &
     ! &* (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(&
     ! &2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - &
     ! &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     ! & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     ! &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     ! &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     ! &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     ! &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     ! & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     ! & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     ! &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     ! &D1)))
     !  k(8,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     ! &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     ! &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     ! &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     ! &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     ! &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     ! &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     ! &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     ! &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     ! &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     ! &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     ! &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     ! &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     ! &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     ! &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     ! &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     ! &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     ! &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     ! &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     ! &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     ! & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     ! &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     ! &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     ! & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     ! & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     ! &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     ! &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     ! &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     ! &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     ! &4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1)&
     ! & - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3)&
     ! & + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1&
     ! & - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1&
     ! &) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) -&
     ! & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     ! & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     ! & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     ! &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     ! &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     ! &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     ! &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     ! &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi &
     ! &* xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + x&
     ! &x(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(8,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 +&
     ! & nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2)&
     ! & + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (x&
     ! &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     ! &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     ! &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     ! &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     ! & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     ! &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     ! & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     ! & yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - e&
     ! &ta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + x&
     ! &x(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     ! & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     ! & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     ! & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     ! &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     ! &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     ! &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     ! & + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi &
     ! &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     ! &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     ! & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     ! & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     ! & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     ! &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     ! &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     ! &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     ! & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     ! &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     ! &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     ! & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     ! & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     ! & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     ! &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     ! & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     ! & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     ! & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * &
     ! &E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * &
     ! &(xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2)&
     ! & + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     ! &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     ! &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     ! &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     ! &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     ! & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     ! & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     ! &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.&
     ! &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     ! & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     ! &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     ! &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     ! &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     ! &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     ! &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     ! &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     ! &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     ! &)))
     !  k(8,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     ! &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     ! &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     ! & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     ! & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     ! & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     ! &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     ! &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     ! & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     ! &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     ! &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     ! & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     ! &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     ! & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     ! & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     ! & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     ! & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     ! &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     ! & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     ! & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     ! & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     ! & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     ! &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     ! &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     ! &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     ! &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     ! & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     ! &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     ! & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     ! &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     ! &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     ! & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     ! &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     ! &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     ! & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     ! &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     ! &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     ! &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     ! &(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     ! &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     ! &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     ! &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     ! &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     ! &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     ! &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     ! &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     ! &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     ! &yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     ! &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     ! &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     ! & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     ! &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     ! &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     ! &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     ! &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     ! & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     ! & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     ! &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     ! &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     ! & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     ! &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     ! &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     ! &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     ! &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     ! &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     ! &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     ! &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
     !  k(8,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     ! &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     ! &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     ! &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     ! &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     ! &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     ! &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     ! &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     ! & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     ! & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     ! &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     ! &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     ! &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     ! &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     ! &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     ! & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     ! & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     ! &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     ! &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     ! &1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     ! &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     ! &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     ! &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     ! &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     ! & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     ! &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     ! &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     ! &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     ! &(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0&
     ! &.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy&
     ! &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     ! &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     ! &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     ! &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     ! &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     ! &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     ! & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     ! &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     ! &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     ! &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     ! &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     ! &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     ! &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     ! &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     ! &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     ! &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     ! &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     ! &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     ! &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))


end function k



end subroutine generateesm


