! part of the simpleFEM Program by Maximilian Hartig
! This is essentially a subroutine that generates the element stiffness matrix for a quadrilateral element
! Input: list  of Nodes involved

! We create the stiffness matrix for a four-node quadrilateral element
! The matrix is a 8x4 matrix. After matrix multiplication with the displacement vector, it should equal the force (or zero)
! k x u = f
! Where k=8x4, u=8x1={ux1,uy1,ux2,uy2,ux3,uy3,ux4,uy4} and f=8x1={fx1,fy1,...,fy4}
!

subroutine generateesm(k, element, property)

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
real, dimension(4,8), intent(inout)::k
real, dimension(8,4) :: kdummy ! I had to create this dummy because I made an error with the orientation of the element stiffness matrix and ended up calculating the transpose of k in Maple
type(propertytype), intent(in):: property

type(nodetype), dimension(4) :: node
real, dimension(8) :: localcoords
real, dimension(2,2) :: Jacoby
real, dimension(4) :: xx, yy
real::placeholder=0


! This pointer helps us select which function to integrate:

procedure(kxx), pointer :: kpntr => null()




!We create the vectors of nodal coordinates that are important for the evaluation of our integrals
do i=1,4
	xx(i)=element%node(i)%x
	yy(i)=element%node(i)%y
end do


do i=1,8
	do j=1,4
		placeholder=0
		call gaussint(placeholder,i,j)
		kdummy(i,j) = placeholder
	end do
end do

! We need to transpose the dummy element stiffness matrix to orient it the right way.
! This is necessary due to my mistake in creating the lengthy terms below
k=transpose(kdummy)


contains




subroutine gaussint(placeholder,ipos,jpos)

	implicit none
	real:: placeholder
	integer, intent(in):: ipos,jpos
	real, dimension(3) :: xi, eta
	real:: w1, w2, w
	integer::i,j
	real:: intermedres=0

	call funcselector(ipos,jpos)

	xi=(/sqrt(3./5.),-sqrt(3./5.),0./)
	eta=(/sqrt(3./5.),-sqrt(3./5.),0./)
	intermedres=0

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
			intermedres=intermedres+w*kpntr(xx,yy,xi(i),eta(j),property%E,property%nu)
		end do
	end do
	placeholder = intermedres
end subroutine gaussint

subroutine funcselector(ipos,jpos) !Selects the correct function from the k-entries below
	integer, intent(in) :: ipos, jpos

	if (ipos==1) then 
		if (jpos==1) then
			kpntr=k11
		else if (jpos==2) then 
			kpntr=k12
		else if (jpos==3) then 
			kpntr=k13
		else if (jpos==4) then 
			kpntr=k14
		else
			call exit()
		endif
	else if (ipos==2) then
		if (jpos==1) then
			kpntr=k21
		else if (jpos==2) then 
			kpntr=k22
		else if (jpos==3) then 
			kpntr=k23
		else if (jpos==4) then 
			kpntr=k24
		else
			call exit()
		endif
	else if (ipos==3) then
		if (jpos==1) then
			kpntr=k31
		else if (jpos==2) then 
			kpntr=k32
		else if (jpos==3) then 
			kpntr=k33
		else if (jpos==4) then 
			kpntr=k34
		else
			call exit()
		endif
	else if (ipos==4) then
		if (jpos==1) then
			kpntr=k41
		else if (jpos==2) then 
			kpntr=k42
		else if (jpos==3) then 
			kpntr=k43
		else if (jpos==4) then 
			kpntr=k44
		else
			call exit()
		endif
	else if (ipos==5) then
		if (jpos==1) then
			kpntr=k51
		else if (jpos==2) then 
			kpntr=k52
		else if (jpos==3) then 
			kpntr=k53
		else if (jpos==4) then 
			kpntr=k54
		else
			call exit()
		endif
	else if (ipos==6) then
		if (jpos==1) then
			kpntr=k61
		else if (jpos==2) then 
			kpntr=k62
		else if (jpos==3) then 
			kpntr=k63
		else if (jpos==4) then 
			kpntr=k64
		else
			call exit()
		endif
	else if (ipos==7) then
		if (jpos==1) then
			kpntr=k71
		else if (jpos==2) then 
			kpntr=k72
		else if (jpos==3) then 
			kpntr=k73
		else if (jpos==4) then 
			kpntr=k74
		else
			call exit()
		endif
	else if (ipos==8) then
		if (jpos==1) then
			kpntr=k81
		else if (jpos==2) then 
			kpntr=k82
		else if (jpos==3) then 
			kpntr=k83
		else if (jpos==4) then 
			kpntr=k84
		else
			call exit()
		endif
	else
		call exit()
	endif
end subroutine funcselector




! Those functions are copy-pasted from the Maple program and represent the entries of the element stiffness matrix before integration
real function k11(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu

	k11 = ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4)&
     & - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 &
     &* (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) -&
     & yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D&
     &1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4&
     &) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & xi / 0.4D1)) ** 2) * (nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi&
     & * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(&
     &2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) &
     &- eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) &
     &+ yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (e&
     &ta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(&
     &2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) +&
     & (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (&
     &xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) &
     &+ yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1&
     &))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 -&
     & xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) &
     &* yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1)&
     & * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) /&
     & 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 &
     &+ xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4&
     &) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(&
     &3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + &
     &xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.&
     &8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3&
     &) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k11

real function k21(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu

	k21 = ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4&
     &D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu &
     &/ 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 /&
     & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(&
     &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 *&
     & (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - &
     &yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1&
     &) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4)&
     & - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &xi / 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi * &
     &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) **&
     & 2) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 -&
     & xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) &
     &* yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1)&
     & * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) /&
     & 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 &
     &+ xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4&
     &) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(&
     &3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + &
     &xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.&
     &8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3&
     &) / 0.8D1) * E / (-nu ** 2 + 0.1D1)
end function k21


real function k31(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k31 = ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4)&
     & - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + et&
     &a * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta *&
     & xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (x&
     &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2&
     &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))&
     & * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4&
     &) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) -&
     & eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 *&
     & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - &
     &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4&
     &D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3&
     &) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy&
     &(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1&
     &) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D&
     &1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(&
     &2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * e&
     &ta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) *&
     & xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 +&
     & xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0&
     &.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(&
     &1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k31

real function k41(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k41= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta&
     & * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2&
     &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi&
     & * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (&
     &eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy&
     &(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) &
     &+ 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) -&
     & xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi&
     & / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(&
     &1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / &
     &0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4&
     &) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - x&
     &i * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 /&
     & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3&
     &) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy&
     &(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1&
     &) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D&
     &1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(&
     &2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * e&
     &ta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) *&
     & xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 +&
     & xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0&
     &.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(&
     &1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k41

real function k51(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k51= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4)&
     & - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + et&
     &a / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - et&
     &a * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 /&
     & 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta&
     & * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * x&
     &x(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (0.2D1 * (xi &
     &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + et&
     &a / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - et&
     &a * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 /&
     & 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (x&
     &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2&
     &D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + &
     &xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)&
     &) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(&
     &4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 +&
     & eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) -&
     & eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1&
     & / 0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * &
     &eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) &
     &* xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / &
     &0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + &
     &xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) *&
     & eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta *&
     & yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi &
     &/ 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(&
     &1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1&
     & + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) /&
     & 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k51

real function k61(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k61= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D&
     &1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3&
     &) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + x&
     &i * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx&
     &(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / &
     &0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta &
     &* yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * x&
     &x(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta&
     & * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta&
     & * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta&
     & * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3&
     &) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) *&
     & xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi&
     & - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) +&
     & xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * &
     &(-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2)&
     & + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi *&
     & xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 * (et&
     &a * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2&
     &) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + &
     &0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - x&
     &x(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta&
     & * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta&
     & * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta&
     & * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1)&
     & * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * &
     &xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(&
     &1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3)&
     & * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi /&
     & 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0&
     &.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) &
     &+ yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4&
     &D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) -&
     & xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi&
     & / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi *&
     & yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) &
     &- xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1&
     & / 0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * &
     &eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) &
     &* xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / &
     &0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + &
     &xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) *&
     & eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta *&
     & yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi &
     &/ 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(&
     &1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1&
     & + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) /&
     & 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k61

real function k71(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k71= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4)&
     & - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + et&
     &a * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * &
     &xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (x&
     &i * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) +&
     & yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D&
     &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))&
     & * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4&
     &) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) -&
     & eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 *&
     & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - &
     &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4&
     &D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D&
     &1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3&
     &) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3&
     &) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy&
     &(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1&
     &) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D&
     &1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(&
     &2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * e&
     &ta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) *&
     & xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 +&
     & xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0&
     &.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(&
     &1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k71

real function k81(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k81= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * &
     &yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx&
     &(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx&
     &(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx&
     &(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2&
     &) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) *&
     & yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy&
     &(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2&
     &) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D&
     &1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3&
     &) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta&
     & * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2&
     &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (&
     &eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy&
     &(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) &
     &+ 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) -&
     & xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi&
     & / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(&
     &1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / &
     &0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4&
     &) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - x&
     &i * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(&
     &3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3&
     &) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy&
     &(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1&
     &) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D&
     &1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(&
     &2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * e&
     &ta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) *&
     & xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 +&
     & xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0&
     &.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(&
     &1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k81

real function k12(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k12= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + &
     &xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) &
     &* (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4)&
     & - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - &
     &eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta &
     &* xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (&
     &xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) &
     &+ yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1&
     &)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy&
     &(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1&
     & - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3)&
     & - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * y&
     &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     &D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1&
     & * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) &
     &- yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4&
     &D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(&
     &4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 &
     &+ xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(&
     &3) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * y&
     &y(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(&
     &1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8&
     &D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx&
     &(2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * &
     &eta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) &
     &* xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 &
     &+ xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / &
     &0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy&
     &(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k12

real function k22(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k22= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) &
     &+ xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi *&
     & xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu &
     &/ 0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * y&
     &y(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(&
     &1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(&
     &3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(&
     &4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2)&
     & * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * &
     &yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(&
     &2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2)&
     & - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D&
     &1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) -&
     & xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta &
     &* xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * &
     &(eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - y&
     &y(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1)&
     & + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) &
     &- xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + x&
     &i / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) -&
     & eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + et&
     &a * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) &
     &* (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy&
     &(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta /&
     & 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(&
     &4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 &
     &+ xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(&
     &3) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * y&
     &y(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(&
     &1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8&
     &D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx&
     &(2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * &
     &eta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) &
     &* xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 &
     &+ xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / &
     &0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy&
     &(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k22

real function k32(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k32= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 &
     &* (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) -&
     & yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D&
     &1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4&
     &) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - &
     &xi / 0.4D1)) ** 2) * (nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi &
     &* yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2&
     &) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - &
     &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta&
     & * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2)&
     & + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0&
     &.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) &
     &+ eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - e&
     &ta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + x&
     &x(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     & + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) *&
     & (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1&
     &) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(&
     &4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi&
     & / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D&
     &1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(&
     &1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * e&
     &ta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * &
     &xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2)&
     & * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 -&
     & xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0&
     &.8D1) * E / (-nu ** 2 + 0.1D1)

end function k32

real function k42(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k42= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta&
     & * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2)&
     & + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (e&
     &ta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(&
     &2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) +&
     & 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - &
     &xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi /&
     & 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi * yy(1&
     &) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) -&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi&
     & * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + &
     &xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2) *&
     & (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1&
     &) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(&
     &4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi&
     & / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D&
     &1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(&
     &1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * e&
     &ta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * &
     &xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2)&
     & * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 -&
     & xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0&
     &.8D1) * E / (-nu ** 2 + 0.1D1)

end function k42

real function k52(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k52= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - et&
     &a / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - et&
     &a * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta &
     &* yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx&
     &(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * (&
     &-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - &
     &yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi &
     &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1&
     & * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx&
     &(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + et&
     &a / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - et&
     &a * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k52

real function k62(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k62= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / &
     &0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta &
     &* yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * x&
     &x(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta&
     & * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta&
     & * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta&
     & * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3&
     &) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) *&
     & xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi&
     & - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) +&
     & xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * &
     &(0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) &
     &+ xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0&
     &.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx&
     &(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     &.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta &
     &* yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2&
     &D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + &
     &yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1&
     &) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - x&
     &x(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta&
     & * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta&
     & * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta&
     & * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1)&
     & * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * &
     &xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(&
     &1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3)&
     & * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / &
     &0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy&
     &(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 &
     &+ eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - x&
     &i * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k62

real function k72(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k72= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta&
     & * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * x&
     &x(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (xi&
     & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1&
     & * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx&
     &(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) &
     &* (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4)&
     & - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - &
     &eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - &
     &eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     &/ 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (&
     &xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) &
     &+ yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)&
     &) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(&
     &4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 -&
     & eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) -&
     & eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k72

real function k82(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k82= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta&
     & * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2)&
     & + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi *&
     & xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (e&
     &ta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(&
     &2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) +&
     & 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - &
     &xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi /&
     & 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (&
     &0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1)&
     & + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.&
     &4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) &
     &- xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi&
     & / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi *&
     & yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) &
     &- xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k82

real function k13(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k13= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + et&
     &a * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * &
     &xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (0.2D1 * (xi&
     & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D&
     &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1))&
     & * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4&
     &) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + &
     &eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - &
     &eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * &
     &(eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - y&
     &y(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1)&
     & + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) &
     &- xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + x&
     &i / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi &
     &* yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3)&
     & - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy&
     &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     &1 / 0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k13

real function k23(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k23= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) &
     &+ xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy&
     &(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1&
     & + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - &
     &xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + et&
     &a * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * &
     &xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 * (e&
     &ta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(&
     &2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) +&
     & 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - &
     &xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi &
     &/ 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (&
     &0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1)&
     & + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.&
     &4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) &
     &- xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + x&
     &i / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi &
     &* yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3)&
     & - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy&
     &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     &1 / 0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k23

real function k33(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k33 = ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - et&
     &a / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - et&
     &a * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta &
     &* yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx&
     &(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (et&
     &a * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2&
     &) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + &
     &0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - x&
     &x(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta&
     & * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta&
     & * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta&
     & * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1)&
     & * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * &
     &xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(&
     &1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3)&
     & * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / &
     &0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy&
     &(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 &
     &+ eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - x&
     &i * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k33

real function k43(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k43 = ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / &
     &0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(&
     &4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 &
     &- eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - x&
     &i * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta &
     &* yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx&
     &(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0&
     &.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx&
     &(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     &.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta &
     &* yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2&
     &D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + &
     &yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1&
     &) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - x&
     &x(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta&
     & * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta&
     & * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta&
     & * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1)&
     & * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * &
     &xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(&
     &1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3)&
     & * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / &
     &0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy&
     &(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 &
     &+ eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - x&
     &i * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 + xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k43

real function k53(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k53= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + xi / 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 *&
     & (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - &
     &yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1)&
     & + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) &
     &- xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi&
     & / 0.4D1)) ** 2) * (nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * &
     &yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta&
     & * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta&
     & * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta&
     & * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3&
     &) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) *&
     & xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi&
     & - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) +&
     & xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * &
     &(0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) +&
     & xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta&
     & * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1&
     & / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + et&
     &a * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * &
     &xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1&
     &) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) -&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi &
     &* xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + x&
     &x(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (-xx(&
     &1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1) * et&
     &a * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(4) * x&
     &i / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8&
     &D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx&
     &(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(1) * e&
     &ta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * eta * y&
     &y(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * xi / 0&
     &.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2) * yy(&
     &1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3)&
     & * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) &
     &* E / (-nu ** 2 + 0.1D1)

end function k53

real function k63(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k63= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx&
     &(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0&
     &.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx&
     &(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1))) * (nu * (-0.2D1 * (eta *&
     & yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) +&
     & yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D&
     &1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1)&
     & - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D&
     &1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi * yy(1) - x&
     &i * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(&
     &1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) &
     &- xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2) * (-xx(&
     &1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1) * et&
     &a * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(4) * x&
     &i / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8&
     &D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx&
     &(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(1) * e&
     &ta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * eta * y&
     &y(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * xi / 0&
     &.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2) * yy(&
     &1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3)&
     & * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) &
     &* E / (-nu ** 2 + 0.1D1)

end function k63

real function k73(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k73= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta *&
     & yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(&
     &2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi &
     &* yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + y&
     &y(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k73

real function k83(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k83= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / &
     &0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta &
     &* yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.&
     &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.&
     &4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2&
     &D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + &
     &yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1)&
     & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0&
     &.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(&
     &4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k83

real function k14(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k14 = ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 +&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
     &(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + x&
     &i * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) &
     &* (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4)&
     & - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + &
     &eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - &
     &eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta *&
     & xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (&
     &xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) &
     &+ yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1&
     &)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy&
     &(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 &
     &- eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) &
     &- eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy&
     &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1&
     &D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1&
     & * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) &
     &- yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4&
     &D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(&
     &4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 &
     &+ xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(&
     &3) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * y&
     &y(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(&
     &1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8&
     &D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx&
     &(2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * &
     &eta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) &
     &* xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 &
     &+ xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / &
     &0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy&
     &(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k14

real function k24(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k24= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) &
     &+ xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * &
     &xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu &
     &/ 0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * y&
     &y(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(&
     &1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(&
     &3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(&
     &4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2)&
     & * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * &
     &yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(&
     &2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2)&
     & - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D&
     &1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) -&
     & xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta *&
     & xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * &
     &(eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - y&
     &y(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1)&
     & + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) &
     &- xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + x&
     &i / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) -&
     & eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1&
     & / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) &
     &* (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy&
     &(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta /&
     & 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(&
     &4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 &
     &+ xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - &
     &xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(&
     &3) * eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * y&
     &y(2) * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(&
     &1) / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8&
     &D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx&
     &(2) * eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * &
     &eta * yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) &
     &* xi / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 &
     &+ xx(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / &
     &0.8D1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy&
     &(1) / 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k24

real function k34(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k34= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
     &3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi&
     & * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) *&
     & (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) &
     &- yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta&
     & * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * x&
     &x(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (xi&
     & * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + &
     &yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D&
     &1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + x&
     &x(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) &
     &* (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4)&
     & - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - e&
     &ta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - e&
     &ta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * &
     &(eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - y&
     &y(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1)&
     & + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) &
     &- xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi&
     & / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi *&
     & yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) &
     &- xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k34

real function k44(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k44 = ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.&
     &4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx&
     &(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy&
     &(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1&
     &) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3&
     &) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4&
     &) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) &
     &* yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * y&
     &y(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2&
     &) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) &
     &- xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1&
     & - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - &
     &xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta&
     & * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * x&
     &x(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (e&
     &ta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(&
     &2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) +&
     & 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - &
     &xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi /&
     & 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
     &a * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) &
     &- xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) &
     &- xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) &
     &- xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - &
     &xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(&
     &3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) &
     &* yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * &
     &yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / &
     &0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * &
     &xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (&
     &0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1)&
     & + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.&
     &4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) &
     &- xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * &
     &eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * &
     &eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * &
     &eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy&
     &(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4)&
     & * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + &
     &xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx&
     &(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi&
     & / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi *&
     & yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) &
     &- xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(&
     &3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(&
     &4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(&
     &1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi&
     & - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - &
     &xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(&
     &1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3)&
     & * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D&
     &1 / 0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) *&
     & eta * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2)&
     & * xi / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) /&
     & 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 +&
     & xx(3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) &
     &* eta * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta &
     &* yy(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi&
     & / 0.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx&
     &(1) * yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D&
     &1 + xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) &
     &/ 0.8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k44

real function k54(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k54= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.&
     &4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     &D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * &
     &(-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) -&
     & yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * e&
     &ta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * e&
     &ta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * e&
     &ta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(&
     &1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) &
     &* xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + x&
     &x(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(&
     &3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta&
     & / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta&
     & * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta *&
     & yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(&
     &2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * (&
     &-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - &
     &yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * et&
     &a * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * et&
     &a * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * et&
     &a * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1&
     &) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) *&
     & xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx&
     &(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3&
     &) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta &
     &/ 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta &
     &* xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (et&
     &a * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2&
     &) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0&
     &.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx&
     &(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0&
     &.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(&
     &4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k54

real function k64(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k64= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + &
     &xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) &
     &* eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) &
     &* eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) &
     &* eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) *&
     & yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy&
     &(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3)&
     & * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy&
     &(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3&
     &)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx&
     &(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / &
     &0.2D1) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(&
     &4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 +&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 + xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta *&
     & yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(&
     &2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.&
     &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.&
     &4D1)) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2&
     &D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + &
     &yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4&
     &) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1&
     &) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2&
     &) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi +&
     & xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx&
     &(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy&
     &(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4&
     &) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1)&
     & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0&
     &.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(&
     &4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1)&
     & * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3)&
     & * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4)&
     & * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) *&
     & yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy&
     &(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2)&
     & + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) -&
     & xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 -&
     & eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi&
     & * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / &
     &0.4D1 - xi / 0.4D1))) * (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta&
     & * yy(1) / 0.8D1 - xx(1) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * x&
     &i / 0.8D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8&
     &D1 - xx(2) * yy(1) * xi / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(&
     &3) * eta * yy(2) / 0.8D1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * et&
     &a * yy(4) / 0.8D1 + xx(1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy&
     &(3) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0&
     &.8D1 + xx(4) * yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) &
     &* yy(4) / 0.8D1 + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + &
     &xx(3) * yy(2) / 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.&
     &8D1 + xx(4) * yy(3) / 0.8D1) * E / (-nu ** 2 + 0.1D1)

end function k64

real function k74(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k74= ((0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4&
     &) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) &
     &* eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) &
     &* eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) &
     &* eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * &
     &yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(&
     &4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) &
     &+ xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - &
     &xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - &
     &eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi &
     &* xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - &
     &xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - &
     &xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - &
     &xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx&
     &(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3)&
     & * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * &
     &yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy&
     &(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0&
     &.4D1 - xi / 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 &
     &* (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) -&
     & yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1&
     &) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4)&
     & - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) *&
     & eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) *&
     & eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) *&
     & eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * y&
     &y(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4&
     &) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) +&
     & xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - x&
     &x(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - &
     &xi / 0.4D1)) ** 2) * (nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi &
     &* yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2)&
     & + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * (-0.2D1 * (eta * yy(1) - &
     &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta &
     &* xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) &
     &+ xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0&
     &.1D1 / 0.2D1 - nu / 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) &
     &+ eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - et&
     &a * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx&
     &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     &+ xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) *&
     & (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1&
     &) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(&
     &4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi&
     & / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D&
     &1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(&
     &1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * e&
     &ta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * &
     &xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2)&
     & * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 -&
     & xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0&
     &.8D1) * E / (-nu ** 2 + 0.1D1)

end function k74

real function k84(xx,yy,xi,eta,E,nu)
	implicit none
	real, dimension(4) :: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu
	k84= ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta *&
     & yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - x&
     &x(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - x&
     &x(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - x&
     &x(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(&
     &2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) &
     &* yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * y&
     &y(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(&
     &2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4&
     &D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(&
     &3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2) +&
     & xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * x&
     &x(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.1D1 / 0.2D1 - nu /&
     & 0.2D1) * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - eta&
     & * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) -&
     & xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) -&
     & xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) -&
     & xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - x&
     &x(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3&
     &) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) *&
     & yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * y&
     &y(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0&
     &.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta * x&
     &x(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta&
     & * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta&
     & * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta&
     & * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3&
     &) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) *&
     & xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi&
     & - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) +&
     & xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * &
     &(-0.1D1 / 0.4D1 - xi / 0.4D1)) * (0.2D1 * (xi * yy(1) - xi * yy(2)&
     & + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi *&
     & xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1))) * (nu * (-0.2D1 * (e&
     &ta * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(&
     &2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + &
     &0.2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - x&
     &x(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta&
     & * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta&
     & * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta&
     & * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1)&
     & * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * &
     &xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(&
     &1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3)&
     & * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi /&
     & 0.4D1)) ** 2 + (0.1D1 / 0.2D1 - nu / 0.2D1) * (0.2D1 * (xi * yy(1&
     &) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) -&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi &
     &* xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + x&
     &x(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2) *&
     & (-xx(1) * yy(3) * xi / 0.8D1 - xx(3) * eta * yy(1) / 0.8D1 - xx(1&
     &) * eta * yy(4) / 0.8D1 + xx(1) * yy(2) * xi / 0.8D1 - xx(3) * yy(&
     &4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx(2) * yy(1) * xi&
     & / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D&
     &1 + xx(2) * yy(4) * xi / 0.8D1 + xx(2) * eta * yy(4) / 0.8D1 + xx(&
     &1) * eta * yy(3) / 0.8D1 - xx(2) * eta * yy(3) / 0.8D1 - xx(4) * e&
     &ta * yy(2) / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * yy(3) * &
     &xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1 + xx(2)&
     & * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) / 0.8D1 -&
     & xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * yy(3) / 0&
     &.8D1) * E / (-nu ** 2 + 0.1D1)

end function k84



end subroutine generateesm


