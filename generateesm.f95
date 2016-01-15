! part of the simpleFEM Program by Maximilian Hartig
! This is essentially a subroutine that generates the element stiffness matrix for a quadrilateral element
! Input: list  of Nodes involved



subroutine generateesm(kk, element, property)

use mytypes


type(elementtype), intent(in) :: element
real, dimension(8,8), intent(inout)::kk
type(propertytype), intent(in):: property

type(nodetype), dimension(4) :: node
real, dimension(8) :: localcoords
real, dimension(2,2) :: Jacoby
real, dimension(4) :: xx, yy


do i=1,4
     xx(i)=element%node(i)%x
     yy(i)=element%node(i)%y
end do



call gaussint(kk)



contains




subroutine gaussint(dummymatrix)

	implicit none
	real, dimension(8,8), intent(inout) :: dummymatrix
	real, dimension(3) :: xi, eta
	real:: w1, w2, w
	integer::i,j
	real, dimension(8,8):: intermedres=0

!	call funcselector(ipos,jpos)

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
			intermedres=intermedres+w*k(xx,yy,xi(i),eta(j),property%E,property%nu)
		end do
	end do
	dummymatrix = intermedres
end subroutine gaussint



! Those functions are copy-pasted from the Maple program and represent the entries of the element stiffness matrix before integration
function k(xx,yy,xi,eta,E,nu)

     !!!!!!!!! THIS IS FOR PLAIN STRAIN !!!!!!!!

	implicit none
     real, dimension(8,8) :: k
	real, dimension(4):: xx, yy ! Those are the coordinates of the points involved
	real::xi, eta, E, nu



	 k(1,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
      k(1,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(1,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
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
     &.4D1)))
      k(1,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
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
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(1,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta &
     &* yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) &
     &+ yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2&
     &D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1&
     &) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     &D1)))
      k(1,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.&
     &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * y&
     &y(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3&
     &) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (&
     &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(1,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.&
     &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(1,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * &
     &(xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2)&
     & + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(2,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1&
     & / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta&
     & * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) *&
     & eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) *&
     & eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) *&
     & eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * &
     &yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(&
     &1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) &
     &* xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(&
     &3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)&
     &) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi *&
     & yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy&
     &(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 &
     &* (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(&
     &2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(2,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (&
     &0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * y&
     &y(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) -&
     & xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx&
     &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     &+ xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(2,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (et&
     &a * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2&
     &) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(2,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
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
     &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1)&
     & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(2,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
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
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu)&
     & / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - &
     &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta &
     &* xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) &
     &+ xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(2,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi&
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
     & * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1&
     & * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy&
     &(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) -&
     & 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1&
     &) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     &D1)))
      k(2,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta&
     & * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2)&
     & + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(2,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
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
     &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) &
     &- 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(&
     &1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(3,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
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
     &.4D1)))
      k(3,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (et&
     &a * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2&
     &) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(3,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu&
     &) / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + &
     &eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta&
     & * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(3,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (&
     &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(3,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * &
     &yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + &
     &yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1&
     & * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) &
     &- xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     &)))
      k(3,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
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
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(&
     &1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) &
     &- yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi&
     & * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + &
     &xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(3,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.&
     &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta *&
     & yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) +&
     & yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D&
     &1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1)&
     & - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     &D1)))
      k(3,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(4,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
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
     &) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 *&
     & (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2&
     &) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(4,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
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
     &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1)&
     & - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx&
     &(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta &
     &* yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta &
     &* yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta &
     &* yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) &
     &* xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * x&
     &i - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1&
     &) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) &
     &* yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(4,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (&
     &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(4,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0&
     &.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(4,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     &D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
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
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) /&
     & (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - et&
     &a * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy&
     &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     &+ xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * &
     &xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + &
     &xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(4,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi &
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
     & E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 *&
     & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0&
     &.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) &
     &+ xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     &)))
      k(4,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
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
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) &
     &/ (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - e&
     &ta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + y&
     &y(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     & + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta *&
     & xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) +&
     & xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(4,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi&
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
     &* E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 &
     &* (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(&
     &2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - &
     &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     &D1)))
      k(5,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta &
     &* yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) &
     &+ yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2&
     &D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1&
     &) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     &D1)))
      k(5,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
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
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu)&
     & / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - &
     &eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta &
     &* xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) &
     &+ xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(5,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * &
     &yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + &
     &yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1&
     & * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) &
     &- xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     &)))
      k(5,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
     &) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * y&
     &y(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * y&
     &y(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * y&
     &y(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * &
     &xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi &
     &- xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - x&
     &x(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(&
     &3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1&
     &D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi &
     &* xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
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
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) /&
     & (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - et&
     &a * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy&
     &(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * &
     &yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * &
     &yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) &
     &* xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * x&
     &i + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi +&
     & xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1&
     &) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) &
     &+ xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * &
     &xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + &
     &xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(5,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu)&
     & / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + e&
     &ta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1)&
     & * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2)&
     & * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4)&
     & * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) &
     &* yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * y&
     &y(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3&
     &) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * y&
     &y(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(&
     &3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta *&
     & xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - &
     &nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(5,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
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
     & (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * yy&
     &(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1)&
     & - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi *&
     & xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx&
     &(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(5,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0&
     &.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
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
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * y&
     &y(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + y&
     &y(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 &
     &* (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) -&
     & xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     &)))
      k(5,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
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
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1)&
     & - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3)&
     & + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1&
     & - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1&
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
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(6,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
     &eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3&
     &) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4&
     &) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1&
     &) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi &
     &- xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - x&
     &x(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1&
     &) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) &
     &* yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 &
     &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta &
     &* xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * &
     &eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * &
     &eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * &
     &eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * y&
     &y(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1&
     &) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) *&
     & xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3&
     &) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3))&
     & * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.&
     &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * y&
     &y(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3&
     &) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (&
     &xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) &
     &+ xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx&
     &(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx&
     &(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx&
     &(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2)&
     & * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * &
     &yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + &
     &xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx&
     &(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(6,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi&
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
     & * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1&
     & * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy&
     &(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) &
     &- xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) &
     &+ xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) &
     &+ xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + x&
     &x(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4&
     &) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4&
     &) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) &
     &- xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) -&
     & 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1&
     &) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * &
     &yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * &
     &yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * &
     &yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * &
     &xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi &
     &- xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) &
     &* yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * &
     &yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4&
     &D1)))
      k(6,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
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
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(&
     &1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) &
     &- yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi&
     & * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + &
     &xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2&
     &) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3&
     &) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1&
     &) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) *&
     & yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy&
     &(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx&
     &(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4&
     &) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(6,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) + (0.2D1 * (xi &
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
     & E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 *&
     & (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2&
     &) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - &
     &xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + &
     &xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + &
     &xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(&
     &2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) &
     &* yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) &
     &+ xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - &
     &xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0&
     &.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) &
     &+ xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy&
     &(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy&
     &(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy&
     &(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi&
     & + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - &
     &xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * &
     &yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy&
     &(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1&
     &)))
      k(6,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
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
     & (0.1D1 / 0.4D1 + xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * yy&
     &(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) /&
     & (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) &
     &+ xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) &
     &+ xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi -&
     & xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx&
     &(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4)&
     & * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx&
     &(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4&
     &) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) &
     &- eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) &
     &+ xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * et&
     &a * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * et&
     &a * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy&
     &(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4)&
     & * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * &
     &xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * &
     &yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy&
     &(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1)&
     & - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - &
     &yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi *&
     & xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx&
     &(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)))
      k(6,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.&
     &1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy(&
     &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi&
     & * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) ** 2 * E * (0.1D1 - &
     &nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(6,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 -&
     & 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3)&
     & - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy&
     &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     &1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi *&
     & xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / &
     &(0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta&
     & * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(6,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 +&
     & nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2)&
     & + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - e&
     &ta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + x&
     &x(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     & + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi &
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
     &E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * &
     &(xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2)&
     & + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     &)))
      k(7,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / &
     &(0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy&
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
     &) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta &
     &* yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (&
     &0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta&
     & * yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2)&
     & + yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.&
     &2D1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(&
     &1) - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(7,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1&
     & - 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - &
     &xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu&
     &) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) -&
     & eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) +&
     & yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta&
     & * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2)&
     & + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(7,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (&
     &0.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(&
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.&
     &1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta *&
     & yy(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) +&
     & yy(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D&
     &1 * (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1)&
     & - xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     &D1)))
      k(7,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 &
     &- 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
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
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(&
     &2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) &
     &/ (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - e&
     &ta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + y&
     &y(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     & + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta *&
     & xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) +&
     & xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(7,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0&
     &.1D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3&
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
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D&
     &1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * y&
     &y(1) - eta * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + y&
     &y(3) + yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 &
     &* (eta * xx(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) -&
     & xx(2) + xx(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     &)))
      k(7,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 -&
     & 0.2D1 * nu) * nu * (0.2D1 * (xi * yy(1) - xi * yy(2) + xi * yy(3)&
     & - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx(1) * eta * yy&
     &(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy&
     &(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy&
     &(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * x&
     &i - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi -&
     & xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx&
     &(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3&
     &) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D&
     &1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi *&
     & xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * et&
     &a * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * et&
     &a * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * et&
     &a * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(&
     &3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) &
     &* xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * x&
     &i - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) &
     &+ xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) *&
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi * yy(1) - xi * yy(2&
     &) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (xx&
     &(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx&
     &(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx&
     &(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(&
     &1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) &
     &* yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * y&
     &y(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) &
     &* yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * &
     &yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi &
     &* xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / &
     &(0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta&
     & * yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(&
     &4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * y&
     &y(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * y&
     &y(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) *&
     & xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi&
     & + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + &
     &xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1)&
     & - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) +&
     & xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * x&
     &x(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + x&
     &x(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2)&
     & * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3)&
     & * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1)&
     & * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * &
     &yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(&
     &2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(&
     &2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4)&
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(7,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0.1D1 + nu&
     &) / (0.1D1 - 0.2D1 * nu) + (-0.2D1 * (eta * yy(1) - eta * yy(2) + &
     &eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta &
     &* xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))
      k(7,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(8,1) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &-0.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 &
     &* nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - &
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
     & * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta *&
     & yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)&
     &) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(&
     &3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(&
     &2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * x&
     &i - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi +&
     & xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx&
     &(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) -&
     & xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + x&
     &x(4) * yy(3)) * (-0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx&
     &(1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx&
     &(3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) &
     &* eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) &
     &* eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) &
     &* yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * y&
     &y(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2&
     &) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2&
     &) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) &
     &* yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0&
     &.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * &
     &yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(&
     &3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * &
     &(xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2)&
     & + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(8,2) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (-0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1&
     & + nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(&
     &2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / &
     &(xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) +&
     & xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) +&
     & xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - &
     &xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(&
     &3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) &
     &* yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(&
     &2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4)&
     & * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) -&
     & eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) +&
     & xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta&
     & * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta&
     & * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(&
     &2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) &
     &* xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * x&
     &i + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * y&
     &y(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(&
     &1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (x&
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
     &) * E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D&
     &1 * (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + y&
     &y(2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4)&
     & - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1)&
     & + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2)&
     & + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + &
     &xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(&
     &4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(&
     &4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4)&
     & - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) &
     &- 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(&
     &1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta *&
     & yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta *&
     & yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta *&
     & yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) *&
     & xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi&
     & - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1)&
     & * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) *&
     & yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0&
     &.4D1)))
      k(8,3) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &/ 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * x&
     &x(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta &
     &* yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta &
     &* yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta &
     &* yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3)&
     & * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * &
     &xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi &
     &- xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + &
     &xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (&
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (-0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(&
     &1) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(&
     &3) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) *&
     & eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) *&
     & eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) *&
     & yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy&
     &(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2)&
     & * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2)&
     & * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) *&
     & yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(8,4) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0&
     &.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) +&
     & eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(&
     &1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(&
     &2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(&
     &4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1&
     &) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) *&
     & yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy&
     &(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) *&
     & yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * y&
     &y(3)) * (0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 &
     &+ nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2&
     &) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - &
     &eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + &
     &xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta &
     &* yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta &
     &* yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2&
     &) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) *&
     & xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi&
     & + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy&
     &(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1&
     &) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi&
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
     &* E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 &
     &* (xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(&
     &2) + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) -&
     & xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) +&
     & xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) +&
     & xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx&
     &(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4)&
     & * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4)&
     & + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) -&
     & xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - &
     &0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1)&
     & + xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * y&
     &y(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * y&
     &y(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * y&
     &y(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * x&
     &i + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi -&
     & xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) *&
     & yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * y&
     &y(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4&
     &D1)))
      k(8,5) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) - 0.2D1 * (xi * xx(1) - xi * xx(2) + xi * xx&
     &(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4)) / (xx(1) * eta *&
     & yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta *&
     & yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta *&
     & yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) &
     &* xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * x&
     &i - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi -&
     & xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + x&
     &x(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0&
     &.1D1 / 0.4D1 + xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * &
     &nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - et&
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
     & (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * y&
     &y(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) &
     &/ (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3)&
     & + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2)&
     & + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi &
     &- xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + x&
     &x(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4&
     &) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - x&
     &x(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(&
     &4) * yy(3)) * (0.1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1)&
     & - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3)&
     & + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1&
     & - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy(1&
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
     & * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(8,6) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (0.1D1 / 0.4D1 + xi / 0.4D1)) * E * (0.1D1 - nu) / (0.1D1 +&
     & nu) / (0.1D1 - 0.2D1 * nu) * (-0.2D1 * (eta * yy(1) - eta * yy(2)&
     & + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (x&
     &x(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + x&
     &x(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + x&
     &x(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx&
     &(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3)&
     & * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * &
     &yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2)&
     & * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) *&
     & yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - e&
     &ta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + x&
     &x(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta *&
     & yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta *&
     & yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2)&
     & * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * &
     &xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi &
     &+ xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(&
     &1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1)&
     & + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (0.2D1 * (xi &
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
     &E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * &
     &(xi * yy(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2)&
     & + yy(3) - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - x&
     &x(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + x&
     &x(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + x&
     &x(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2&
     &) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) *&
     & yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) +&
     & xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - x&
     &x(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.&
     &2D1 * (xi * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) +&
     & xx(2) + xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(&
     &4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(&
     &1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(&
     &2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi &
     &+ xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - x&
     &x(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * y&
     &y(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(&
     &4) - xx(4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1&
     &)))
      k(8,7) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
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
     &0.1D1 / 0.4D1 - xi / 0.4D1)) * E / (0.1D1 + nu) / (0.1D1 - 0.2D1 *&
     & nu) * nu * (-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3) - e&
     &ta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * yy(3)&
     & - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * yy(4)&
     & - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * yy(1)&
     & - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) * xi -&
     & xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi - xx&
     &(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - xx(1)&
     & * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx(3) *&
     & yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + eta *&
     & xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1) * e&
     &ta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * e&
     &ta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * e&
     &ta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy&
     &(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1)&
     & * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * &
     &xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3)&
     & + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)) + (-0.2D1 * (eta * yy(1) - eta * &
     &yy(2) + eta * yy(3) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4))&
     & / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3&
     &) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2&
     &) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi&
     & - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + &
     &xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(&
     &4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - &
     &xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx&
     &(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1&
     &) - eta * xx(2) + eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3&
     &) + xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * &
     &eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * &
     &eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * &
     &yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(&
     &4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) &
     &* xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) &
     &* yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * &
     &yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) * E * (0.1&
     &D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu) * (0.2D1 * (xi * yy&
     &(1) - xi * yy(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3)&
     & - yy(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * e&
     &ta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * e&
     &ta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * y&
     &y(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4&
     &) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) *&
     & xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) *&
     & yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * y&
     &y(1) + xx(4) * yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (x&
     &i * xx(1) - xi * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) +&
     & xx(3) - xx(4)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(&
     &2) * eta * yy(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(&
     &3) * eta * yy(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(&
     &1) * yy(2) * xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) &
     &* yy(4) * xi + xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * y&
     &y(2) * xi + xx(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + x&
     &x(2) * yy(1) - xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(&
     &4) * yy(1) + xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)))
      k(8,8) = (-xx(3) * eta * yy(1) / 0.8D1 - xx(2) * eta * yy(3) / 0.8&
     &D1 - xx(3) * yy(4) * xi / 0.8D1 + xx(4) * eta * yy(1) / 0.8D1 - xx&
     &(1) * eta * yy(4) / 0.8D1 - xx(4) * eta * yy(2) / 0.8D1 + xx(1) * &
     &yy(2) * xi / 0.8D1 + xx(3) * eta * yy(2) / 0.8D1 + xx(2) * eta * y&
     &y(4) / 0.8D1 + xx(3) * yy(1) * xi / 0.8D1 + xx(1) * eta * yy(3) / &
     &0.8D1 + xx(2) * yy(4) * xi / 0.8D1 - xx(2) * yy(1) * xi / 0.8D1 - &
     &xx(1) * yy(3) * xi / 0.8D1 - xx(4) * yy(2) * xi / 0.8D1 + xx(4) * &
     &yy(3) * xi / 0.8D1 - xx(1) * yy(2) / 0.8D1 + xx(1) * yy(4) / 0.8D1&
     & + xx(2) * yy(1) / 0.8D1 - xx(2) * yy(3) / 0.8D1 + xx(3) * yy(2) /&
     & 0.8D1 - xx(3) * yy(4) / 0.8D1 - xx(4) * yy(1) / 0.8D1 + xx(4) * y&
     &y(3) / 0.8D1) * ((-0.2D1 * (eta * yy(1) - eta * yy(2) + eta * yy(3&
     &) - eta * yy(4) - yy(1) - yy(2) + yy(3) + yy(4)) / (xx(1) * eta * &
     &yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2) * eta * &
     &yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4) * eta * &
     &yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1) * yy(3) *&
     & xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * yy(1) * xi&
     & - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(3) * xi - &
     &xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * yy(3) + xx&
     &(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy(3)) * (0.&
     &1D1 / 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * xx(1) - eta * xx(2) + &
     &eta * xx(3) - eta * xx(4) - xx(1) - xx(2) + xx(3) + xx(4)) / (xx(1&
     &) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + xx(2&
     &) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + xx(4&
     &) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - xx(1)&
     & * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3) * &
     &yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) * yy(&
     &3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2) * &
     &yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) * yy&
     &(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 - nu) / (0&
     &.1D1 + nu) / (0.1D1 - 0.2D1 * nu) + (0.2D1 * (xi * yy(1) - xi * yy&
     &(2) + xi * yy(3) - xi * yy(4) - yy(1) + yy(2) + yy(3) - yy(4)) / (&
     &xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy(3) + &
     &xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy(2) + &
     &xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * xi - x&
     &x(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi + xx(3&
     &) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + xx(4) *&
     & yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) - xx(2&
     &) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + xx(4) &
     &* yy(3)) * (0.1D1 / 0.4D1 - eta / 0.4D1) - 0.2D1 * (xi * xx(1) - x&
     &i * xx(2) + xi * xx(3) - xi * xx(4) - xx(1) + xx(2) + xx(3) - xx(4&
     &)) / (xx(1) * eta * yy(3) - xx(1) * eta * yy(4) - xx(2) * eta * yy&
     &(3) + xx(2) * eta * yy(4) - xx(3) * eta * yy(1) + xx(3) * eta * yy&
     &(2) + xx(4) * eta * yy(1) - xx(4) * eta * yy(2) + xx(1) * yy(2) * &
     &xi - xx(1) * yy(3) * xi - xx(2) * yy(1) * xi + xx(2) * yy(4) * xi &
     &+ xx(3) * yy(1) * xi - xx(3) * yy(4) * xi - xx(4) * yy(2) * xi + x&
     &x(4) * yy(3) * xi - xx(1) * yy(2) + xx(1) * yy(4) + xx(2) * yy(1) &
     &- xx(2) * yy(3) + xx(3) * yy(2) - xx(3) * yy(4) - xx(4) * yy(1) + &
     &xx(4) * yy(3)) * (-0.1D1 / 0.4D1 - xi / 0.4D1)) ** 2 * E * (0.1D1 &
     &- nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu))

end function k



end subroutine generateesm


