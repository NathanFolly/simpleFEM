! Subroutine that permits to recover the stress from the displacements that have 
! been calculated.
! Uses the more rudimentary approach (1) according to:
! http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf 
subroutine recoverstress(element, property, displacementvector, stressvector)
use mytypes
implicit none

type(elementtype), target, intent(inout) :: element
type(propertytype), target, intent(inout) :: propety
real, pointer :: ux(:), uy(:)


real, dimension(3,3) :: Dee=0
real, pointer :: E, nu


!construct pointer vectors from : 4-element ux vector and 4-element uy vector to -element u vector

ux => element 

E=> property%E
nu => property%nu
! The strain-stress relation matrix:
      Dee(1,1) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      Dee(1,2) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      Dee(2,1) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      Dee(2,2) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      Dee(3,3) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu)



function strainmatrix(ux, uy, E, nu, xi, eta)

	real, dimension(8,8)

strainmatrix(1,1) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 / 0.4&
     &D1 + eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) &
     &- xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(&
     &3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(&
     &4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(&
     &1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3)&
     & - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - &
     &xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(&
     &1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3)&
     & * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D&
     &1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(1,3) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 / 0.4&
     &D1 - eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) &
     &- xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(&
     &3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(&
     &4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(&
     &1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3)&
     & - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - &
     &xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(&
     &1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3)&
     & * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1&
     & / 0.4D1 - xi / 0.4D1)
      strainmatrix(1,5) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 / 0.4D&
     &1 + eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) -&
     & xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(3&
     &) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4&
     &) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1&
     &) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) &
     &- xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - x&
     &i * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1&
     &) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) &
     &* uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)
      strainmatrix(1,7) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 / 0.4D&
     &1 - eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) -&
     & xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(3&
     &) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4&
     &) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1&
     &) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) &
     &- xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - x&
     &i * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1&
     &) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) &
     &* uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1&
     & / 0.4D1 - xi / 0.4D1)
      strainmatrix(2,2) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 &
     &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta &
     &* ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux&
     &(1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux&
     &(2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux&
     &(4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1&
     &) * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) *&
     & uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy&
     &(3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3&
     &) + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3))&
     & * (-0.1D1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(2,4) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta &
     &* ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux&
     &(1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux&
     &(2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux&
     &(4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1&
     &) * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) *&
     & uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy&
     &(3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3&
     &) + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3))&
     & * (0.1D1 / 0.4D1 - xi / 0.4D1)
      strainmatrix(2,6) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta *&
     & ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux(&
     &1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(&
     &2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(&
     &4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1)&
     & * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * &
     &uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(&
     &3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3)&
     & + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) &
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(2,8) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta *&
     & ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux(&
     &1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(&
     &2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(&
     &4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1)&
     & * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * &
     &uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(&
     &3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3)&
     & + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)
      strainmatrix(3,1) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 &
     &/ 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta &
     &* ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux&
     &(1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux&
     &(2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux&
     &(4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1&
     &) * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) *&
     & uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy&
     &(3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3&
     &) + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3))&
     & * (-0.1D1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(3,2) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 / 0.4&
     &D1 + eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) &
     &- xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(&
     &3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(&
     &4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(&
     &1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3)&
     & - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - &
     &xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(&
     &1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3)&
     & * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D&
     &1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(3,3) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 &
     &/ 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta &
     &* ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux&
     &(1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux&
     &(2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux&
     &(4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1&
     &) * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) *&
     & uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy&
     &(3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3&
     &) + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3))&
     & * (0.1D1 / 0.4D1 - xi / 0.4D1)
      strainmatrix(3,4) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1 / 0.4&
     &D1 - eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) &
     &- xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(&
     &3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(&
     &4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(&
     &1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3)&
     & - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - &
     &xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(&
     &1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3)&
     & * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1&
     & / 0.4D1 - xi / 0.4D1)
      strainmatrix(3,5) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 /&
     & 0.4D1 + eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta *&
     & ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux(&
     &1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(&
     &2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(&
     &4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1)&
     & * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * &
     &uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(&
     &3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3)&
     & + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) &
     &* (0.1D1 / 0.4D1 + xi / 0.4D1)
      strainmatrix(3,6) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 / 0.4D&
     &1 + eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) -&
     & xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(3&
     &) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4&
     &) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1&
     &) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) &
     &- xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - x&
     &i * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1&
     &) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) &
     &* uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 &
     &/ 0.4D1 + xi / 0.4D1)
      strainmatrix(3,7) = -0.2D1 * (eta * uy(1) - eta * uy(2) + eta * uy(3) - e&
     &ta * uy(4) - uy(1) - uy(2) + uy(3) + uy(4)) / (eta * ux(1) * uy(3)&
     & - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4)&
     & - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1)&
     & - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) -&
     & xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi&
     & * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1)&
     & * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) *&
     & uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 /&
     & 0.4D1 - eta / 0.4D1) + 0.2D1 * (eta * ux(1) - eta * ux(2) + eta *&
     & ux(3) - eta * ux(4) - ux(1) - ux(2) + ux(3) + ux(4)) / (eta * ux(&
     &1) * uy(3) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(&
     &2) * uy(4) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(&
     &4) * uy(1) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1)&
     & * uy(3) - xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * &
     &uy(1) - xi * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(&
     &3) - ux(1) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3)&
     & + ux(3) * uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) &
     &* (-0.1D1 / 0.4D1 - xi / 0.4D1)
      strainmatrix(3,8) = 0.2D1 * (xi * uy(1) - xi * uy(2) + xi * uy(3) - xi * &
     &uy(4) - uy(1) + uy(2) + uy(3) - uy(4)) / (eta * ux(1) * uy(3) - et&
     &a * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4) - et&
     &a * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1) - et&
     &a * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) - xi *&
     & ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - xi * ux&
     &(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1) * uy&
     &(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) * uy(2&
     &) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (0.1D1 / 0.4D&
     &1 - eta / 0.4D1) - 0.2D1 * (xi * ux(1) - xi * ux(2) + xi * ux(3) -&
     & xi * ux(4) - ux(1) + ux(2) + ux(3) - ux(4)) / (eta * ux(1) * uy(3&
     &) - eta * ux(1) * uy(4) - eta * ux(2) * uy(3) + eta * ux(2) * uy(4&
     &) - eta * ux(3) * uy(1) + eta * ux(3) * uy(2) + eta * ux(4) * uy(1&
     &) - eta * ux(4) * uy(2) + xi * ux(1) * uy(2) - xi * ux(1) * uy(3) &
     &- xi * ux(2) * uy(1) + xi * ux(2) * uy(4) + xi * ux(3) * uy(1) - x&
     &i * ux(3) * uy(4) - xi * ux(4) * uy(2) + xi * ux(4) * uy(3) - ux(1&
     &) * uy(2) + ux(1) * uy(4) + ux(2) * uy(1) - ux(2) * uy(3) + ux(3) &
     &* uy(2) - ux(3) * uy(4) - ux(4) * uy(1) + ux(4) * uy(3)) * (-0.1D1&
     & / 0.4D1 - xi / 0.4D1)
end function strainmatrix

end subroutine recoverstress


