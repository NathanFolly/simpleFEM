! Subroutine that permits to recover the stress from the displacements that have 
! been calculated.
! Uses the more rudimentary approach (1) according to:
! http://www.colorado.edu/!engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf 
subroutine recoverstress(element, property)
use mytypes
implicit none

type(elementtype),intent(inout) :: element
type(propertytype),intent(inout) :: propety
real, allocatable:: ux(:)=0, uy(:)=0
real, allocatable:: stressvector(:)=0

real :: eta=0, xi=0

real, dimension(3,3) :: Dee=0
real:: E, nu


!construct pointer vectors from : 4-element ux vector and 4-element uy vector to -element u vector
allocate(ux(size(element%node)))
allocate(uy(size(element%node)))

do i = 1, size(ux)
     ux(i)= element%node(i)%state%ux
     uy(i)= element%node(i)%state%uy
end do







print *, ux

E= property%E
nu = property%nu
! The strain-stress relation matrix:
      Dee(1,1) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      Dee(1,2) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      Dee(2,1) = E / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu) * nu
      Dee(2,2) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.1D1 - 0.2D1 * nu)
      Dee(3,3) = E * (0.1D1 - nu) / (0.1D1 + nu) / (0.2D1 - 0.2D1 * nu)

! Calculating the stresstensor( in this case a vector with the form sigmax1, sigmay1, tauxy1, sigmax2, sigmay2, ...):


stressvector= 



end subroutine recoverstress


