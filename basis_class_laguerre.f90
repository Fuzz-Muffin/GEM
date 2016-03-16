!---------------------------------------------------------------!
! Large sections of this code was written by Liam Scarlett 
! of Curtin University, cheers mate! 	
!---------------------------------------------------------------!
!
module basis_class_laguerre
  use  para
  use  grid_class
  implicit none
  private

  type, public :: SturmianObject
    real*8 :: energy
	real*8 :: norm !The normalisation constant for each function 
    integer :: l  !angular momentum
    integer :: n !order of sturmian function. Redundant for one-electron states
    integer :: minf, maxf
    real*8, allocatable :: f(:) !the function
  end type SturmianObject

  type, public :: BasisObject
    type(SturmianObject), allocatable :: b(:) !the basis
    type(SturmianObject), allocatable :: chi(:) !the orthonormal basis
    real*8, allocatable :: overlap(:,:) !the overlap matrix
    real*8, allocatable :: K(:,:) !the kinetic energy matrix
    real*8, allocatable :: H(:,:) !the Hamiltonian matrix
    real*8, allocatable :: V(:,:) !the potential matrix
    integer :: N	
  contains
		procedure :: constructLaguerreBasis, constructPotential, &
			constructOverlap, constructHamiltonian, constructOneElecFunctions, basis_cleanup
  end type BasisObject

contains

  subroutine constructLaguerreBasis(self, N, l, alpha, grid)
    class(BasisObject), intent(inout) :: self
    integer, intent(in) :: N !number of Laguerre functions
    integer, intent(in) :: l !angular momentum
    real*8, intent(in) :: alpha !falloff
    type(GridObject), intent(in) :: grid
    integer :: i, j, k
    real*8 :: x

    self%N = N
    allocate (self%b(N), self%chi(N))
    do i=1, N
      allocate (self%b(i)%f(grid%nr))
	 self%b(i)%l = l
	 self%b(i)%norm = sqrt((alpha*gamma(real(i)))/(real(i+l)*gamma(real(i+2*l+1,kind=8))))
	 !self%b(i)%norm = sqrt((alpha*fact(i-1))/(real((i+l)*fact(i+2*l),kind=8)))
	 !print*, self%b(i)%norm
    enddo
    call laguerrePoly(2*l+1, 2*alpha, self, grid)
    do j=1, N
      self%b(j)%f(:) = self%b(j)%norm*self%b(j)%f(:)*&
		(2*alpha*grid%rgrid(:))**(l+1)*exp(-alpha*grid%rgrid(:))
	 call minmaxi(self%b(j)%f, grid, self%b(j)%minf, self%b(j)%maxf)
    enddo
  end subroutine constructLaguerreBasis

  subroutine laguerrePoly(m, lambda, self, grid)
	implicit none
    !Populates the basis with laguerre polynomials
    integer, intent(in) :: m !parameter of Laguere polynomial
    real*8, intent(in) :: lambda
    type(GridObject), intent(in) :: grid
    type(BasisObject), intent(inout) :: self
    integer :: i, n
    real*8 :: r, x, L0, L1
    self%b(1)%f(:) = 1 !zeroth order laguerre poly is 1
    self%b(1)%n = 0 !set degree of poly
    if (size(self%b) .gt. 1) then !if more than one function needed
      do i=1, grid%nr !iterate over grid

       r = grid%rgrid(i)
	   x = lambda*r
	   self%b(2)%f(i) = 1+m-x
	   self%b(2)%n = 1 !set degree of poly

	   do n=2, size(self%b)-1
          !using recurrance relation:
	    self%b(n+1)%f(i)=(((2*n+1+m)-x)*self%b(n)%f(i)-(n+m) &
			*self%b(n-1)%f(i))/dble(n+1)
		self%b(n)%n = n !set degree of poly
	   enddo
      enddo
    endif
  end subroutine laguerrePoly

  subroutine constructPotential(self, alpha, V, grid)
    class(BasisObject), intent(inout) :: self
    type(GridObject), intent(in) :: grid
	integer, intent(in) :: V
	real*8, intent(in) :: alpha
    integer :: i, j, k, N, minf, maxf

    N = self%N
    allocate(self%V(N,N))

    self%V = 0.0
	if (V==4 .or. V==2) then
		STOP 'These potentials are not implimented for the Laguerre basis'
	else if (V==3) then 
	    do i=1, N
			self%V(i,i)= -(1.0/(self%b(i)%norm**2))*(alpha/(real(i+self%b(i)%l)))
		enddo
	end if
  end subroutine constructPotential

  subroutine constructOverlap(self, l, alpha, grid)
    class(BasisObject), intent(inout) :: self
    integer, intent(in) :: l !angular momentum
    real*8, intent(in) :: alpha
    integer :: i, j, minf, maxf
    type(GridObject), intent(in) :: grid
    real*8 :: tmp, tol
    integer :: N

    N = self%N
	allocate(self%overlap(N,N)); self%overlap=0.0
	!print*, N
	do i=1,N
		self%overlap(i,i) = 1.0/(self%b(i)%norm**2)
	end do
	do i=1,N-1
		self%overlap(i,i+1) = -0.5*(1.0/(self%b(i)%norm*self%b(i+1)%norm))*&
			   	sqrt(1.0-(real(l*(l+1))/real((i+l)*(i+l+1))))
		self%overlap(i+1,i) = self%overlap(i,i+1)
	end do 
	!print*,'Overlap', self%overlap
  end subroutine constructOverlap

  subroutine constructHamiltonian(self, l, alpha, grid)
    class(BasisObject), intent(inout) :: self
    integer, intent(in) :: l !angular momentum
    type(GridObject), intent(in) :: grid
    real*8, intent(in) :: alpha
    integer :: i, j, minf, maxf, i1, i2, ri
    real*8 :: k, eta, Vsum, al
    integer*16, external :: gamma_ratio
    integer*16 :: tmp
    integer :: N

    N = self%N
	allocate(self%H(N,N), self%K(N,N))
	
	self%K=0.0; self%H=0.0
	do i=1,N
		self%K(i,i)= (1.0/(self%b(i)%norm**2))*alpha**2&
				-(((alpha**2)/2.0)*(1.0/(self%b(i)%norm**2)))
	end do
	do i=1, N-1
		self%K(i,i+1)= ((alpha**2)/4.0)*sqrt(1.0-(real(l*(l+1))/real((i+l)*(i+l+1))))*&
				(1.0/(self%b(i)%norm*self%b(i+1)%norm))
		self%k(i+1,i)= self%K(i,i+1)
	end do
	!print*, 'K',self%K
	!print*, 'V',self%V
    self%H = self%K + self%V
	!print*, 'H',self%H
  end subroutine constructHamiltonian

  subroutine constructOneElecFunctions(self, w, CI, nst, grid)
    class(BasisObject), intent(inout) :: self
    integer, intent(in) :: nst
    type(GridObject), intent(in) :: grid
    integer :: matz, ierr, n, m, i, i1, i2
    real*8, intent(in) :: w(0:nst-1), CI(0:nst-1,0:nst-1)
    real*8 :: temp(grid%nr), norm
    
    do i=0, nst
      allocate(self%chi(i)%f(grid%nr))
    enddo

    do n=0, nst-1
      temp = 0
      do m=0, nst-1
	   do i=self%b(m)%minf, self%b(m)%maxf
	     temp(i) = temp(i) + CI(m,n)*self%b(m)%f(i)
	   enddo
	 enddo
	 
	 call minmaxi(temp, grid, i1, i2)
      if(sign(1d0, temp(i1+1)) .lt. 0) then
        do i=i1, i2
	     !temp(i) = -temp(i)
	   enddo
      endif
	 self%chi(n)%f = temp
	 self%chi(n)%energy = w(n)
	 self%chi(n)%minf = i1
	 self%chi(n)%maxf = i2
    enddo
  end subroutine constructOneElecFunctions

  subroutine basis_cleanup(self)
    class(BasisObject) :: self

    deallocate(self%b)
    deallocate(self%chi)
    deallocate(self%H)
    deallocate(self%V)
    deallocate(self%K)
    deallocate(self%overlap)
  end subroutine basis_cleanup


end module basis_class_laguerre
