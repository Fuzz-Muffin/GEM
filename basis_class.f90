module basis_class
	use para
	use grid_class
	implicit none
	private

	! This is the Gaussian function object for the basis
	type, public :: gaussian
		real(kind=range) :: N, nu
		real(kind=range), allocatable,dimension(:) :: f
		integer :: l, minf, maxf
	end type Gaussian

	type, public :: complex_gaussian
		real(kind=range) :: a, nu
		complex(kind=range) :: eta, N, N_c, C1, C2
		real(kind=range), allocatable,dimension(:) :: f
		integer :: l, minf, maxf
	end type complex_gaussian

	! Regular Gaussian basis object
	type, public :: basis
		type(gaussian), allocatable,dimension(:) :: G
!		real(kind=range), allocatable,dimension(:,:) :: H, T, V, S 
	integer :: nfancy
	contains
		procedure 	:: new_basis, get_nu, get_N_real, get_ell_r,& 
						calc_G, print_nu, print_N, print_f_reg
		generic 	:: new => new_basis 
		generic		:: calc => calc_G
		generic		:: get_N => get_N_real
		generic 	:: get_ell => get_ell_r
		generic 	:: print_f => print_f_reg
	end type basis

	! Complex-range Gaussian basis object
	type, public :: basis_c
		type(complex_gaussian), allocatable,dimension(:) :: G
!		real(kind=range), allocatable,dimension(:,:) :: H, T, V, S 
		integer :: nfancy
	contains
		procedure 	:: new_basis_complex, get_eta, get_N_c, &
						get_N_complex, get_C1, get_C2, get_a, &
						get_ell_c, calc_G_c 
		procedure 	:: print_eta, print_N_c, print_C1, print_C2, & 
						print_f_complex, print_a
		generic		:: new => new_basis_complex
	   	generic 	:: calc => calc_G_c	
		generic		:: get_N => get_N_complex
		generic 	:: print_f => print_f_complex
		generic 	:: get_ell => get_ell_c
	end type basis_c

contains

!-----------------------------------------------!
! regular basis subroutines 					!
!-----------------------------------------------!
	subroutine new_basis(self, l, nfancy, r1, rN, grid)
		class(basis), intent(inout) :: self
		type(GridObject), intent(in) :: grid
		integer, intent(in) :: l, nfancy
		integer :: i, j, alstat
		real(kind=range), intent(in) :: r1, rN
		real(kind=range) :: tmp

		self%nfancy = nfancy
		if (.not. allocated(self%G)) then
			allocate(self%G(nfancy),STAT = alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate basis'
		end if
		do i=1,nfancy
			allocate (self%G(i)%f(grid%nr))
			self%G(i)%l = l
			self%G(i)%nu= (r1*(rN/r1)**((real(i-1))/(real(nfancy-1))))**(-2)
			self%G(i)%N	= sqrt(((2.0d0**(l+2))*(2.0d0*self%G(i)%nu)**&
							(dble(l)+1.5d0))/(sqrt(pi)*dble(doublefact(2*l+1))))
			do j=1,grid%nr
				tmp = self%G(i)%nu*grid%rgrid(j)**2 
				if (tmp > exptol) then
					self%G(i)%f(j)=0.0
				else
					self%G(i)%f(j)= self%G(i)%N*(grid%rgrid(j)**(self%G(i)%l))*exp(-tmp)
				end if 
			end do
			call minmaxi(self%G(i)%f, grid, self%G(i)%minf, self%G(i)%maxf)
		end do
	end subroutine new_basis

	! N getter
	subroutine get_N_real(self, i, N)
		class(basis), intent(in) :: self
		integer, intent(in) :: i
		real(kind=range), intent(out) :: N
		N=self%G(i)%N
	end subroutine get_N_real
	
	! nu getter
   	subroutine get_nu(self, i, nu)
		class(basis), intent(in) :: self
		integer, intent(in) :: i
		real(kind=range), intent(out) :: nu
		nu=self%G(i)%nu
	end subroutine get_nu

	! l getter
	subroutine get_ell_r(self,i, l)
		class(basis), intent(in) :: self
		integer :: i
		integer, intent(out) :: l
		l=self%G(i)%l
	end subroutine get_ell_r
	
	! Calculate value of the ith G function at radius r
	subroutine calc_G(self, i, r, res)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		real(kind=range) :: tmp
		real(kind=range), intent(in) :: r
		real(kind=range), intent(out) :: res
		tmp=self%G(i)%nu*r*r
		if (tmp>exptol) then
			res=0.0
		else
			res = real(self%G(i)%N)*(r**(self%G(i)%l))*exp(-tmp)
		end if
	end subroutine calc_g

	! Print nu value
	subroutine print_nu(self, i)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'nu_',i,'=',self%G(i)%nu
	end subroutine print_nu

	! Print N value
	subroutine print_N(self, i)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'N_',i,'=', self%G(i)%N
	end subroutine print_N

	! Print the ith function on the grid
	subroutine print_f_reg(self, i, grid)
		class(basis), intent(in) :: self
		type(GridObject) :: grid
		integer, intent(in) :: i
		integer :: j, N
		N = grid%nr
		print*,'(a,I3.3,a)', 'r 		G(',i,')'
		do j=1, N
			print*, grid%rgrid(j), '', self%G(i)%f(j)
		end do
	end subroutine print_f_reg
!-----------------------------------------------!

!-----------------------------------------------!
! complex basis subroutines 					!
!-----------------------------------------------!
! basis goes: cos, cos, cos, ... ,sin,sin,sin.

	subroutine new_basis_complex(self, l, nfancy, r1, rN, a, grid)
		class(basis_c), intent(inout) :: self
		type(GridObject), intent(in) :: grid
		integer, intent(in) :: l, nfancy
		integer :: i,j, alstat, length
		real(kind=range), intent(in) :: r1, rN, a
		real(kind=range) :: p, locut
		complex(kind=range) :: tmp, tmp_a
		
		self%nfancy = nfancy
		if (.not. allocated(self%G)) then
			allocate(self%G(nfancy),STAT = alstat)
				if (alstat/=0) STOP 'Insufficient space to allocate basis'
		end if

		p =-(dble(l)+1.5d0)

!-------------------------------------!
! This is for one at a time functions !
!-------------------------------------!
		do i=1, nfancy/2
			self%G(i)%l = l
			self%G(i)%a = a
			self%G(i)%nu 	= (r1*(rN/r1)**((real(i-1))/(real(nfancy/2-1))))**(-2)
			self%G(i)%eta 	= self%G(i)%nu+a*self%G(i)%nu*imag
			self%G(i)%N		= sqrt(cmplx(((2.0d0**(l+2))*(2.0d0*self%G(i)%eta)**(real(l)+1.5d0))&
								/(sqrt(pi)*real(doublefact(2*l+1))),kind=8))
			self%G(i)%N_c  	= (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
								((2.0d0*self%G(i)%eta)**p +(2.0d0*conjg(self%G(i)%eta))**p &
								+2.0d0*(2.0d0*self%G(i)%nu)**p))**(-0.5d0) 
			self%G(i)%C1	=self%G(i)%N_c/(2.d0*self%G(i)%N)
			self%G(i)%C2	=self%G(i)%N_c/(2.d0*conjg(self%G(i)%N))
			allocate(self%G(i)%f(grid%nr))
			do j=1, grid%nr
				tmp=self%G(i)%nu*grid%rgrid(j)*grid%rgrid(j); locut=1.0d-300
				tmp_a=self%G(i)%eta*grid%rgrid(j)*grid%rgrid(j)
				if (abs(tmp)>exptol .or. abs(tmp_a)>exptol) then
					self%G(i)%f(j)=0.d0
				else
					self%G(i)%f(j)=real(self%G(i)%N*((grid%rgrid(j)**self%G(i)%l)/2.0)*(exp(-tmp_a)+exp(-conjg(tmp_a))),kind=range)
				end if
			end do
			call minmaxi(self%G(i)%f, grid, self%G(i)%minf, self%G(i)%maxf)
		end do
			self%G(nfancy/2+1:nfancy)%nu=self%G(1:nfancy/2)%nu
		do i=nfancy/2+1,nfancy
			self%G(i)%eta 	= self%G(i)%nu+a*self%G(i)%nu*imag
			self%G(i)%N		= sqrt(cmplx(((2.0d0**(l+2))*(2.0d0*self%G(i)%eta)**(real(l)+1.5d0))&
								/(sqrt(pi)*real(doublefact(2*l+1))),kind=8))
			self%G(i)%N_c 	= (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
								(-(2.0d0*self%G(i)%eta)**p -(2.0d0*conjg(self%G(i)%eta))**p &
								+2.0d0*(2.0d0*self%G(i)%nu)**p))**(-0.5d0) 
			self%G(i)%C1	= (-imag*self%G(i)%N_c)/(2.d0*self%G(i)%N)
			self%G(i)%C2	= -(-imag*self%G(i)%N_c)/(2.d0*conjg(self%G(i)%N))
			do j=1, grid%nr
				tmp=self%G(i)%nu*grid%rgrid(j)*grid%rgrid(j); locut=1.0d-300
				tmp_a=self%G(i)%eta*grid%rgrid(j)*grid%rgrid(j)
				if (abs(tmp)>exptol .or. abs(tmp_a)>exptol) then
					self%G(i)%f(j)=0.d0
				else
					self%G(i)%f(j)=real(self%G(i)%N*((grid%rgrid(j)**self%G(i)%l)/(2.0*imag))*(exp(-tmp_a)-exp(-conjg(tmp_a))),kind=range)
				end if
			end do
			call minmaxi(self%G(i)%f, grid, self%G(i)%minf, self%G(i)%maxf)
		end do
	end subroutine new_basis_complex

	subroutine get_N_c(self, i, N_c)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: N_c
		N_c=self%G(i)%N
	end subroutine get_N_c

	subroutine get_a(self, i, a)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		real(kind=range), intent(out) :: a
		a=self%G(i)%a
	end subroutine get_a

	subroutine get_eta(self, i, eta)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: eta
		eta=self%G(i)%eta
	end subroutine get_eta

	subroutine get_N_complex(self, i, N)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: N
		N=self%G(i)%N
	end subroutine get_N_complex

	subroutine get_C1(self, i, C)
		class(basis_c), intent(in) :: self
		integer :: i
		complex(kind=range), intent(out) :: C
		C=self%G(i)%C1
	end subroutine get_C1

	subroutine get_C2(self, i, C)
		class(basis_c), intent(in) :: self
		integer :: i
		complex(kind=range), intent(out) :: C
		C=self%G(i)%C2
	end subroutine get_C2
	
	subroutine get_ell_c(self,i, l)
		class(basis_c), intent(in) :: self
		integer :: i
		integer, intent(out) :: l
		l=self%G(i)%l
	end subroutine get_ell_c
	
	! Calculate value of complex G_i function at radius r
	subroutine calc_G_c(self, i, r, res)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		integer :: length 
		real(kind=range), intent(in) :: r
	   	real(kind=range) ::	tmp, locut
		complex(kind=range) :: tmp_a
		complex(kind=range), intent(out) :: res	
		tmp=self%G(i)%nu*r*r; locut=1.0d-300; length=size(self%G(:))
		tmp_a=self%G(i)%eta*r*r
		if (abs(tmp)>exptol .or. abs(tmp_a)>exptol) then
			res=0
		else
			if (i<=(length/2)) then
				res = self%G(i)%N*((r**self%G(i)%l)/2.0)*(exp(-tmp_a)+exp(-conjg(tmp_a)))
			else if (i>=(length/2+1)) then 
				res = self%G(i)%N*((r**self%G(i)%l)/(2.0*imag))*(exp(-tmp_a)-exp(-conjg(tmp_a)))
			end if
		end if
		if (abs(dble(res))<locut) then 
			res=(0.d0,0.d0)
		end if
	end subroutine calc_g_c

	subroutine print_a(self, i)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		print*, 'alpha=', self%G(i)%a
	end subroutine print_a	

	subroutine print_eta(self, i)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'eta_',i,'=',real(self%G(i)%eta),'i*',aimag(self%G(i)%eta)
	end subroutine print_eta

	subroutine print_N_c(self, i)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'N_c_',i,'=',real(self%G(i)%N),'i*',aimag(self%G(i)%N)
	end subroutine print_N_c

	subroutine print_C1(self, i)
		class(basis_c), intent(in) :: self
		integer :: i
		print*, 'C1_',i,'=', self%G(i)%C1
	end subroutine print_C1

	subroutine print_C2(self, i)
		class(basis_c), intent(in) :: self
		integer :: i
		print*, 'C2_',i,'=', self%G(i)%C2
	end subroutine print_C2
	
	! Print the ith function on the grid
	subroutine print_f_complex(self, i, grid)
		class(basis_c), intent(in) :: self
		type(GridObject) :: grid
		integer, intent(in) :: i
		integer :: j, N
		N = grid%nr
		print*,'(a,I3.3,a)', 'r 		G(',i,')'
		do j=1, N
			print*, grid%rgrid(j), '', self%G(i)%f(j)
		end do
	end subroutine print_f_complex

end module basis_class
