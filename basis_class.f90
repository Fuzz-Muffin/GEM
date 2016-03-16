module basis_class
	use para
	implicit none
	private

	! This is the Gaussian function object for the basis
	type, public :: basis
		real(kind=range), allocatable,dimension(:) :: nu
		complex(kind=range), allocatable,dimension(:) :: N
		integer :: l
	contains
		procedure 	:: new_basis
		procedure 	:: get_nu
		procedure 	:: get_N_real
		procedure 	:: get_ell
		procedure 	:: calc_G
		procedure	:: print_nu, print_N
		generic 	:: new => new_basis 
		generic		:: calc => calc_G
		generic		:: get_N => get_N_real
	end type basis

	! Complex gaussian function object
	type, extends(basis), public :: basis_c
		complex(kind=range), allocatable,dimension(:) :: eta, N_c, C1, C2
		real(kind=range) :: a
	contains
		procedure 	:: new_basis_complex
		procedure 	:: get_eta
		procedure 	:: get_N_c
		procedure 	:: get_N_complex
		procedure 	:: get_C1
		procedure 	:: get_C2
		procedure 	:: get_a
		procedure 	:: calc_G_c
		procedure 	:: print_eta, print_N_c, print_C1, print_C2
		generic		:: new => new_basis_complex
	   	generic 	:: calc => calc_G_c	
		generic		:: get_N => get_N_complex
	end type basis_c

contains

!-----------------------------------------------!
! regular basis subroutines 					!
!-----------------------------------------------!
	subroutine new_basis(self, l, nfancy, r1, rN)
		class(basis), intent(inout) :: self
		integer, intent(in) :: l, nfancy
		integer :: i, alstat
		real(kind=range), intent(in) :: r1, rN	
		if (.not. allocated(self%nu)) then
			allocate(self%nu(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%N)) then
			allocate(self%N(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		self%l 	= l
		do i=1,nfancy
			self%nu(i)	= (r1*(rN/r1)**((real(i-1))/(real(nfancy-1))))**(-2)
			self%N(i)= sqrt(((2.0d0**(l+2))*(2.0d0*self%nu(i))**(dble(l)+1.5d0))/(sqrt(pi)*dble(doublefact(2*l+1))))
		end do
	end subroutine new_basis

	! N getter
	subroutine get_N_real(self, i, N)
		class(basis), intent(in) :: self
		integer, intent(in) :: i
		real(kind=range), intent(out) :: N
		N=real(self%N(i))
	end subroutine get_N_real
	
	! nu getter
   	subroutine get_nu(self, i, nu)
		class(basis), intent(in) :: self
		integer, intent(in) :: i
		real(kind=range), intent(out) :: nu
		nu=self%nu(i)
	end subroutine get_nu

	! l getter
	subroutine get_ell(self, l)
		class(basis), intent(in) :: self
		integer, intent(out) :: l
		l=self%l
	end subroutine get_ell
	
	! Calculate value of G function at radius r
	subroutine calc_G(self, i, r, res)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		real(kind=range) :: tmp
		real(kind=range), intent(in) :: r
		real(kind=range), intent(out) :: res
		tmp=self%nu(i)*r*r
		if (tmp>exptol) then
			res=0.0
		else
			res = real(self%N(i))*(r**(self%l))*exp(-tmp)
		end if
	end subroutine calc_g

	! Print nu value
	subroutine print_nu(self, i)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'nu_',i,'=',self%nu(i)
	end subroutine print_nu

	! Print N value
	subroutine print_N(self, i)
		class(basis), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'N_',i,'=',real(self%N(i))
	end subroutine print_N
!-----------------------------------------------!

!-----------------------------------------------!
! complex basis subroutines 					!
!-----------------------------------------------!
! basis goes cos, sin, cos, sin... etc

	subroutine new_basis_complex(self, l, nfancy, r1, rN, a)
		class(basis_c), intent(inout) :: self
		integer, intent(in) :: l, nfancy
		integer :: i,j, alstat
		real(kind=range), intent(in) :: r1, rN, a
		real(kind=range) :: p
		complex(kind=range) :: tmp	
		
		if (.not. allocated(self%nu)) then
			allocate(self%nu(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%N)) then
			allocate(self%N(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%eta)) then
			allocate(self%eta(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%N_c)) then
			allocate(self%N_c(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%C1)) then
			allocate(self%C1(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if
		if (.not. allocated(self%C2)) then
			allocate(self%C2(nfancy),STAT = alstat)
				if (alstat/=0) then
					print*, 'Insufficient space to allocate basis'
					STOP
				end if
		end if

		self%l 	= l
		self%a 	= a
		p		=-(dble(l)+1.5d0)

!-------------------------------------!
! This is for one at a time functions !
!-------------------------------------!
		do i=1, nfancy/2
			self%nu(i)	= (r1*(rN/r1)**((real(i-1))/(real(nfancy/2-1))))**(-2)
			self%eta(i) = self%nu(i)+a*self%nu(i)*imag
			self%N(i)	= sqrt(cmplx(((2.0d0**(l+2))*(2.0d0*self%eta(i))**(real(l)+1.5d0))/(sqrt(pi)*real(doublefact(2*l+1))),kind=8))
			self%N_c(i) = (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
				((2.0d0*self%eta(i))**p +(2.0d0*conjg(self%eta(i)))**p +2.0d0*(2.0d0*self%nu(i))**p))**(-0.5d0) 
			self%C1(i)	=self%N_c(i)/(2.d0*self%N(i))
			self%C2(i)	=self%N_c(i)/(2.d0*conjg(self%N(i)))
		end do
			self%nu(nfancy/2+1:nfancy)=self%nu(1:nfancy/2)
			!self%eta(nfancy/2+1:nfancy)=self%eta(1:nfancy/2)
			!self%N(nfancy/2+1:nfancy)=self%N(1:nfancy/2)
		do i=nfancy/2+1,nfancy
			self%eta(i) = self%nu(i)+a*self%nu(i)*imag
			self%N(i)	= sqrt(cmplx(((2.0d0**(l+2))*(2.0d0*self%eta(i))**(real(l)+1.5d0))/(sqrt(pi)*real(doublefact(2*l+1))),kind=8))
			self%N_c(i) = (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
				(-(2.0d0*self%eta(i))**p -(2.0d0*conjg(self%eta(i)))**p +2.0d0*(2.0d0*self%nu(i))**p))**(-0.5d0) 
			self%C1(i)	=(-imag*self%N_c(i))/(2.d0*self%N(i))
			self%C2(i)	=-(-imag*self%N_c(i))/(2.d0*conjg(self%N(i)))
		end do

!-----------------------------------!
! This is for alternating functions !
!-----------------------------------!
!		do i=1,nfancy/2
!			self%nu((i*2)-1)= (r1*(rN/r1)**((real(i-1))/(real((nfancy/2)-1))))**(-2)
!			self%nu(i*2)= self%nu((i*2)-1)
!		end do
!		do i=1, nfancy, 2
!			self%eta(i) = self%nu(i)+a*self%nu(i)*imag
!			self%eta(i+1) = self%eta(i)
!			self%N(i)	= sqrt(cmplx(((2.0d0**(l+2))*(2.0d0*self%eta(i))**(real(l)+1.5d0))/(sqrt(pi)*real(doublefact(2*l+1))),kind=8))
!			self%N(i+1)	= self%N(i)
!		end do
!		do i=1, nfancy
!			if (mod(i,2)/=0) then
!				self%N_c(i) = (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
!					((2.0d0*self%eta(i))**p +(2.0d0*conjg(self%eta(i)))**p +2.0d0*(2.0d0*self%nu(i))**p))**(-0.5d0) 
!				self%C1(i)	=self%N_c(i)/(2.d0*self%N(i))
!				self%C2(i)	=self%N_c(i)/(2.d0*conjg(self%N(i)))
!			else
!				self%N_c(i) = (((sqrt(pi)*doublefact(2*l+1))/((2.0d0**l)*16.0d0))*&
!					(-(2.0d0*self%eta(i))**p -(2.0d0*conjg(self%eta(i)))**p +2.0d0*(2.0d0*self%nu(i))**p))**(-0.5d0) 
!				self%C1(i)	=(-imag*self%N_c(i))/(2.d0*self%N(i))
!				self%C2(i)	=-(-imag*self%N_c(i))/(2.d0*conjg(self%N(i)))
!			end if
!		end do

	end subroutine new_basis_complex

	subroutine get_N_c(self, i, N_c)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: N_c
		N_c=self%N_c(i)
	end subroutine get_N_c

	subroutine get_a(self, a)
		class(basis_c), intent(in) :: self
		real(kind=range), intent(out) :: a
		a=self%a
	end subroutine get_a

	subroutine get_eta(self, i, eta)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: eta
		eta=self%eta(i)
	end subroutine get_eta

	subroutine get_N_complex(self, i, N)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i
		complex(kind=range), intent(out) :: N
		N=self%N(i)
	end subroutine get_N_complex

	subroutine get_C1(self, i, C)
		class(basis_c), intent(in) :: self
		integer :: i
		complex(kind=range), intent(out) :: C
		C=self%C1(i)
	end subroutine get_C1

	subroutine get_C2(self, i, C)
		class(basis_c), intent(in) :: self
		integer :: i
		complex(kind=range), intent(out) :: C
		C=self%C2(i)
	end subroutine get_C2
	
	! Calculate value of complex G_i function at radius r
	subroutine calc_G_c(self, i, r, res)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		integer :: length 
		real(kind=range), intent(in) :: r
	   	real(kind=range) ::	tmp, locut
		complex(kind=range) :: tmp_a
		complex(kind=range), intent(out) :: res	
		tmp=self%nu(i)*r*r; locut=1.0d-300; length=size(self%N_c)
		tmp_a=self%eta(i)*r*r
		if (abs(tmp)>exptol .or. abs(tmp_a)>exptol) then
			res=0
		else
			!if(mod(i,2)/=0) then
			if (i<=(length/2)) then
				res = self%N_c(i)*((r**self%l)/2.0)*(exp(-tmp_a)+exp(-conjg(tmp_a)))
				!res = self%N_c(i)*(r**self%l)*exp(-tmp)*cos(self%a*tmp)
			else if (i>=(length/2+1)) then 
				res = self%N_c(i)*((r**self%l)/(2.0*imag))*(exp(-tmp_a)-exp(-conjg(tmp_a)))
				!res = self%N_c(i)*(r**self%l)*exp(-tmp)*sin(self%a*tmp)
			end if
		end if
		if (abs(dble(res))<locut) then 
			res=(0.d0,0.d0)
		end if
	end subroutine calc_g_c

	subroutine print_eta(self, i)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'eta_',i,'=',real(self%eta(i)),aimag(self%eta(i)),'i'
	end subroutine print_eta

	subroutine print_N_c(self, i)
		class(basis_c), intent(in) :: self
		integer, intent(in) :: i 
		print*, 'N_c_',i,'=',real(self%N_c(i)),aimag(self%N_c(i)),'i'
	end subroutine print_N_c

	subroutine print_C1(self, i)
		class(basis_c), intent(in) :: self
		integer :: i
		print*, self%C1(i)
	end subroutine print_C1

	subroutine print_C2(self, i)
		class(basis_c), intent(in) :: self
		integer :: i
		print*, self%C2(i)
	end subroutine print_C2

end module basis_class
