module HTVS_class
	use para
	use basis_class
	implicit none
	private
	
	! This object contains the H, T, V and S matrices 
	type, public :: HTVS
		real(kind=range), allocatable,dimension(:,:) :: H, T, V, S
	contains
		procedure 	:: new_HTVS, new_HTVS_complex, print_HTVS, destruct
		generic		:: new => new_HTVS, new_HTVS_complex
		final 		:: delete_HTVS
	end type HTVS

contains

!**********************************************************************
! HTVS constructor, regular Gaussians								  *
!**********************************************************************
	subroutine new_HTVS(self, mu, nfancy, Gbasis, Vstat)
		class(HTVS), intent(inout) :: self
		type(basis), intent(in) :: Gbasis
		logical :: verbose
		integer, intent(in) :: Vstat, nfancy
	   	integer :: alstat, i
		real(kind=range), intent(in) :: mu
		real(kind=range) :: testH
		
		! We now allocate space for all the required matrices
		if (allocated(self%H)) deallocate(self%H) 
		allocate(self%H(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate H'
		
		if (allocated(self%T)) deallocate(self%T) 
		allocate(self%T(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate T'
		
		if (allocated(self%V)) deallocate(self%V) 
		allocate(self%V(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate V'

		if (allocated(self%S)) deallocate(self%S) 
		allocate(self%S(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate S'
		
		! Let's build a HTVS object
		call new_T(self%T, nfancy, Gbasis, mu)
		call new_V(self%V, nfancy, Gbasis, mu, Vstat)
		call new_S(self%S, nfancy, Gbasis)
		
		self%H = self%T+self%V

		!check to make sure H is not totally 0
		testH=SUM(self%H)
		if (testH==0) STOP 'ERROR in HTVS construction, H=0'
	end subroutine new_HTVS

!**********************************************************************
! HTVS constructor, complex-range Gaussians							  *
!**********************************************************************
	subroutine new_HTVS_complex(self, mu, nfancy, Gbasis, Vstat)
		class(HTVS), intent(inout) :: self
		type(basis_c), intent(in) :: Gbasis
		logical :: verbose
		integer, intent(in) :: Vstat, nfancy
	   	integer :: alstat, i
		real(kind=range), intent(in) :: mu
		real(kind=range) :: testH
		
		! We now allocate space for all the required matrices
		if (allocated(self%H)) deallocate(self%H) 
		allocate(self%H(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate H'
		
		if (allocated(self%T)) deallocate(self%T) 
		allocate(self%T(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate T'
		
		if (allocated(self%V)) deallocate(self%V) 
		allocate(self%V(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate V'

		if (allocated(self%S)) deallocate(self%S) 
		allocate(self%S(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'Insufficient space to allocate S'
		
		! Let's build a HTVS object
		call new_T_c(self%T, nfancy, Gbasis, mu)
		call new_V_c(self%V, nfancy, Gbasis, mu, Vstat)
		call new_S_c(self%S, nfancy, Gbasis)
		
		self%H = self%T+self%V

		!check to make sure H is not totally 0
		testH=SUM(self%H)
		if (testH==0) STOP 'ERROR in HTVS construction, H=0'
	end subroutine new_HTVS_complex

!**********************************************************************
! Matrix elements, regular Gaussians								  *
!**********************************************************************

	subroutine new_T(self, nfancy, Gbasis, mu)
		type(basis) :: Gbasis
		integer, intent(in) :: nfancy
		integer :: i, j, l
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		real(kind=range), intent(in) :: mu
		real(kind=range) :: nui, nuj
		call Gbasis%get_ell(l)
		do i=1, nfancy
			do j=1, nfancy
				call Gbasis%get_nu(i, nui)
				call Gbasis%get_nu(j, nuj)
				self(i,j)= (hbar**2/mu)*(((2.0*real(l)+3.0)*nui*nuj)/(nui+nuj))&
				*((2.0*sqrt(nui*nuj))/(nui+nuj))**(real(l)+1.5)
			end do 
		end do 
	end subroutine new_T 

	subroutine new_V(self, nfancy, Gbasis, mu, Vstat)
		type(basis) :: Gbasis
		integer, intent(in) :: nfancy, Vstat
		integer :: i, j, l
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		real(kind=range), intent(in) :: mu
		real(kind=range) :: nui, nuj		
		! Using the Vstat parameter to select the desired potential
		call Gbasis%get_ell(l)
		select case (Vstat)
		case(1)
			self=0
		case(2)
			do i=1, nfancy
				do j=1, nfancy
					call Gbasis%get_nu(i, nui)
					call Gbasis%get_nu(j, nuj)
					self(i,j)= ((real(l)+1.5)/(nui &
						+nuj))*((2.0*sqrt(nui &	
						*nuj))/(nui+nuj))**(real(l)+1.5)
				end do 
			end do	
		case(3)
			do i=1, nfancy
				do j=1, nfancy 
					call Gbasis%get_nu(i, nui)
					call Gbasis%get_nu(j, nuj)
					self(i,j)= ((2.0)/(sqrt(pi)))* &
						((2.0**l*real(fact(l)))/(real(doublefact(2*l+1))))*&
						sqrt(nui+nuj)*((2.0*sqrt(nui*nuj))/(nui+nuj))**(real(l)+1.5)
				end do 
			end do
		case(4)
			do i=1, nfancy
				do j=1, nfancy
					call Gbasis%get_nu(i, nui)
					call Gbasis%get_nu(j, nuj)
					self(i,j)= ((2.0*sqrt(nui* &
						nuj))/(nui+nuj+mu))**(real(l)+1.5)
				end do 
			end do
		case default 
			self=0
		end select		
		! Attractive potential
		self=-self
	end subroutine new_V

	subroutine new_S(self, nfancy, Gbasis)
		type(basis) :: Gbasis
		integer :: i, j, l
		integer, intent(in) :: nfancy
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		real(kind=range) :: nui, nuj
		call Gbasis%get_ell(l)
		do i=1,nfancy
			do j=1,nfancy
				call Gbasis%get_nu(i, nui)
				call Gbasis%get_nu(j, nuj)
				self(i,j) = ((2.0*sqrt(nui*nuj))/(nui+nuj))**(real(l)+1.5)
			end do 
		end do
	end subroutine new_S

!**********************************************************************
! Matrix elements, complex-range Gaussians							  *
!**********************************************************************

	subroutine new_T_c(self, nfancy, Gbasis, mu)
		type(basis_c) :: Gbasis
		integer, intent(in) :: nfancy
		integer :: i, j, l
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		real(kind=range), intent(in) :: mu
		complex(kind=range) :: etai, etaj, tmp_a, tmp_b, tmp_c, tmp_d,&
				C1i, C2i, C1j, C2j 
		call Gbasis%get_ell(l)
		do i=1, nfancy
			call Gbasis%get_eta(i, etai)
			call Gbasis%get_C1(i, C1i)
			call Gbasis%get_C2(i, C2i)
			do j=1, nfancy
				call Gbasis%get_eta(j, etaj)
				call Gbasis%get_C1(j, C1j)
				call Gbasis%get_C2(j, C2j)

				tmp_a=conjg(C1i)*C1j*((hbar**2)/(mu))*(((2.0d0*dble(l)+3.0d0)*conjg(etai)*etaj)/(conjg(etai)+etaj))&
					*((2.d0*sqrt(conjg(etai)*etaj))/(conjg(etai)+etaj))**(dble(l)+1.5d0)
				tmp_b=conjg(C2i)*C1j*((hbar**2)/(mu))*(((2.0d0*dble(l)+3.0d0)*etai*etaj)/(etai+etaj))&
					*((2.d0*sqrt(etai*etaj))/(etai+etaj))**(dble(l)+1.5d0)
				tmp_c=conjg(C1i)*C2j*((hbar**2)/(mu))*(((2.0d0*dble(l)+3.0d0)*conjg(etai)*conjg(etaj))/(conjg(etai)+conjg(etaj)))&
					*((2.d0*sqrt(conjg(etai)*conjg(etaj)))/(conjg(etai)+conjg(etaj)))**(dble(l)+1.5d0)
				tmp_d=conjg(C2i)*C2j*((hbar**2)/(mu))*(((2.0d0*dble(l)+3.0d0)*conjg(etaj)*etai)/(conjg(etaj)+etai))&
					*((2.d0*sqrt(conjg(etaj)*etai))/(conjg(etaj)+etai))**(dble(l)+1.5d0)
				self(i,j)=real(tmp_a+tmp_b+tmp_c+tmp_d,kind=range)
			end do 
		end do 
	end subroutine new_T_c 

	subroutine new_V_c(self, nfancy, Gbasis, mu, Vstat)
		type(basis_c) :: Gbasis
		integer, intent(in) :: nfancy, Vstat
		integer :: i, j, l
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		real(kind=range), intent(in) :: mu
		complex(kind=range) :: etai, etaj, tmp_a, tmp_b, tmp_c, tmp_d,&
				C1i, C2i, C1j, C2j 
		! Using the Vstat parameter to select the desired potential
		call Gbasis%get_ell(l)
		select case (Vstat)
		case(1)
			self=0
		case(2)
			do i=1, nfancy
				call Gbasis%get_eta(i, etai)
				call Gbasis%get_C1(i, C1i)
				call Gbasis%get_C2(i, C2i)
				do j=1, nfancy
					call Gbasis%get_eta(j, etaj)
					call Gbasis%get_C1(j, C1j)
					call Gbasis%get_C2(j, C2j)
					tmp_a = conjg(C1i)*C1j*((real(l)+1.5)/(conjg(etai) &
								+etaj))*((2.d0*sqrt(conjg(etai) &	
								*etaj))/(conjg(etai)+etaj))**(real(l)+1.5)
					tmp_b = conjg(C2i)*C1j*((real(l)+1.5)/(etai &
								+etaj))*((2.d0*sqrt(etai &	
								*etaj))/(etai+etaj))**(real(l)+1.5)
					tmp_c = conjg(C1i)*C2j*((real(l)+1.5)/(conjg(etai) &
								+conjg(etaj)))*((2.d0*sqrt(conjg(etai) &	
								*conjg(etaj)))/(conjg(etai)+conjg(etaj)))**(real(l)+1.5)
					tmp_d = conjg(C2i)*C2j*((real(l)+1.5)/(conjg(etaj) &
					       		+etai))*((2.d0*sqrt(conjg(etaj) &	
				         		*etai))/(conjg(etaj)+etai))**(real(l)+1.5)
					self(i,j)=real(tmp_a+tmp_b+tmp_c+tmp_d,kind=range)
				end do
			end do
		case(3)
			print*, 'Hi!'
			do i=1, nfancy
				call Gbasis%get_eta(i, etai)
				call Gbasis%get_C1(i, C1i)
				call Gbasis%get_C2(i, C2i)
				do j=1, nfancy
					call Gbasis%get_eta(j, etaj)
					call Gbasis%get_C1(j, C1j)
					call Gbasis%get_C2(j, C2j)
					tmp_a= conjg(C1i)*C1j*((2.d0)/(sqrt(pi)))* &
								((2.d0**l*dble(fact(l)))/(dble(doublefact(2*l+1))))*&
								sqrt(conjg(etai)+etaj)*((2.0*sqrt(conjg(etai)*etaj))/(conjg(etai)+etaj))**(dble(l)+1.5d0)
					tmp_b= conjg(C2i)*C1j*((2.d0)/(sqrt(pi)))* &
								((2.d0**l*dble(fact(l)))/(dble(doublefact(2*l+1))))*&
								sqrt(etaj+etai)*((2.0*sqrt(etaj*etai))/(etaj+etai))**(real(l)+1.5)
					tmp_c= conjg(C1i)*C2j*((2.d0)/(sqrt(pi)))* &
								((2.d0**l*real(fact(l)))/(real(doublefact(2*l+1))))*&
								sqrt(conjg(etaj)+conjg(etai))*((2.0d0*sqrt(conjg(etaj)*conjg(etai)))/(conjg(etaj)+conjg(etai)))**(dble(l)+1.5d0)
					tmp_d= conjg(C2i)*C2j*((2.d0)/(sqrt(pi)))* &
								((2.d0**l*dble(fact(l)))/(dble(doublefact(2*l+1))))*&
								sqrt(conjg(etaj)+etai)*((2.0*sqrt(conjg(etaj)*etai))/(conjg(etaj)+etai))**(dble(l)+1.5d0)
					self(i,j)=real(tmp_a+tmp_b+tmp_c+tmp_d,kind=range)
				end do 
			end do
		case(4)
			do i=1, nfancy
				call Gbasis%get_eta(i, etai)
				call Gbasis%get_C1(i, C1i)
				call Gbasis%get_C2(i, C2i)
				do j=1, nfancy
					call Gbasis%get_eta(j, etaj)
					call Gbasis%get_C1(j, C1j)
					call Gbasis%get_C2(j, C2j)
					tmp_a = conjg(C1i)*C1j*((2.0*sqrt(conjg(etai)* &
								etaj))/(conjg(etai)+etaj+mu))**(real(l)+1.5)
					tmp_b = conjg(C2i)*C1j*((2.0*sqrt(etai* &
								etaj))/(etai+etaj+mu))**(real(l)+1.5)
					tmp_c = conjg(C1i)*C2j*((2.0*sqrt(conjg(etai)* &
								conjg(etaj)))/(conjg(etai)+conjg(etaj)+mu))**(real(l)+1.5)
					tmp_d = conjg(C2i)*C2j*((2.0*sqrt(conjg(etaj)* &
								etai))/(conjg(etaj)+etai+mu))**(real(l)+1.5)
					self(i,j)=real(tmp_a+tmp_b+tmp_c+tmp_d,kind=range)
				end do 
			end do
		case default 
			self=0
		end select		
		! Attractive potential
		self=-self
	end subroutine new_V_c

	subroutine new_S_c(self, nfancy, Gbasis)
		type(basis_c) :: Gbasis
		integer :: i, j, l
		integer, intent(in) :: nfancy
		real(kind=range), intent(inout) :: self(nfancy,nfancy)
		complex(kind=range) :: etai, etaj, tmp_a, tmp_b, tmp_c, tmp_d,&
				C1i, C2i, C1j, C2j
		real(kind=range) :: p
		call Gbasis%get_ell(l)
		!print*,'S'
		do i=1, nfancy
			call Gbasis%get_C1(i, C1i)
			call Gbasis%get_C2(i, C2i)
			call Gbasis%get_eta(i, etai)
			do j=1, nfancy
				call Gbasis%get_eta(j, etaj)
				call Gbasis%get_C1(j, C1j)
				call Gbasis%get_C2(j, C2j)
				tmp_a = conjg(C1i)*C1j*((2.0d0*sqrt(conjg(etai)*etaj))/(conjg(etai)+etaj))**(dble(l)+1.5d0)
				tmp_b = conjg(C2i)*C1j*((2.0d0*sqrt(etaj*etai))/(etaj+etai))**(dble(l)+1.5d0)
				tmp_c = conjg(C1i)*C2j*((2.0d0*sqrt(conjg(etai)*conjg(etaj)))/(conjg(etai)+conjg(etaj)))**(dble(l)+1.5d0)
				tmp_d = conjg(C2i)*C2j*((2.0d0*sqrt(conjg(etaj)*etai))/(conjg(etaj)+etai))**(dble(l)+1.5d0)
				self(i,j)=real(tmp_a+tmp_b+tmp_c+tmp_d,kind=range)
			end do 
		end do
	end subroutine new_S_c

!**********************************************************************
! Aux functions and deconstructor									  *
!**********************************************************************
! We construct some getters to access the elements
	subroutine print_HTVS(self)
		class(HTVS), intent(in) :: self
		print*, 'H'
			print*, self%H
		print*, 'T'
			print*, self%T
		print*, 'V'
			print*, self%V
		print*, 'S'
			print*, self%S
	end subroutine
	
	subroutine destruct(self)
		class(HTVS) :: self
		if (allocated(self%H)) deallocate(self%H)
		if (allocated(self%T)) deallocate(self%T)
		if (allocated(self%V)) deallocate(self%V)
		if (allocated(self%S)) deallocate(self%S)
	end subroutine destruct

	subroutine delete_HTVS(self)
		type(HTVS) :: self
		if (allocated(self%H)) deallocate(self%H)
		if (allocated(self%T)) deallocate(self%T)
		if (allocated(self%V)) deallocate(self%V)
		if (allocated(self%S)) deallocate(self%S)
	end subroutine delete_HTVS

end module HTVS_class
