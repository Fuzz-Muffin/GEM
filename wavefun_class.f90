module wavefun_class
	use para
	use grid_class
	implicit none 
	private
	
	type, public :: wavefunction
		real(kind=range), allocatable,dimension(:) :: psi
	contains
		procedure 	:: new => new_wavefun
		procedure	:: destruct
		final		:: delete_wavefun
	end type wavefun

contains
	
	! Construct wave function set
	subroutine new_wavefun(self, eigenvec, basis, grid, V)
		class(wavefun), intent(inout) :: self
		integer :: i, j, alstat, nfancy
	   	integer, intent(in) :: Rlength, V
		real(kind=range) :: tmp
		real(kind=range), dimension(:,:), intent(in) :: eigenvec, basis
		
		nfancy=size(eigenvec,1)
		print*,'nfancy in wavefun_class: ',nfancy

		if (allocated(self%psi)) deallocate(self%psi)
		allocate(self%psi(0:Rlength,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'ERROR: check Psi allocation'

		! Setting the first column of basisval to be the radial values 
		self%psi=0
		do i=1, Rlength
			do j=1,nfancy
				self%psi(i,j)=sum(eigenvec(:,j)*basis(i,:))
			end do 
		end do 
	
		! Here we check the physicallity of the wave function
		! For hydrogenic potential we want the first bit to be +	
		if (V==3) then
			do i=1, nfancy
				if (self%psi(1,i) <0) then
					if (self%psi(2,i) <0) then
						self%psi(:,i)=-self%psi(:,i)
					end if
				end if
			end do
		else
		! Otherwise we sum and check if it is + dominant  
			do i=1, nfancy
				tmp=sum(self%psi(:,i))
			   	if (tmp<=0) then
					self%psi(:,i)=-self%psi(:,i)
				end if
			end do 	
		end if
	end subroutine new_wavefun

	subroutine destruct(self)
		class(wavefun) :: self
		if (allocated(self%psi)) deallocate(self%psi)
	end subroutine destruct
	
	subroutine delete_wavefun(self)
		type(wavefun) :: self
		if (allocated(self%psi)) deallocate(self%psi)
	end subroutine delete_wavefun

end module wavefun_class
