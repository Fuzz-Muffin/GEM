module eigen_class
	use para 
	implicit none
	private

	type, public :: eigen
		real(kind=range), allocatable,dimension(:) :: eigenval
		real(kind=range), allocatable,dimension(:,:) :: eigenvec
	contains
		procedure 	:: new => new_eigen
		procedure 	:: get_eigenval, get_eigenvec, print_eigenvec, print_eigenval, destruct
		final		:: delete_eigen
	end type eigen

contains
	
	! Constructor
	subroutine new_eigen(self, H, S, verbose)
		class(eigen), intent(inout) :: self
		logical :: verbose
		integer :: lwork, info, alstat, i, j, nfancy
		real(kind=range), dimension(:,:), intent(in) :: H, S
		real(kind=range), allocatable,dimension(:,:) :: B
		real(kind=range), allocatable,dimension(:) :: work
		
		nfancy=size(H,1)

		if (allocated(self%eigenvec)) deallocate(self%eigenvec)
		allocate(self%eigenvec(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'ERROR: check eigenvec'
		self%eigenvec=0.0

		if (allocated(B)) deallocate(B)
		allocate(B(nfancy,nfancy), STAT=alstat)
			if (alstat/=0) STOP 'ERROR: check B'
		B=0
		 
		if (allocated(self%eigenval)) deallocate(self%eigenval)
		allocate(self%eigenval(nfancy), STAT=alstat)
			if (alstat/=0) STOP 'ERROR: check eigenval'

		! Build A from H
		do i=1, nfancy
			do j=i, nfancy
				self%eigenvec(i,j)=H(i,j)
			end do
		end do

		! Build B from S
		do i=1, nfancy
			do j=i, nfancy
				B(i,j)=S(i,j)
			end do
		end do 

		if (allocated(work)) deallocate(work)
		allocate(work(1))

		!print*,'nfancy', nfancy

		! We probe DSYGV to find out the dimension of the work array 
		call DSYGV(1, 'V', 'U', nfancy, self%eigenvec, nfancy, B, nfancy, self%eigenval, work, -1, info)
			
		if (info/=0) then
			print*, 'ERROR: problem with the first DSYGV call'
			print*, 'WORK(1)=', work(1), ', INFO=', info
			STOP
		end if	

		lwork=work(1)
		deallocate(work)
		allocate(work(lwork), STAT=alstat)
			if (alstat/=0) then
				print*, 'ERROR: check work allocation'
				STOP
			end if
		!print*, H
		!print*, self%eigenvec
		!print*, B
		
		! This is literally diagonalising the Hamiltonian, amazing. 
		call DSYGV(1, 'V', 'U', nfancy, self%eigenvec, nfancy, B, nfancy, self%eigenval, work, lwork, info)

		! Troubleshooting the DSYGV call, if it fails we halt the program
		if (verbose) then
			if (info<0) then
				print*, 'Problem with SGEEV perameter'
				select case (info)
					case(-1)
						print*, 'ITYPE is of illeagal value'
					case(-2)
						print*, 'JOBZ is of illeagal value'
					case(-3)
						print*, 'UPLO is of illeagal value'
					case(-4)
						print*, 'N is of illeagal value'
					case(-5)
						print*, 'A is of illeagal value'
					case(-6)
						print*, 'LDA is of illeagal value'
					case(-7)
						print*, 'B is of illeagal value'
					case(-8)
						print*, 'LDB is of illeagal value'
					case(-9)
						print*, 'W is of illeagal value'
					case(-10)
						print*, 'WORK is of illeagal value'
					case(-11)
						print*, 'LWORK is of illeagal value'
					case default
						print*, 'ERROR with INFO of SGEEV'
				end select
				STOP
			else if (INFO>0) then
				print*, 'INFO=', info
				print*, 'The factorization of B failed, and the eigenvalues, eigenvectors have not been computed'
				STOP
			else
				print*, 'DYSGV subroutine was successful'
			end if
		end if
	end subroutine new_eigen

	! Getters
	subroutine get_eigenval(self, i, res)
		class(eigen), intent(in) :: self
		integer :: i
		real(kind=range), intent(out) :: res
		res=self%eigenval(i)
	end subroutine get_eigenval

	subroutine get_eigenvec(self, i, v)
		class(eigen), intent(in) :: self
		integer :: i, n		
		real(kind=range), allocatable,dimension(:),intent(out) :: v	
		n=size(self%eigenvec,1)
		allocate(v(n))
		v=self%eigenvec(:,i)
	end subroutine get_eigenvec

	subroutine print_eigenval(self, i)
		class(eigen), intent(in) :: self
		integer, intent(in) :: i
		print'(a,I3,a,2E15.5)', 'eigen value #',i,'=',self%eigenval(i)
	end subroutine print_eigenval

	subroutine print_eigenvec(self, i)
		class(eigen), intent(in) :: self
		integer, intent(in) :: i
		integer :: j, n
		n=size(self%eigenvec(:,i))
		print*, 'Eigenvector',i
		do j=1,n
			print'(2E15.5)', self%eigenvec(j,i)
		end do
	end subroutine print_eigenvec
	
	subroutine destruct(self)
		class(eigen) :: self
		if (allocated(self%eigenvec)) deallocate(self%eigenvec)
		if (allocated(self%eigenval)) deallocate(self%eigenval)
	end subroutine destruct

	subroutine delete_eigen(self)
		type(eigen) :: self
		if (allocated(self%eigenvec)) deallocate(self%eigenvec)
		if (allocated(self%eigenval)) deallocate(self%eigenval)
	end subroutine delete_eigen

end module eigen_class
