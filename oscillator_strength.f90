module oscillator_strength
	use para
	use basis_class
	use basis_class_laguerre 
	use grid_class
	use eigen_class
	implicit none 
	public

interface os_calc_numerical
	module procedure os_calc_reg_n, os_calc_comp_n, os_calc_lag_n
end interface os_calc_numerical

interface os_calc_analytic
	module procedure os_calc_reg_a, os_calc_comp_a
end interface os_calc_analytic

contains

	!------------------------------------------------------------------------!
	! Calculation of the oscillator strength and static atomic polarizability!
	! Via Simpson's method numerical integration 							 !
	!------------------------------------------------------------------------!

	subroutine os_calc_reg_n(Gbasis, energies, psi, grid)
		type(basis), intent(in), dimension(:) :: Gbasis
		type(GridObject), intent(in) :: grid
		real(kind=range), intent(in), dimension(:,:,:) :: psi
		real(kind=range), intent(in), dimension(:,:) :: energies 
		real(kind=range), allocatable, dimension(:) :: f, a
		real(kind=range) :: tmp
		integer :: l_i, l_f, k, max_tmp, min_tmp, tmp_a, tmp_b, tmp_c, tmp_d
		integer :: i, counter=0
		open(2, file='oscillator_strength_numerical.txt', action='write', status='replace')
		write(2,*) '	1s -> np'
		call Gbasis(1)%get_ell(1, l_i)
		do i=1,size(energies(:,1))
			if (energies(i,2)<0.0) counter=counter+1
		end do
	   	allocate(f(size(Energies(:,1)))); allocate(a(size(Energies(:,1)))); a=0.0; f=0.0
		do i=1, counter
			tmp=0.0
			call Gbasis(2)%get_ell(1, l_f)
			! Simpson integration
			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
					grid%rgrid(min_tmp:max_tmp)**3*grid%weight(min_tmp:max_tmp))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
		end do
		write(2,*) 'Continuum'
		do i=counter+1, size(energies(:,1))
			tmp=0.0
			call Gbasis(2)%get_ell(1, l_f)
			! Simpson integration
			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
					grid%rgrid(min_tmp:max_tmp)**3*grid%weight(min_tmp:max_tmp))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(2E15.5)') f(i)
		end do
			
		write(2,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:counter))
		write(2,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(counter+1:size(f)))
		write(2,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,size(Energies(:,1))
			a(i)=f(i)*(energies(i,2)-energies(1,1))**(-2)
			write(2,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', sum(a(1:i))
		end do
		write(2,*) ''
		write(2,'(a,E15.5,a)') 'static dipole polarizability is thus a=', sum(a)*0.148184*10.0**24,' cm^3'
		close(2)
		deallocate(f); deallocate(a)
	end subroutine os_calc_reg_n

	subroutine os_calc_comp_n(Gbasis, energies, psi, grid)
		type(basis_c), intent(in), dimension(:) :: Gbasis
		type(GridObject), intent(in) :: grid
		real(kind=range), intent(in), dimension(:,:,:) :: psi
		real(kind=range), intent(in), dimension(:,:) :: energies 
		real(kind=range), allocatable, dimension(:) :: f, a
		real(kind=range) :: tmp
		integer :: l_i, l_f, k, tmp_a, tmp_b, tmp_c, tmp_d, max_tmp, min_tmp
		integer :: i, counter=0
		open(2, file='oscillator_strength_numerical.txt', action='write', status='replace')
		write(2,*) '	1s -> np'
		call Gbasis(1)%get_ell(1, l_i)
		do i=1,size(energies(:,1))
			if (energies(i,2)<0.0) counter=counter+1
		end do
	   	allocate(f(size(Energies(:,1)))); allocate(a(size(Energies(:,1)))); a=0.0; f=0.0
		do i=1, counter
			tmp=0.0
			call Gbasis(2)%get_ell(1, l_f)
			! Simpson integration
			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
					grid%rgrid(min_tmp:max_tmp)**3*grid%weight(min_tmp:max_tmp))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
		end do
		write(2,*) 'Continuum'
		do i=counter+1, size(energies(:,1))
			tmp=0.0
			call Gbasis(2)%get_ell(1, l_f)
			! Simpson integration
			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
					grid%rgrid(min_tmp:max_tmp)**3*grid%weight(min_tmp:max_tmp))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(2E15.5)') f(i)
		end do
			
		write(2,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:counter))
		write(2,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(counter+1:size(f)))
		write(2,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,size(Energies(:,1))
			a(i)=f(i)*(energies(i,2)-energies(1,1))**(-2)
			write(2,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', sum(a(1:i))
		end do
		write(2,*) ''
		write(2,'(a,E15.5,a)') 'static dipole polarizability is thus a=', sum(a)*0.148184*10.0**24,' cm^3'
		close(2)
		deallocate(f); deallocate(a)
	end subroutine os_calc_comp_n
	
	subroutine os_calc_lag_n(Lbasis, energies, psi, grid)
		type(BasisObject), intent(in), dimension(:) :: Lbasis
		type(GridObject), intent(in) :: grid
		real(kind=range), intent(in), dimension(:,:,:) :: psi
		real(kind=range), intent(in), dimension(:,:) :: energies 
		real(kind=range), allocatable, dimension(:) :: f, a
		real(kind=range) :: tmp
		integer :: l_i, l_f, k, max_tmp, min_tmp, tmp_a, tmp_b, tmp_c, tmp_d
		integer :: i, counter=0
		open(2, file='oscillator_strength_numerical.txt', action='write', status='replace')
		write(2,*) '	1s -> np'
		l_i=Lbasis(1)%b(1)%l
		do i=1,size(energies(:,1))
			if (energies(i,2)<0.0) counter=counter+1
		end do
	   	allocate(f(size(Energies(:,1)))); allocate(a(size(Energies(:,1)))); a=0.0; f=0.0
		do i=1, counter
			tmp=0.0
			l_f=Lbasis(2)%b(1)%l
			! Simpson integration
!			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
!			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
!			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
!			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
!					grid%rgrid(min_tmp:max_tmp)**2*grid%weight(min_tmp:max_tmp))
			tmp=sum(psi(1,:,1)*psi(2,:,i)*grid%rgrid(:)**2*grid%weight(:))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			print*, f(i)	
			write(2,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
		end do
		write(2,*) 'Continuum'
		do i=counter+1, size(energies(:,1))
			tmp=0.0
			!l_f=Lbasis(2)%b(1)%l
			! Simpson integration
!			call minmaxi(psi(1,:,i), grid, tmp_b, tmp_a)
!			call minmaxi(psi(2,:,i), grid, tmp_d, tmp_c)
!			max_tmp=min(tmp_a, tmp_c); min_tmp=max(tmp_b, tmp_d)
!			tmp=sum(psi(1,min_tmp:max_tmp,1)*psi(2,max_tmp:min_tmp,i)*&
!					grid%rgrid(min_tmp:max_tmp)**2*grid%weight(min_tmp:max_tmp))
			tmp=sum(psi(1,:,1)*psi(2,:,i)*grid%rgrid(:)**2*grid%weight(:))
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(2E15.5)') f(i)
		end do
			
		write(2,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:counter))
		write(2,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(counter+1:size(f)))
		write(2,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,size(Energies(:,1))
			a(i)=f(i)*(energies(i,2)-energies(1,1))**(-2)
			write(2,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', sum(a(1:i))
		end do
		write(2,*) ''
		write(2,'(a,E15.5,a)') 'static dipole polarizability is thus a=', sum(a)*0.148184*10.0**24,' cm^3'
		close(2)
		deallocate(f); deallocate(a)
	end subroutine os_calc_lag_n
	!------------------------------------------------------------------------!

	!------------------------------------------------------------------------!
	! Calculation of the oscillator strength and static atomic polarizability!
	! Via analytic method 				 									 !
	!------------------------------------------------------------------------!
	subroutine os_calc_reg_a(Gbasis, eigenstuff)
		type(basis), intent(in), dimension(:) :: Gbasis
	   	type(eigen), intent(in), dimension(:) :: eigenstuff
		integer :: l_i, l_f, N, boundstates, i, j, k, counter=0
		real(kind=range), allocatable, dimension(:) :: f, a
		real(kind=range) :: nu_j, nu_k, Nj, Nk, tmp_a, tmp_b,&
							Ei, Ef, tmp=0.0
		open(3, file='oscillator_strength_analytic.txt', action='write', status='replace')
		write(3,*) '	1s -> np'
		call Gbasis(1)%get_ell(1, l_i)
		call Gbasis(2)%get_ell(1, l_f)
		call Eigenstuff(1)%get_eigenval(1, Ei)
		N=size(eigenstuff(1)%eigenval(:))
		do i=1,N
			if (eigenstuff(1)%eigenval(i)<0.0) counter=counter+1
		end do
		boundstates=counter
		counter=0
	   	allocate(f(N)); allocate(a(N)); a=0.0; f=0.0
		do i=1, N
			call eigenstuff(2)%get_eigenval(i, Ef)
			tmp=0.0
			do j=1, N
					call Gbasis(1)%get_N(j,Nj)
					call Gbasis(1)%get_nu(j, nu_j)
				do k=1, N
					call Gbasis(2)%get_N(k,Nk)
					call Gbasis(2)%get_nu(k,nu_k)
					tmp_a=eigenstuff(1)%eigenvec(j,1)
					tmp_b=eigenstuff(2)%eigenvec(k,i)
					tmp=tmp+(tmp_a*tmp_b*Nj*Nk*&
							((sqrt(pi)*real(doublefact(l_i+l_f+2)))/&
							(2.d0**(0.5*(l_i+l_f+5))*(nu_j+nu_k)**(0.5*(l_i+l_f+4)))))
					end do
				end do
			f(i)=(2.0/3.0)*(Ef-Ei)*(tmp)**2
			if (i<=boundstates) then
				write(3,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
			else
				write(3,'(2E15.5)') f(i)
			end if 
		end do
		write(3,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:boundstates))
		write(3,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(boundstates+1:size(f)))
		write(3,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,N
			call eigenstuff(1)%get_eigenval(1,Ei)
			call eigenstuff(2)%get_eigenval(i,Ef)
			a(i)=f(i)*(Ef-Ei)**(-2)
			write(3,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', sum(a(1:i))
		end do
		write(3,*) ''
		write(3,'(a,E15.5,a)') 'static dipole polarizability is thus a=',	sum(a)*0.148184*10.0**(24),' cm^3'
		close(3)
		deallocate(f); deallocate(a)
	end subroutine os_calc_reg_a

	subroutine os_calc_comp_a(Gbasis, eigenstuff)
		type(basis_c), intent(in), dimension(:) :: Gbasis
	   	type(eigen), intent(in), dimension(:) :: eigenstuff
		integer :: l_i, l_f, N, boundstates, i, j, k, counter=0
		real(kind=range), allocatable, dimension(:) :: f, a
		complex(kind=range) :: C1j, C2j, C1k, C2k,Nj,Nk,w,x,y,z,eta_j,eta_k
		real(kind=range) :: tmp_a, tmp_b, Ei, Ef, tmp=0.0
		open(3, file='oscillator_strength_analytic.txt', action='write', status='replace')
		write(3,*) '	1s -> np'
		call Gbasis(1)%get_ell(1, l_i)
		call Gbasis(2)%get_ell(1, l_f)
		call Eigenstuff(1)%get_eigenval(1, Ei)
		N=size(eigenstuff(1)%eigenval(:))
		do i=1,N
			if (eigenstuff(1)%eigenval(i)<0.0) counter=counter+1
		end do
		boundstates=counter
		counter=0
	   	allocate(f(N)); allocate(a(N)); a=0.0; f=0.0
		do i=1, N
			call eigenstuff(2)%get_eigenval(i, Ef)
			tmp=0.0
			do j=1, N
					call Gbasis(1)%get_N(j, Nj)
					call Gbasis(1)%get_eta(j, eta_j)
					call Gbasis(1)%get_C1(j, C1j)
					call Gbasis(1)%get_C2(j, C2j)
				do k=1, N
					call Gbasis(2)%get_N(k,Nk)
					call Gbasis(2)%get_eta(k,eta_k)
					call Gbasis(2)%get_C1(k, C1k)
					call Gbasis(2)%get_C2(k, C2k)
					tmp_a=eigenstuff(1)%eigenvec(j,1)
					tmp_b=eigenstuff(2)%eigenvec(k,i)
					w=conjg(C1j)*C1k*(conjg(Nj)*Nk*&
							((sqrt(pi)*real(doublefact(l_i+l_f+2)))/&
							(2.d0**(0.5*(l_i+l_f+5))*(conjg(eta_j)+eta_k)**(0.5*(l_i+l_f+4)))))
					x=conjg(C2j)*C1k*(Nj*Nk*&
							((sqrt(pi)*real(doublefact(l_i+l_f+2)))/&
							(2.d0**(0.5*(l_i+l_f+5))*(eta_j+eta_k)**(0.5*(l_i+l_f+4)))))
					y=conjg(C1j)*C2k*(conjg(Nj)*conjg(Nk)*&
							((sqrt(pi)*real(doublefact(l_i+l_f+2)))/&
							(2.d0**(0.5*(l_i+l_f+5))*(conjg(eta_j)+conjg(eta_k))**(0.5*(l_i+l_f+4)))))
					z=conjg(C2j)*C2k*(Nj*conjg(Nk)*&
							((sqrt(pi)*real(doublefact(l_i+l_f+2)))/&
							(2.d0**(0.5*(l_i+l_f+5))*(eta_j+conjg(eta_k))**(0.5*(l_i+l_f+4)))))
					tmp=tmp+(tmp_a*tmp_b*real(w+x+y+z,kind=range))
				end do
			end do
			f(i)=(2.0/3.0)*(Ef-Ei)*(tmp)**2
			if (i<=boundstates) then
				write(3,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
			else
				write(3,'(2E15.5)') f(i)
			end if 
		end do
		write(3,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:boundstates))
		write(3,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(boundstates+1:size(f)))
		write(3,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,N
			call eigenstuff(1)%get_eigenval(1,Ei)
			call eigenstuff(2)%get_eigenval(i,Ef)
			a(i)=f(i)*(Ef-Ei)**(-2)
			write(3,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', sum(a(1:i))
		end do
		write(3,*) ''
		write(3,'(a,E15.5,a)') 'static dipole polarizability is thus a=',	sum(a)*0.148184*10.0**(24),' cm^3'
		close(3)
		deallocate(f); deallocate(a)
	end subroutine os_calc_comp_a
end module oscillator_strength				
