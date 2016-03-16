module oscillator_strength
	use para
	use basis_class
	use wavefun_class
	use eigen_class
	implicit none 
	public

interface os_calc_numerical
	module procedure os_calc_reg_n, os_calc_comp_n
end interface os_calc_numerical

interface os_calc_analytic
	module procedure os_calc_reg_a, os_calc_comp_a
end interface os_calc_analytic

contains

	!------------------------------------------------------------------------!
	! Calculation of the oscillator strength and static atomic polarizability!
	! Via Simpson's method numerical integration 							 !
	!------------------------------------------------------------------------!

	subroutine os_calc_reg_n(Gbasis, energies, wavefuns, rgrid)
		type(basis), intent(in), dimension(:) :: Gbasis
		type(wavefun), intent(in), dimension(:) :: wavefuns
		real(kind=range), intent(in),dimension(:,:) :: energies 
		real(kind=range), allocatable, dimension(:) :: rgrid, f, a
		real(kind=range) :: tmp, tmp_a, tmp_b, h, x0, x1, x2, x
		integer :: l_i, l_f, rmax, k
		integer :: i, counter=0
		open(2, file='oscillator_strength_numerical.txt', action='write', status='replace')
		write(2,*) '	1s -> np'
		call Gbasis(1)%get_ell(l_i)
		do i=1,size(energies(:,1))
			if (energies(i,2)<0.0) counter=counter+1
		end do
	   	allocate(f(size(Energies(:,1)))); allocate(a(size(Energies(:,1)))); a=0.0; f=0.0
		do i=1, counter
			tmp=0.0; x0=0.0
			call Gbasis(2)%get_ell(l_f)
			! Simpson integration
			h=rgrid(2)-rgrid(1)	
			rmax=size(rgrid)
			x1=0.0; x2=0.0
			x0=(0.0)+(wavefuns(1)%psi(rmax,1)*wavefuns(2)%psi(rmax,i)*rmax**3)
			do k=1, (rmax-1)
				x=real(k*h)
				tmp_a=wavefuns(1)%psi(k,1)
				tmp_b=wavefuns(2)%psi(k,i)
				if (mod(k,2)==0) then
					x2=x2+(tmp_a*tmp_b*x**3)
				else
					x1=x1+(tmp_a*tmp_b*x**3)
				end if
			end do
			tmp=h*(x0+2.0*x2+4.0*x1)/3.0
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
		end do
		write(2,*) 'Continuum'
		do i=counter+1, size(energies(:,1))
			tmp=0.0; x0=0.0
			call Gbasis(2)%get_ell(l_f)
			! Simpson integration
			h=rgrid(2)-rgrid(1)	
			rmax=size(rgrid)
			x1=0.0; x2=0.0
			x0=(0.0)+(wavefuns(1)%psi(rmax,1)*wavefuns(2)%psi(rmax,i)*rmax**3)
			do k=1, (rmax-1)
				x=real(k*h)
				tmp_a=wavefuns(1)%psi(k,1)
				tmp_b=wavefuns(2)%psi(k,i)
				if (mod(k,2)==0) then
					x2=x2+(tmp_a*tmp_b*x**3)
				else
					x1=x1+(tmp_a*tmp_b*x**3)
				end if
			end do
			tmp=h*(x0+2.0*x2+4.0*x1)/3.0
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(2E15.5)') f(i)
		end do
			
		write(2,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:counter))
		write(2,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(counter+1:size(f)))
		write(2,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,size(Energies(:,1))
			a(i)=a(i-1)+f(i)*(energies(i,2)-energies(1,1))**(-2)
			write(2,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', a(i)
		end do
		write(2,*) ''
		write(2,'(a,E15.5,a)') 'static dipole polarizability is thus a=', (a(size(Energies(:,1))))*0.148184*10.0**24,' cm^3'
		close(2)
		deallocate(f); deallocate(a)
	end subroutine os_calc_reg_n

	subroutine os_calc_comp_n(Gbasis, energies, wavefuns, rgrid)
		type(basis_c), intent(in), dimension(:) :: Gbasis
		type(wavefun), intent(in), dimension(:) :: wavefuns
		real(kind=range), intent(in),dimension(:,:) :: energies 
		real(kind=range), allocatable, dimension(:) :: rgrid, f, a
		real(kind=range) :: tmp, tmp_a, tmp_b, h, x0, x1, x2, x
		integer :: l_i, l_f, rmax, k
		integer :: i, counter=0
		open(2, file='oscillator_strength_numerical.txt', action='write', status='replace')
		write(2,*) '	1s -> np'
		call Gbasis(1)%get_ell(l_i)
		do i=1,size(energies(:,1))
			if (energies(i,2)<0.0) counter=counter+1
		end do
	   	allocate(f(size(Energies(:,1)))); allocate(a(size(Energies(:,1)))); a=0.0; f=0.0
		do i=1, counter
			tmp=0.0; x0=0.0
			call Gbasis(2)%get_ell(l_f)
			! Simpson integration
			h=rgrid(2)-rgrid(1)	
			rmax=size(rgrid)
			x1=0.0; x2=0.0
			x0=(0.0)+(wavefuns(1)%psi(rmax,1)*wavefuns(2)%psi(rmax,i)*rmax**3)
			do k=1, (rmax-1)
				x=real(k*h)
				tmp_a=wavefuns(1)%psi(k,1)
				tmp_b=wavefuns(2)%psi(k,i)
				if (mod(k,2)==0) then
					x2=x2+(tmp_a*tmp_b*x**3)
				else
					x1=x1+(tmp_a*tmp_b*x**3)
				end if
			end do
			tmp=h*(x0+2.0*x2+4.0*x1)/3.0
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(I2,a,a,2E15.5)') i+l_f,'p','',f(i)
		end do
		write(2,*) 'Continuum'
		do i=counter+1, size(energies(:,1))
			tmp=0.0; x0=0.0
			call Gbasis(2)%get_ell(l_f)
			! Simpson integration
			h=rgrid(2)-rgrid(1)	
			rmax=size(rgrid)
			x1=0.0; x2=0.0
			x0=(0.0)+(wavefuns(1)%psi(rmax,1)*wavefuns(2)%psi(rmax,i)*rmax**3)
			do k=1, (rmax-1)
				x=real(k*h)
				tmp_a=wavefuns(1)%psi(k,1)
				tmp_b=wavefuns(2)%psi(k,i)
				if (mod(k,2)==0) then
					x2=x2+(tmp_a*tmp_b*x**3)
				else
					x1=x1+(tmp_a*tmp_b*x**3)
				end if
			end do
			tmp=h*(x0+2.0*x2+4.0*x1)/3.0
			f(i) = (2.0/3.0)*(energies(i,2)-energies(1,1))*(tmp)**2		
			!print*, f(i)	
			write(2,'(2E15.5)') f(i)
		end do
			
		write(2,'(a,2E15.5)') 'Sum of discrete spectrum :',sum(f(1:counter))
		write(2,'(a,2E15.5)') 'Sum of continous spectrum :',sum(f(counter+1:size(f)))
		write(2,'(a,2E15.5)') 'Sum of oscillator strengths :',sum(f)

		! Dipole polarizability
		do i=1,size(Energies(:,1))
			a(i)=a(i-1)+f(i)*(energies(i,2)-energies(1,1))**(-2)
			write(2,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', a(i)
		end do
		write(2,*) ''
		write(2,'(a,E15.5,a)') 'static dipole polarizability is thus a=',a(size(Energies(:,1)))*0.148184*10.0**(24),' cm^3'
		close(2)
		deallocate(f); deallocate(a)
	end subroutine os_calc_comp_n
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
		call Gbasis(1)%get_ell(l_i)
		call Gbasis(2)%get_ell(l_f)
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
			a(i)=a(i-1)+f(i)*(Ef-Ei)**(-2)
			write(3,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', a(i)
		end do
		write(3,*) ''
		write(3,'(a,E15.5,a)') 'static dipole polarizability is thus a=',	a(N)*0.148184*10.0**(24),' cm^3'
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
		call Gbasis(1)%get_ell(l_i)
		call Gbasis(2)%get_ell(l_f)
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
			a(i)=a(i-1)+f(i)*(Ef-Ei)**(-2)
			write(3,'(a,I3.3,a,2E15.5)') 'Static dipole polarizability using',i,'terms, a=', a(i)
		end do
		write(3,*) ''
		write(3,'(a,E15.5,a)') 'static dipole polarizability is thus a=',	a(N)*0.148184*10.0**(24),' cm^3'
		close(3)
		deallocate(f); deallocate(a)
	end subroutine os_calc_comp_a
end module oscillator_strength				
