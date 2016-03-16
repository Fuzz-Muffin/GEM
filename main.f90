!************************************************************************!
! This is the main code for GEM, can use either the Gaussian, complex- 	 !
! range, or Laguerre basis. 											 !
!------------------------------------------------------------------------!
! Filip Vukovic 														 !
! Curtin University, Bently 											 !
!************************************************************************!
program main
	use para				! These are constants and useful generic functions
	use toms707 			! The hypergeometric function, TOMS algorithm #707
	use input_class 		! This reads the data file and makes input object
	use grid_class 			! This is the grid object, written by D.Fursa
	use basis_class			! Constructs the basis required
	use basis_class_laguerre! Laguerre basis
	use HTVS_class			! This makes the matrices and elements
	use eigen_class			! This diags the Hamil and gives eigenval/vectors
	use gnufor2 			! This is to plot on the fly using gnuplot, smart.
	use oscillator_strength	! Subroutines for calculating oscilator strength
   	
	implicit none
	type(input) 		:: indata
	type(GridObject) 	:: grid 
	type(HTVS) 			:: HTVSmats
	type(eigen), 		allocatable, dimension(:) :: eigenstuff
	type(basis), 		allocatable, dimension(:) :: Gbasis
	type(BasisObject), 	allocatable, dimension(:) :: Lbasis
	type(basis_c), 		allocatable, dimension(:) :: Gbasis_c
	character(len=50) 	:: outfile
	logical 			:: ex
	integer 			:: i, j, k, h, p, counter, Rlength, alstat, N, l, ell, min_tmp, max_tmp
	real(kind=range) 	:: tmp,tmp_a,tmp_b,e_exp,lambda,a,t1,t2,hg
	real(kind=range), allocatable, dimension(:) 	:: Rgrid, tmp_fun
	real(kind=range), allocatable, dimension(:,:) 	:: basisval, expectation, overlap, energies
	real(kind=range), allocatable, dimension(:,:,:) :: psi
	complex(kind=range) :: tmp_c
	complex(kind=range), allocatable, dimension(:,:) :: basisval_c

	!---------------------------------------------------------------------!
	! Building the input object and setting things up				 	
	!---------------------------------------------------------------------!		
	! OPT's: input object, set to 10, verbose file writing =1/ silent=0 
	call indata%new(10, 1)

	! Some setting for the complex basis functions, if used
	a=indata%alpha

	! We use twice nfancy for complex basis, as we have pairs of functions
	if (indata%BasisType==0 .or. indata%BasisType==2) then
		N=indata%nfancy
	else if (indata%BasisType==1) then 
		N=2*indata%nfancy
	end if

	!---------------------------------------------------------------------!
	! Allocating the memory required and building the grid object		
	!---------------------------------------------------------------------!	
	call grid%setgrid(indata)

	select case (indata%BasisType) 
	case(0)
		allocate(Gbasis((indata%l_f-indata%l_i)+1),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check Gbasis_c'
	case(1)
		allocate(Gbasis_c((indata%l_f-indata%l_i)+1),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check Gbasis'
	case(2)
		allocate(Lbasis((indata%l_f-indata%l_i)+1),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check Lbasis'
	end select
	
	allocate(energies(N,indata%l_f-indata%l_i+1),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check energies allocation'

	allocate(eigenstuff(indata%l_f-indata%l_i+1),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check eigenstuff allocation'

	allocate(psi(indata%l_f-indata%l_i+1,grid%nr,N),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check psi allocation'

	!---------------------------------------------------------------------!
	! Here we loop over the angular momenta from l_i -> l_f 
	!---------------------------------------------------------------------!
	do l=1, (indata%l_f-indata%l_i+1)
		ell= indata%l_i+(l-1)
		print*,l

		!---------------------------------------------------------------------!
		! Let's build the desired basis 	
		!---------------------------------------------------------------------!
		select case (indata%BasisType)
		case(0) 
			call Gbasis(l)%new(ell, N, indata%r1, indata%rN, grid)
			if (indata%verbose) then
				do i=1, N
					print*,'Using regulair Gaussians'
					call Gbasis(l)%print_nu(i)
					call Gbasis(l)%print_N(i)
					call Gbasis(l)%print_f(i, grid)
				end do 
			end if 
		case(1)
			call Gbasis_c(l)%new(ell, N, indata%r1, indata%rN, a, grid)
			if (indata%verbose) then
			print*,'Using complex range Gaussians'
				do i=1, N
					call Gbasis_c(l)%print_a(i)
					call Gbasis_c(l)%print_eta(i) 
					call Gbasis_c(l)%print_N_c(i)
					call Gbasis_c(l)%print_f(i, grid)
				end do
			end if
		case(2)
			if (indata%verbose) print*, 'Using the Laguerre basis'
			call Lbasis(l)%constructLaguerreBasis(N, ell, indata%alpha_lag, grid)	
		end select
		!---------------------------------------------------------------------!
		
		!---------------------------------------------------------------------!
		! Here we plot some basis functions to see what is going on	 !
		!---------------------------------------------------------------------!
		if (indata%basisplots) then
			do i=1,N
			write(outfile, '(a,I3.3,a)') 'G',i,'.png'
			select case (indata%BasisType)
			case(0)
				call plot(real(grid%rgrid(:),kind=8),real(Gbasis(l)%G(i)%f(:),kind=8),filename=outfile,terminal='png')
			case(1) 		
				call plot(real(grid%rgrid(:),kind=8),real(Gbasis_c(l)%G(i)%f(:),kind=8),filename=outfile,terminal='png')
			case(2)
				call plot(real(grid%rgrid(:),kind=8),real(Lbasis(l)%b(i)%f(:),kind=8),filename=outfile,terminal='png')
			end select
			end do
		end if
		!---------------------------------------------------------------------!
	
		!---------------------------------------------------------------------!
		! Building H,T,V,S and diagonalising H	
		!---------------------------------------------------------------------!
		select case (indata%BasisType)
		case(0)
			call HTVSmats%new(indata%mu, N, Gbasis(l), indata%Vstat) 
		case(1)
			call HTVSmats%new(indata%mu, N, Gbasis_c(l), indata%Vstat) 
		case(2)
			call Lbasis(l)%constructOverlap(ell, indata%alpha_lag, grid)
			call Lbasis(l)%constructPotential(indata%alpha_lag, indata%vstat, grid)
			call Lbasis(l)%constructHamiltonian(ell, indata%alpha_lag, grid)
		end select
	
		if (indata%verbose .and. indata%BasisType<=1) call HTVSmats%print_HTVS()
		! Time to diagonalise the Hamiltonian...
		if (indata%BasisType==2) then
			call eigenstuff(l)%new(Lbasis(l)%H, Lbasis(l)%overlap, indata%verbose)
		else
			call eigenstuff(l)%new(HTVSmats%H, HTVSmats%S, indata%verbose)
		end if
		
		! Print eigenvectors to terminal in verbose mode for debug
		if (indata%verbose) then
			do i=1, N
				call eigenstuff(l)%print_eigenval(i)
				call eigenstuff(l)%print_eigenvec(i)	
			end do		
		end if
		!---------------------------------------------------------------------!
			
		!---------------------------------------------------------------------!
		! Building the wavefunctions of calculated eigen energies		
		!---------------------------------------------------------------------!
		do i=1, grid%nr
			do j=1, N
				select case (indata%BasisType)
				case(0)
					tmp=0.0
					do k=1,N
						tmp=tmp+eigenstuff(l)%eigenvec(k,j)*Gbasis(l)%G(k)%f(i)
					end do
					psi(l,i,j)=tmp
				case(1)
					tmp=0.0
					do k=1,N
						tmp=tmp+eigenstuff(l)%eigenvec(k,j)*Gbasis_c(l)%G(k)%f(i)
					end do
					psi(l,i,j)=tmp
				case(2)
					tmp=0.0
					do k=1,N
						tmp=tmp+eigenstuff(l)%eigenvec(k,j)*Lbasis(l)%b(k)%f(i)
					end do
					psi(l,i,j)=tmp
				end select
				tmp=0.0
			end do 
		end do
		!---------------------------------------------------------------------!
		! We check the physicallity of the wave functions 	
		!---------------------------------------------------------------------!
		do i=1, N
			if (psi(l,1,i)<0.0 .and. psi(l,2,i)<0.0) then
				psi(l,:,i) = -psi(l,:,i)
			end if
		end do
		!---------------------------------------------------------------------!

		!---------------------------------------------------------------------!
		! Write the wavefunctions and other data of the system to a file
		!---------------------------------------------------------------------!
		write(outfile, '(a,I1,a,I2.2,a,I3.3,a,F5.3,a,I4.4,a,I4.4,a)')&
			   'result_',indata%Vstat,'V_',ell,'l_',N,'n_', indata%r1,&
			   'r1_',int(indata%rN),'rmax_',int(indata%Rmax),'rad.txt'
		open(2, file=outfile, action='write', status='replace')
			write(2,*) '# Parameters of run: Potential=',indata%Vstat,&
					', Angular momentum=',ell,', Basis size=', N, ', r1=',&
				   	indata%r1,', r_max=', indata%rN
			write(2,*) '# r		wavefunctions 1->N			energy levels(eigenvalues)	analytic_energies'
			do i=0, grid%nr
				if (i<=N) then
					tmp=i+ell
					e_exp=-1d0/(2d0*tmp**2)
					call eigenstuff(l)%get_eigenval(i, tmp_a)
					if (tmp_a<0.0) then
						write(2,*) grid%rgrid(i), psi(l,i,:), tmp_a, e_exp
					else
						write(2,*) grid%rgrid(i), psi(l,i,:), tmp_a
					end if
				else
					write(2,*) grid%rgrid(i), psi(l,i,:)
				end if
			end do 
		close(2)	
		!---------------------------------------------------------------------!
	  
		!---------------------------------------------------------------------!
		! Save the Energy levels to a file and to energies array 
		!---------------------------------------------------------------------!
		energies(:,l)=eigenstuff(l)%eigenval(:)
		write(outfile, '(a,I2.2,a)') 'Energy_levels_l',ell,'.txt' 
		open(3, file=outfile, action='write', status='replace')
			if (indata%BasisType==3) then
				write(3,*) 'Laguerre Basis	Analytic'
				print*, 'Laguerre Basis	Analytic'
			else
				write(3,*) 'GEM method	Analytic'
				print*, 'GEM method	Analytic'
			end if 
			do i=1, N
				tmp=i+ell
				e_exp=-1d0/(2d0*tmp**2)
				call eigenstuff(l)%get_eigenval(i, lambda) 
				if (indata%Vstat==3 .and. lambda<0.0) then
					write(3,'(1P,2E15.5)') lambda, e_exp
					print'(I3, ES20.9,ES20.9)', ell+i, lambda, e_exp
				else
					write(3,'(2E15.5)') lambda
					print'(I3.3,ES20.9)', i, lambda
				end if
			end do
		close(3)
		!---------------------------------------------------------------------!
	
		!---------------------------------------------------------------------!
		! Calculate the overlap for the basis functions for free H	
		! <psi_nl|sin(kr)>											
		! Print the overlap values to a file 						
		!---------------------------------------------------------------------!
		if (indata%Vstat==1 .and. ell==0) then
			allocate (overlap(N,N))
			allocate (tmp_fun(grid%nr))
			do i=1, N
				call eigenstuff(l)%get_eigenval(i, tmp)
				do k=1, grid%nr
					tmp_fun(k)=sin(sqrt(2.0*tmp)*grid%rgrid(k))
				end do
				!print*, tmp_fun(:)
				do j=1, N
					call eigenstuff(l)%get_eigenval(i, tmp)
					do k=1, grid%nr
						tmp_fun(k)=sin(sqrt(2.0*tmp)*grid%rgrid(k))
					end do
					call minmaxi(psi(l,:,j), grid, min_tmp, max_tmp)
					tmp=0.0
					!print*, grid%rgrid(:)
					do k=min_tmp, max_tmp
						tmp=tmp+(psi(l,k,j)*tmp_fun(k)*grid%rgrid(k)*grid%weight(k))
					end do
					!print*, tmp
					overlap(j,i) = tmp
				end do
				if (indata%verbose) print*, 'E',i,' overlaps',overlap(i,:)
			end do 	
			deallocate(tmp_fun)
		!---------------------------------------------------------------------!
			write(outfile, '(a,I2.2,a)') 'Free_H_overlaps_l',ell,'.txt' 
			open(4, file=outfile, action='write', status='replace')
				write(4,*) '*------------------------------------------------&
						------------------------------------------------*'
				write(4,*) 'Parameters of run:'
				write(4,'(a,I1,a,I2.2,a,I1,a,I3.3,a,F4.3,a,F5.1,a,F5.1)') ' Potential=',&
						indata%Vstat,', Angular momentum=',ell,', Basis:',indata%BasisType,&
						', Basis size=', N, ', r1=', indata%r1,', r_max=',indata%rN,&
						', max radius:', indata%Rmax 
				if (indata%BasisType==1) then
					write(4,'(a,F4.3)'),'alhpa=',a
				end if
				write(4,*) '**************************************************&
						**************************************************'
				do i=1, N
					write(4,*) '*------------------------------------------------*'
					call eigenstuff(l)%get_eigenval(i, lambda)
					write(4,'(a, I3, a, 2E15.5)') ' Overlaps for state n=',i, ', E=', lambda
					write(4,*) '*------------------------------------------------*'
					write(4,*) ' n	<psi_n|sin(kr)>'
					write(4,*) ' ------------------'
					do j=1, N
						write(4, '(I3,2E15.5)') j, overlap(j,i)
					end do
					write(4,*) '**************************************************'
					write(4,*) ' '
				end do
			close(4)
			write(outfile, '(a,I2.2,a)') 'Energy_levels_overlaps_l',ell,'.txt' 
			open(5, file=outfile, action='write', status='replace')
				write(5,*) 'GEM energy	<psi_n|sin(kr)>'
				if (indata%overlaps) print*, 'GEM method 	<psi_n|sin(kr)>'
				do i=1,N
					call eigenstuff(l)%get_eigenval(i, lambda)
					write(5,'(1P,2E15.5)') lambda, overlap(i,i)
					if (indata%overlaps) print'(1p,2E15.5)', lambda, overlap(i,i)
				end do 
			close(5)
		end if
		!---------------------------------------------------------------------!
		
		
		!---------------------------------------------------------------------!
		! Calculate expectation values <r^p> = <psi_nl|r^p|psi_nl>		
		!---------------------------------------------------------------------!
		! p_max is the largest absolute power of r needed, set as required
	  	! We only want this calculation for the BOUND states for hydrogenic V
		!---------------------------------------------------------------------!
		if (indata%expect_r .and. indata%Vstat==3) then
			indata%p_max=2
			allocate(overlap(N,(indata%p_max*2)))
			counter=0
			do p=-indata%p_max,indata%p_max
				if (p/=0) then 
					counter=counter+1
					do i=1, N
						call eigenstuff(l)%get_eigenval(i, lambda)
						if (lambda<0) then
						call minmaxi(psi(l,:,i), grid, min_tmp, max_tmp)
						overlap(i,counter) = sum(psi(l,min_tmp:max_tmp,j)**2*&
								grid%rgrid(min_tmp:max_tmp)**(p+2))
						print'(a,I2.2,a,I2.2,a,I2,a,2E15.5)', 'n=',i+ell,&
									' l=',ell,' p=',p,' <r^p>=',overlap(i,counter)
						end if
					end do
				end if
			end do	
		end if	
		!---------------------------------------------------------------------!
		
		!---------------------------------------------------------------------!
		! Plotting the wavefunctions of calculated eigen energies	 
		!---------------------------------------------------------------------!
		if (indata%wavefunplots) then
			allocate(tmp_fun(grid%nr))
			select case (indata%BasisType)
			case(2)
	 		if (indata%Vstat==3) then
				do i=1, N
					!do j=0, Rlength
					!	tmp_fun(j)=R(i,ell,1,rgrid(j))
					!end do
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_r_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(rgrid(:)*tmp_fun(:),kind=8),&
							real(grid%rgrid(:),kind=8), real(grid%rgrid(:)*psi(l,:,i)/grid%rgrid(:),kind=8),&
							filename=outfile, terminal='png')
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(tmp_fun(:),kind=8),&
							real(grid%rgrid(:),kind=8), real(psi(l,:,i)/grid%rgrid(:),kind=8),&
							filename=outfile, terminal='png')
				end do
			else if (indata%Vstat==1 .and. ell==0) then
				do i=1, N
					do j=1, grid%nr 
						call eigenstuff(l)%get_eigenval(i,tmp)
						tmp_a=sqrt(2.0*tmp)
						tmp_fun(j)=sin(tmp_a*grid%rgrid(j))
					end do	
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'',ell,'_freeH.png'
					call plot(real(grid%rgrid(:),kind=range), real(psi(l,:,i)*real(grid%rgrid(:)),kind=range),&
								!real(grid%rgrid(:),kind=range), real(overlap(i,i)*psi(l,:,i),kind=range),&
						   		real(grid%rgrid(:),kind=range), tmp_fun(:),filename=outfile, terminal='png')
				end do
			end if

			case default
	 		if (indata%Vstat==3) then
				do i=1, N
					!do j=0, Rlength
					!	tmp_fun(j)=R(i,ell,1,rgrid(j))
					!end do
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_r_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(rgrid(:)*tmp_fun(:),kind=8),&
							real(grid%rgrid(:),kind=8), real(grid%rgrid(:)*psi(l,:,i),kind=8),&
							filename=outfile, terminal='png')
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(tmp_fun(:),kind=8),&
							real(grid%rgrid(:),kind=8), real(psi(l,:,i),kind=8),&
							filename=outfile, terminal='png')
				end do
			else if (indata%Vstat==1 .and. ell==0) then
				do i=1, N
					do j=1, grid%nr 
						call eigenstuff(l)%get_eigenval(i,tmp)
						tmp_a=sqrt(2.0*tmp)
						tmp_fun(j)=sin(tmp_a*grid%rgrid(j))
					end do	
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'',ell,'_freeH.png'
					call plot(real(grid%rgrid(:),kind=range), real(grid%rgrid(:)*psi(l,:,i),kind=range),&
								real(grid%rgrid(:),kind=range), real(overlap(i,i)*grid%rgrid(:)*psi(l,:,i),kind=range),&
						   		real(grid%rgrid(:),kind=range), tmp_fun(:),filename=outfile, terminal='png')
				end do
			end if
			end select
		deallocate(tmp_fun)
		end if 
		!---------------------------------------------------------------------!

		!---------------------------------------------------------------------!
		! Destructing the HTVS object so we can reuse them 		
		! all the important stuff has been extracted by now		
		!---------------------------------------------------------------------!
		call HTVSmats%destruct()

	end do
	!---------------------------------------------------------------------!
	! End of the angular momenta loop
	!---------------------------------------------------------------------!
	
	!---------------------------------------------------------------------!
	! Calculation of oscillator strength & static atomic polarizability
	!---------------------------------------------------------------------!
	if (indata%OScalc) then
		if (l<2) STOP 'ERROR: Need at least two bases to calc os'
		select case (indata%BasisType)
		case(0)
			if (indata%OSmeth==1 .or. indata%OSmeth==3) then
				call cpu_time(t1)
					call os_calc_numerical(Gbasis, energies, psi, grid)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for numerical oscillator strength calc', t2-t1, 'seconds'
			end if
			if (indata%OSmeth==0 .or. indata%OSmeth==3) then
				call cpu_time(t1)
					call os_calc_analytic(Gbasis, eigenstuff)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for analytic oscillator strength calc', t2-t1, 'seconds'
			end if
		case(1)
			if (indata%OSmeth==1 .or. indata%OSmeth==3) then
				call cpu_time(t1)
					call os_calc_numerical(Gbasis_c, energies, psi, grid)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for numerical oscillator strength calc', t2-t1, 'seconds'
			end if
			if (indata%OSmeth==0 .or. indata%OSmeth==3) then
				call cpu_time(t1)
					call os_calc_analytic(Gbasis_c, eigenstuff)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for analytic oscillator strength calc', t2-t1, 'seconds'
			end if
		case(2)
			if (indata%OSmeth==1 .or. indata%OSmeth==3) then
				call cpu_time(t1)
					call os_calc_numerical(Lbasis, energies, psi, grid)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for numerical oscillator strength calc', t2-t1, 'seconds'
			else if (indata%OSmeth==0) then
				STOP 'ERROR:No analytic method for Laguerre basis oscillator strength implimented'
			end if
		end select
	end if

end program main
