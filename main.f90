! Main code for the GEM method 
program main
	use para				! These are constants and useful generic functions
	use toms707 			! The hypergeometric function, TOMS algorithm #707
	use input_class 		! This reads the data file and makes input object
	use basis_class			! Constructs the basis required
	use HTVS_class			! This makes the matrices and elements
	use eigen_class			! This diags the Hamil and gives eigenval/vectors
	use wavefun_class		! This makes the wave-function set 
	use gnufor2 			! This is to plot on the fly using gnuplot
	use oscillator_strength	! Subroutines for calculating oscilator strength
   	
	implicit none
	type(input) 	:: indata
	type(HTVS) 		:: HTVSmats
	type(eigen), 	allocatable, dimension(:) 	:: eigenstuff
	type(wavefun), 	allocatable, dimension(:)	:: wavefuns
	type(basis), 	allocatable, dimension(:) 	:: Gbasis
	type(basis_c), 	allocatable, dimension(:) 	:: Gbasis_c
	character(len=50) :: outfile
	logical :: ex
	integer :: i, j, k, p, counter, Rlength, alstat, N, l, ell
	real(kind=range) :: stepsize, tmp,tmp_a,tmp_b,e_exp,lambda,a,h,x,x0,x1,x2,t1,t2,t3,hg
	real(kind=range), allocatable, dimension(:) :: Rgrid, tmp_fun
	real(kind=range), allocatable, dimension(:,:) :: basisval, expectation, overlap, energies
	complex(kind=range) :: tmp_c
	complex(kind=range), allocatable, dimension(:,:) :: basisval_c

	!------------------------------------------------------------------------!
	! Building the input object and setting things up				 		 !
	!------------------------------------------------------------------------!		
	! OPT's: input object, set to 10, verbose file writing =1/ silent=0 
	call indata%new(10, 1)

	! Some setting for the complex basis functions, if used
	a=indata%alpha

	! We use twice nfancy for complex basis, as we have pairs of fun's
	! Trap code here if BasisType is not 0 or 1
	if (indata%BasisType==0) then
		N=indata%nfancy
	else if (indata%BasisType==1) then 
		N=2*indata%nfancy
	else
		STOP 'ERROR: BasisType must be 0 or 1'
	end if

	!------------------------------------------------------------------------!
	! Allocating the memory required and setting the grids 					 !
	!------------------------------------------------------------------------!		
	! Build the grid to record the radial values of the basis functions
	Rlength=int(indata%Rmax/indata%stepsize)
		!if(indata%verbose) 
		print*, 'Rlength= ', Rlength

	! We use rgrid to store the radial points
	allocate(rgrid(0:Rlength), STAT=alstat)
		if (alstat/=0) STOP 'ERROR: check rgrid allocation'
	allocate(tmp_fun(0:Rlength), STAT=alstat)
		if (alstat/=0) STOP 'ERROR: check rgrid allocation'

	! Basisval stores the values of the G functions at the grid points
	if (indata%BasisType==1) then
		allocate(basisval_c(0:Rlength,N),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check basisval_c'
		allocate(Gbasis_c((indata%l_f-indata%l_i)+1),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check Gbasis'
	else if (indata%BasisType==0) then
		allocate(basisval(0:Rlength,N),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check basisval'
		allocate(Gbasis((indata%l_f-indata%l_i)+1),STAT= alstat)
			if (alstat/=0) STOP 'ERROR: check Gbasis_c'
	end if
	
	allocate(energies(N,indata%l_f-indata%l_i+1),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check energies allocation'

	allocate(wavefuns((indata%l_f-indata%l_i)+1),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check wavefuns allocation'

	allocate(eigenstuff(indata%l_f-indata%l_i+1),STAT= alstat)
		if (alstat/=0) STOP 'ERROR: check eigenstuff allocation'

	! Setting the radial values 
	do i=0, Rlength
		rgrid(i)=real(i)*indata%stepsize
	end do

	!------------------------------------------------------------------------!
	! Here we loop over the angular momenta from l_i -> l_f 				 !
	!------------------------------------------------------------------------!
	do l=1, (indata%l_f-indata%l_i+1)
		ell= indata%l_i+(l-1)

		print*, 'ell=',ell
		print*, 'l=', l

		!------------------------------------------------------------------------!
		! Let's build the basis onto the radial grid							 !
		!------------------------------------------------------------------------!
		if (indata%BasisType==1) then
			if (indata%verbose) print*,'Using complex range Gaussians'
			call Gbasis_c(l)%new(ell, N, indata%r1, indata%rN, a)
			do i=0, Rlength
				do j=1, N
					tmp_c=(0.0,0.0)
					call Gbasis_c(l)%calc(j, rgrid(i), tmp_c)
					basisval_c(i,j)=tmp_c
					if (indata%verbose) print*, 'r=',rgrid(i),'	G value=',basisval_c(i,j)
				end do 
			end do
			if (indata%verbose) then
				do j=1, N
					call Gbasis_c(l)%print_nu(j)
					call Gbasis_c(l)%print_eta(j) 
					call Gbasis_c(l)%print_N_c(j)
					call Gbasis_c(l)%print_N(j)
				end do
			end if
		else if (indata%BasisType==0) then
			if (indata%verbose) print*,'Using regular Gaussians'
			call Gbasis(l)%new(ell, N, indata%r1, indata%rN)
			do i=0, Rlength
				do j=1, N
					if (indata%verbose) then
						call Gbasis(l)%print_nu(j)
						call Gbasis(l)%print_N(j)
					end if 
					call Gbasis(l)%calc(j, rgrid(i), tmp)
					basisval(i,j)=tmp; tmp=0
					if (indata%verbose) print*, 'r=',rgrid(i),'	G value=',basisval(i,j)
				end do 
			end do 
		end if
		!------------------------------------------------------------------------!
		
		!------------------------------------------------------------------------!
		! Here we plot some Gaussian basis functions to see what is going on	 !
		!------------------------------------------------------------------------!
		if (indata%basisplots) then
			do i=1,N
			write(outfile, '(a,I3.3,a)') 'G',i,'.png'
				if (indata%BasisType==0) then
					call plot(real(rgrid(:),kind=8),real(basisval(:,i),kind=8),filename=outfile,terminal='png')
				else if (indata%BasisType==1) then		
					call plot(real(rgrid(:),kind=8),real(basisval_c(:,i),kind=8),filename=outfile,terminal='png')
				end if
			end do
		end if
		!------------------------------------------------------------------------!
	
		!------------------------------------------------------------------------!
		! Building H,T,V,S and diagonalising H									 !
		!------------------------------------------------------------------------!
		if (indata%BasisType==0) then
			call HTVSmats%new(indata%mu, N, Gbasis(l), indata%Vstat) 
		else if (indata%BasisType==1) then
			call HTVSmats%new(indata%mu, N, Gbasis_c(l), indata%Vstat) 
		end if
	
		if (indata%verbose) call HTVSmats%print_HTVS()
	
		! Time to diagonalise the Hamiltonian...
		call eigenstuff(l)%new(HTVSmats%H, HTVSmats%S, indata%verbose)
		
		! Print eigenvectors to terminal in verbose mode for debug
		if (indata%verbose) then
			do i=1, N
				call eigenstuff(l)%print_eigenval(i)
				call eigenstuff(l)%print_eigenvec(i)	
			end do		
		end if
		!------------------------------------------------------------------------!
	
		!------------------------------------------------------------------------!
		! Building the wavefunctions of calculated eigen energies				 !
		!------------------------------------------------------------------------!		
		if (indata%BasisType==0) then
			call wavefuns(l)%new(eigenstuff(l)%eigenvec, basisval, Rlength, indata%Vstat)
		else if (indata%BasisType==1) then
			call wavefuns(l)%new(eigenstuff(l)%eigenvec, real(basisval_c), Rlength, indata%Vstat)
		end if
		!------------------------------------------------------------------------!
		
		!------------------------------------------------------------------------!
		! Write the wavefunctions and other data of the system to a file		 !
		!------------------------------------------------------------------------!		
		write(outfile, '(a,I1,a,I2.2,a,I3.3,a,F5.3,a,I4.4,a,I4.4,a)') 'result_',indata%Vstat,'V_',ell,&
				'l_',N,'n_', indata%r1,'r1_',int(indata%rN),'rmax_',int(indata%Rmax),'rad.txt'
		open(2, file=outfile, action='write', status='replace')
			write(2,*) '# Parameters of run: Potential=',indata%Vstat,', Angular momentum=',ell,&
				', Basis size=', N, ', r1=', indata%r1,', r_max=', indata%rN
			write(2,*) '# r		wavefunctions 1->N			energy levels(eigenvalues)	analytic_energies'
	
			do i=0, Rlength
				if (i<=N) then
					tmp=i+ell
					e_exp=-1d0/(2d0*tmp**2)
					call eigenstuff(l)%get_eigenval(i, tmp_a)
					write(2,*) rgrid(i), wavefuns(l)%psi(i,:), tmp_a, e_exp
				else
					write(2,*) rgrid(i), wavefuns(l)%psi(i,:)
				end if
			end do 
		close(2)	
		!------------------------------------------------------------------------!
	
		!------------------------------------------------------------------------!
		! Save the Energy levels to a file and to energies array for safekeeping !
		!------------------------------------------------------------------------!		
		energies(:,l)=eigenstuff(l)%eigenval(:)
		write(outfile, '(a,I2.2,a)') 'Energy_levels_l',ell,'.txt' 
		open(3, file=outfile, action='write', status='replace')
			write(3,*) 'GEM method	Analytic'
			print*, 'GEM method	Analytic'
			do i=1, N
				tmp=i+ell
				e_exp=-1d0/(2d0*tmp**2)
				call eigenstuff(l)%get_eigenval(i, lambda) 
				if (indata%Vstat==3) then
					write(3,'(1P,2E15.5)') lambda, e_exp
					print'(I3, ES20.9,ES20.9)', ell+i, lambda, e_exp
				else
					write(3,'(2E15.5)') lambda
					print'(I3.3,ES20.9)', i, lambda
				end if
			end do
		close(3)
		!------------------------------------------------------------------------!
	
		!------------------------------------------------------------------------!
		! Calculate the overlap for the basis functions for free H				 !
		! <psi_nl|sin(kr)>														 !
		! Print the overlap values to a file 									 !
		!------------------------------------------------------------------------!
		if (indata%Vstat==1 .and. ell==0) then
			allocate (overlap(N,N))
			do i=1, N
				do j=1, N
					if (indata%OVmeth==1) then
						call eigenstuff(l)%get_eigenval(i, tmp)
						do k=0, Rlength
							tmp_fun(k)=sin(sqrt(2.0*tmp)*rgrid(k))
						end do	
						! Simpsons method for overlaps, n=Rlength, a=0, b=Rmax
						h=indata%stepsize; x1=0d0; x2=0d0
						tmp_b=wavefuns(l)%psi(Rlength,j)*rgrid(Rlength)
						x0=(0)+(tmp_b*tmp_fun(Rlength))
						do k=1, (Rlength-1)
							x=real(k*h)
							tmp=wavefuns(l)%psi(k,j)*rgrid(k)
							if (mod(k,2)==0) then
								x2=x2+(tmp*tmp_fun(k))
							else
								x1=x1+(tmp*tmp_fun(k))
							end if
						end do
						overlap(j,i)=h*(x0+2.0*x2+4.0*x1)/3.0
					! Analytic method 
					else if (indata%OVmeth==0 .and. indata%BasisType==0) then
						tmp=0.0
						do k=1,N
							call eigenstuff(l)%get_eigenval(i, tmp_a)
							x0=sqrt(2.0*tmp_a); tmp_a=0.0
							tmp_b=eigenstuff(l)%eigenvec(k,j)
							call Gbasis(l)%get_nu(k,x)
							x2=(5.d0+real(ell))/2.d0; h=-(x0*x0)/(4.0*x)
							!print*, 'u=', tmp_b
							!print*, 'a=', x2
							!print*, 'b=', 3.0/2.0
							!print*, 'z=', h
							if (h<-700.0) then
								hg=0.d0
							else
								hg = real(CONHYP(complex(x2,0.d0),complex(3.d0/2.d0,0.d0),complex(h,0.d0),0,10),kind=range)
							end if
							!print*, '1F1(h)=',hg
							x1 = gamma((5.d0+real(ell))/2.d0)
							tmp=tmp+(tmp_b*(0.5*x0*x**(-5.0/2.0-ell/2.0)*x1*hg))
							!print*, tmp
							tmp_a=0.0; tmp_b=0.0; x=0.0; x0=0.0; x1=0.0
						end do
						overlap(j,i)=tmp
						print*, '<',j,'|sin(',i,'r)>=',overlap(i,j)
						tmp_a=0.0; tmp_b=0.0; x=0.0; x0=0.0; x1=0.0
						h=0.0
					end if
				end do
				if (indata%verbose) print*, 'E',i,' overlaps',overlap(i,:)
			end do 	
		!------------------------------------------------------------------------!
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
		!------------------------------------------------------------------------!
		
		
		!------------------------------------------------------------------------!
		! Calculate expectation values <r^p> = <psi_nl|r^p|psi_nl>				 !
		!------------------------------------------------------------------------!
		! p_max is the largest absolute power of r needed, set as required
	  	! We only want this calculation for the BOUND states for hydrogenic potential
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
							! Simpsons method for expectation values, n=Rlength, a=0, b=Rmax
							h=indata%stepsize; x1=0d0; x2=0d0
							tmp_b=wavefuns(l)%psi(Rlength,i)
							x0=(0)+(tmp_b**2 *(indata%Rmax**p)*indata%Rmax**2)
							do k=1, (Rlength-1)
								x=real(k*h)
								tmp=wavefuns(l)%psi(k,i)
								if (mod(k,2)==0) then
									x2=x2+(tmp**2*(x**p)*x**2)
								else
									x1=x1+(tmp**2*(x**p)*x**2)
								end if
							end do
							overlap(i,counter)=h*(x0+2.0*x2+4.0*x1)/3.0
							print'(a,I2.2,a,I2.2,a,I2,a,2E15.5)', 'n=',i+ell,&
									' l=',ell,' p=',p,' <r^p>=',overlap(i,counter)
						end if
					end do
				end if
			end do	
		end if	
		!------------------------------------------------------------------------!
		
		!------------------------------------------------------------------------!
		! Plotting the wavefunctions of calculated eigen energies	 			 !
		!------------------------------------------------------------------------!		
		if (indata%wavefunplots) then
	 		if (indata%Vstat==3) then
				do i=1, N
					!do j=0, Rlength
					!	tmp_fun(j)=R(i,ell,1,rgrid(j))
					!end do
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_r_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(rgrid(:)*tmp_fun(:),kind=8),&
							real(rgrid(:),kind=8), real(rgrid(:)*wavefuns(l)%psi(:,i),kind=8),&
							filename=outfile, terminal='png')
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'_',ell,'.png'
					call plot(&
							!real(rgrid(:),kind=8), real(tmp_fun(:),kind=8),&
							real(rgrid(:),kind=8), real(wavefuns(l)%psi(:,i),kind=8),&
							filename=outfile, terminal='png')
				end do
			else if (indata%Vstat==1 .and. ell==0) then
				do i=1, N
					do j=0, Rlength
						call eigenstuff(l)%get_eigenval(i,tmp)
						tmp_a=sqrt(2.0*tmp)
						tmp_fun(j)=sin(tmp_a*rgrid(j))
					end do	
					write(outfile, '(a,I3.3,a,I2.2,a)') 'wavefunction_',i+ell,'',ell,'_freeH.png'
					call plot(real(rgrid(:),kind=range), real(rgrid(:)*wavefuns(l)%psi(:,i),kind=range),&
								real(rgrid(:),kind=range), real(overlap(i,i)*rgrid(:)*wavefuns(l)%psi(:,i),kind=range),&
						   		real(rgrid(:),kind=range),tmp_fun(:),filename=outfile, terminal='png')
				end do
			end if
		end if 
		!------------------------------------------------------------------------!

		!------------------------------------------------------------------------!
		! Destructing the HTVS object so we can reuse them 						 !
		! all the important stuff has been extracted by now						 !
		!------------------------------------------------------------------------!
		call HTVSmats%destruct()

	end do
	!------------------------------------------------------------------------!
	! End of the angular momenta loop				 						 !
	!------------------------------------------------------------------------!
	
	!------------------------------------------------------------------------!
	! Calculation of the oscillator strength and static atomic polarizability!
	!------------------------------------------------------------------------!
	if (indata%OScalc) then
		if (l<2) STOP 'ERROR: Need at least two bases to calc os'
		if (indata%BasisType==0) then
			if (indata%OSmeth==1 .or. indata%OSmeth==3) then
				call cpu_time(t1)
				call os_calc_numerical(Gbasis, energies, wavefuns, rgrid)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for numerical oscillator strength calc', t2-t1, 'seconds'
			end if
			if (indata%OSmeth==0 .or. indata%OSmeth==3) then
				call cpu_time(t1)
				call os_calc_analytic(Gbasis, eigenstuff)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for analytic oscillator strength calc', t2-t1, 'seconds'
			end if
		else 
			if (indata%OSmeth==1 .or. indata%OSmeth==3) then
				call cpu_time(t1)
				call os_calc_numerical(Gbasis_c, energies, wavefuns, rgrid)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for numerical oscillator strength calc', t2-t1, 'seconds'
			end if
			if (indata%OSmeth==0 .or. indata%OSmeth==3) then
				call cpu_time(t1)
				call os_calc_analytic(Gbasis_c, eigenstuff)
				call cpu_time(t2)
				if (indata%time) print*, 'time taken for analytic oscillator strength calc', t2-t1, 'seconds'
			end if
		end if
	end if

end program main
