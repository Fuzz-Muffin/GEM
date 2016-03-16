module input_class
	use para
	implicit none 
	private

	type, public :: input
	! These are the input things for GEM method
	   	integer :: l_i, l_f	
		integer :: nfancy, Vstat, BasisType, OSmeth, p_max, OVmeth
		logical :: verbose, basisplots,wavefunplots, overlaps, expect_r, OScalc, time
	   	real(kind=range) :: r1, rN, alpha, mu, lowtol, hightol, Rmax, stepsize, alpha_lag
	! And these are for the Laguerre basis
		integer :: npwave, npdbl, ndouble, ltmax
	   	real(kind=range) :: formcut, regcut, expcut, qmax

	contains	
		procedure :: new => new_input
	end type input
	
contains
	! Read the read_in file and construct an input object
	subroutine new_input(self, infile, iwrite)
		class(input), intent(inout) :: self
		integer, intent(in) :: infile, iwrite
		logical :: ex
		integer :: iostat_data
	  	! infile checks for the existence of data.in
	   	! iwrite =1 then we print the data to the terminal, =0 then not
	   	
		if(infile .eq. 10) then
       		inquire(file='data.in',exist=ex)
       		if(.not. ex) STOP 'ERROR: File data.in does not exists'
        	open(infile,file='data.in',iostat=iostat_data)
       		if(iostat_data.ne.0) STOP 'ERROR: input_class cannot open file data.in'
   		endif

		!------------------------------------------------------------------------!
		! Here we read the values from data.in and assigns values to input object!
		!------------------------------------------------------------------------!
		read(infile,*) self%verbose, self%time
		read(infile,*) self%npwave, self%npdbl, self%ndouble, self%ltmax, self%qmax
		read(infile,*) self%formcut, self%regcut, self%expcut
		read(infile,*) self%basisplots, self%wavefunplots
		read(infile,*) self%overlaps, self%OVmeth
		read(infile,*) self%expect_r, self%p_max
		read(infile,*) self%OScalc, self%OSmeth 
		read(infile,*) self%Rmax, self%stepsize
		read(infile,*) self%BasisType
		read(infile,*) self%l_i, self%l_f
		read(infile,*) self%alpha_lag
		read(infile,*) self%r1, self%rN, self%nfancy
		read(infile,*) self%alpha
		read(infile,*) self%mu 
		read(infile,*) self%lowtol, self%hightol
		read(infile,*) self%Vstat
		!------------------------------------------------------------------------!
		! For iwrite=1 print all of the input to terminal 						 !
		!------------------------------------------------------------------------!
		if (iwrite==1) then
			print*, 'Verbose mode ', self%verbose, 'Timing code?', self%time
			print*, self%npwave, self%npdbl, self%ndouble, self%ltmax, self%qmax,'npwave, npdbl, ndouble, ltmax, qmax'
			print*, self%formcut, self%regcut, self%expcut, 'formcut, regcut, expcut'
			print*, 'Save plots? ', self%basisplots, self%wavefunplots
			print*, 'Print calculated overlaps? (only for free H)', self%overlaps
			if (self%OVmeth==0) then
				print*, 'Using analytic overlap subroutine'
			else if (self%OVmeth==1) then
				print*, 'Using numerical overlap strength subroutine'
			end if
			print*, 'calculate expectation values? (only for H atom)', self%expect_r
			print*, 'For <n|r^p|n> use abs val of p=', self%p_max 
			print*, 'Calculate oscillator strengths for 1s->np', self%OScalc
			if (self%OSmeth==0) then
				print*, 'Using analytic oscillator strength subroutine'
			else if (self%OSmeth==1) then
				print*, 'Using numerical oscillator strength subroutine'
			else
				print*, 'Using both analytic and numerical OS subroutines'
			end if
			print*, 'Rmax, stepsize:', self%Rmax, self%stepsize
			if (self%BasisType==2) then
				print*, 'Using the Laguerre basis'
				print*, 'Parameters, alpha_lag: ', self%alpha_lag
			else if (self%BasisType==1) then
				print*, 'Using the complex range Gaussian basis'
				print*, 'Parameters, r1, rN, N: ',self%r1, self%rN, self%nfancy  
				print*, 'Parameters, alpha: ', self%alpha
			else if (self%BasisType==0) then
				print*, 'Using the regular Gaussian basis'
				print*, 'Parameters, r1, rN, N: ',self%r1, self%rN, self%nfancy  
			end if
			print*,'Angular momentum, calc from l=', self%l_i, 'to l=', self%l_f
			print*, 'Reduced mass mu = ', self%mu
			print*, 'lowtol, hightol: ', self%lowtol, self%hightol
			print*, 'Vstat = ', self%Vstat
		end if 
		!------------------------------------------------------------------------!

	end subroutine new_input

end module input_class
