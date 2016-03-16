module para
	use toms707
	implicit none
	integer, parameter:: range=SELECTED_REAL_KIND(p=15,r=307)
	real(kind=range), parameter:: epsilon = 1.0E-9_range
	real(kind=range), parameter:: pi   = atan(1.d0)*4.d0
	real(kind=range), parameter:: exptol = 200.00
	real(kind=range), parameter:: tol	= 1d-60
	! Defining the imaginary unit
	complex(kind=range), parameter:: imag = (0.0,1.0)
	! By setting hbar to unity, we are implicitly using atomic units 
	real(kind=range), parameter:: hbar = 1

contains    

	integer function fact(x) result(tot)
		implicit none
		integer, intent(in) :: x
		integer :: i
		tot=x
		if (tot/=0) then
			do i=(x-1),1,-1
				tot=tot*i
			end do
		else
			tot=1
		end if
	end function fact

	integer function doublefact(x) result(tot)
		implicit none
		integer, intent(in) :: x
		integer :: i
		if (x/=0) then
			tot=1
			if ((mod(x,2)/=0)) then
				do i=1,((x+1)/2)
					tot=tot*(2*i-1)
				end do 
			else
				do i=1,(x/2)
					tot=tot*2*i
				end do
			end if
		else
			tot=1
		end if
	end function doublefact
	
	real(kind=range) function sinkr(e, r) 
		implicit none
		real(kind=range), intent(in) :: r, e
		real(kind=range) :: k
		k=sqrt(2*e)
		sinkr=sin(k*r)
	end function sinkr

!	real(kind=range) function R(n,l,Z,x)
!		implicit none
!		real(kind=range), intent(in) :: x
!		real(kind=range) :: rho, tmp_a, tmp_b, a,b
!		integer, intent(in) :: n,l,Z
!		rho=((2.0*Z)/real(n))*x
!		tmp_a=1.0/(real(fact(2*l+1)))*sqrt(((2.0*Z)/(real(n)))**3*&
!				(real(fact(n+l)))/(2.0*real(n*fact(n-l-1))))
!		a=real(l+1-n,kind=range); b=real(2*l+2,kind=range)
!		if () then 
!			R=0.0
!		else
!			tmp_b=real(conhyp(complex(a,0.0),complex(b,0.0),complex(rho,0.0),0,10),kind=range)
!			R=tmp_a*exp(-rho/2.0)*rho**l*tmp_b
!		end if
!	end function R
		

end module para
