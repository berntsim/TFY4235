FUNCTION Ur(x)
use parameters
implicit none
real(kind=wp) :: Ur,tmp = 1._wp
real(kind=wp),intent(in) :: x

if (x .gt. 0._wp) then
	if (MODULO(x,tmp) .lt. alpha) then							!This funcion is defined in 
		Ur = MODULO(x,tmp)/(alpha)							!http://web.phys.ntnu.no/~ingves/Teaching/TFY4235/Assignments/TFY4235_Assignment_01.pdf
		return										!where it is a piecewise smooth funcion defined within given limits. 
	else if ((MODULO(x,tmp) .gt. alpha) .and. (MODULO(x,tmp) .lt. tmp)) then		!This function tests which regime one is in, and the returns the appropriate value.
		Ur = (tmp-MODULO(x,tmp))/(tmp-alpha)
		return
	else 
		write (*,*) "Invalid argument, out of range of potential function"
	end if
else if (x .eq. 0._wp) then
	Ur = 0._wp
	return
else
	if ((-tmp .lt. (MODULO(x,tmp)-tmp)) .and. ((MODULO(x,tmp)-tmp) .lt. (alpha-tmp))) then
		Ur = MODULO(x,tmp)/alpha
		return
	else if ((alpha-tmp .lt. (MODULO(x,tmp)-tmp)) .and. ((MODULO(x,tmp)-tmp) .lt. 0._wp)) then
		Ur = (tmp - MODULO(x,tmp))/(tmp-alpha)
		return
	else
		write (*,*) "Invalid argument, out of range of potential function"
	end if
end if
END FUNCTION

