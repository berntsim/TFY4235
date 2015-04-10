FUNCTION f(t)
use parameters
implicit none
real(kind=wp)	:: f,tmp = 1._wp
real(kind=wp),intent(in) :: t 

if (((MODULO(t,tmp) .gt. 0._wp) .or. (MODULO(t,tmp) .eq. 0._wp)) .and. (MODULO(t,tmp) .lt. 3._wp*tau/4._wp)) then			!should be limit 3*pi/2
	f = 0._wp
	return
else if ((MODULO(t,tmp) .gt. 3._wp*tau/4._wp) .and. (MODULO(t,tmp) .lt. tau)) then			! should be limit 3*pi/2 and 2*pi || or 3*tau/4 and tau
	f = 1._wp
	return
else 
	write (*,*) "Invalid argument, out of range of f function"
end if

END FUNCTION
