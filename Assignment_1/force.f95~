FUNCTION force(x,t)
use parameters
implicit none
real(kind=wp)		:: force,tmp = 1._wp
real(kind=wp),intent(in):: x,t


dc = MODULO(t/(tau),tmp)
if (flashing .eqv. .true.) then

	if (((dc .gt. 0._wp) .or. (dc .eq. 0._wp)) .and. (dc .lt. 3._wp/4._wp)) then
		force = 0
		return
	else if ((dc .gt. 3._wp/4._wp) .and. (dc .lt. 1)) then
		if (((MODULO(x,tmp) .gt. 0._wp) .or. (MODULO(x,tmp) .eq. 0._wp)) .and. (MODULO(x,tmp) .lt. alpha)) then
			force = -1._wp/(alpha)
			!print*,"jonas"
			return
		else if ((MODULO(x,tmp) .gt. alpha) .and. (MODULO(x,tmp) .lt. 1._wp)) then
			force = (1._wp)/((1._wp-alpha))
			!print*, "Per"
			return
		end if
	else
		write (*,*) "Invalid argument, out of range of force function"
		print*, "decimalpart = ", dc
		print*, "tau = ",tau
		print*, "x = " ,x
		print*, 
	end if
else
	if (((MODULO(x,tmp) .gt. 0._wp) .or. (MODULO(x,tmp) .eq. 0._wp)) .and. (MODULO(x,tmp) .lt. alpha)) then
		force = -1._wp/(alpha)
		return
	else if ((MODULO(x,tmp) .gt. alpha) .and. (MODULO(x,tmp) .lt. 1._wp)) then
		force = (1._wp)/(1._wp-alpha)
		return		
	end if
end if

if (ensemble .eqv. .true.) then
force = 0
end if
END FUNCTION

