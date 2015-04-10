FUNCTION potential(x,t)
use parameters
implicit none
real(kind=wp) :: potential, Ur, f
real(kind=wp),intent(in) :: x,t

potential = Ur(x)!*f(t)
return
END FUNCTION
