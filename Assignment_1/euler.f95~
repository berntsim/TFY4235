FUNCTION euler(xIn,tIn,xi)
use parameters
use random
implicit none
real(kind=wp),intent(in)	:: xIn,tIn,xi
real(kind=wp)			:: euler,force

euler = xIn + force(xIn,tIn)*deltaT + const*xi			!See http://web.phys.ntnu.no/~ingves/Teaching/TFY4235/Assignments/TFY4235_Assignment_01.pdf for
										!where this comes from. 
return
END FUNCTION
