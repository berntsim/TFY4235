FUNCTION Pu(U)
use parameters
implicit none
!real(kind=wp) :: KbT,deltaU
real(kind=wp),intent(in) :: U
real(kind=wp) ::Pu

Pu = dexp(-(U/KbT))/(KbT*(1._wp - dexp(-deltaU/KbT)))

END FUNCTION
