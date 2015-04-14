MODULE diffusion
use parameters
IMPLICIT NONE
contains

SUBROUTINE plot_gnuplot()
        Call SYSTEM('gnuplot -p "plotfortran.p"')
END SUBROUTINE


SUBROUTINE optimized_theta_scheme_absorbing_boundaries(vector_)
real(kind=wp), dimension(N), intent(inout) :: vector_
real(kind=wp), dimension(N) :: tmp_
integer :: i_
integer :: info_

do i_=1,N                                       !Initializing the diagonals
        ADiag(i_) = 1._wp+2._wp*theta*K
        BDiag(i_) = 1._wp-2._wp*(1._wp-theta)*K
end do

do i_ =1,N-1                                    !Initializing the upper and lower diagonals
        AUpDiag(i_) = -theta*k
        ALowDiag(i_) = -theta*k
        BUpDiag(i_) = (1._wp-theta)*K
        BLowDiag(i_) = (1._wp-theta)*K
end do

do i_=1,N                                       !This section is a matrix multiplication for a tridiagonal matrix with a vector, which takes less computating time than matmul
        if (i_ .eq. 1) then
                tmp_(i_) = BDiag(i_)*vector_(i_) + BUpDiag(i_)*vector_(i_+1)
        else if (i_ .eq. N) then
                tmp_(i_) = BLowDiag(i_-1)*vector_(i_-1) + BDiag(i_)*vector_(i_)
        else 
                tmp_(i_) = BLowDiag(i_-1)*vector_(i_-1) + BDiag(i_)*vector_(i_) + BUpDiag(i_)*vector_(i_+1)
        end if
end do
vector_=tmp_
call dgtsv(N,1,ALowDiag,ADiag,AUpDiag,vector_,N,info_)          !Using the lapack routine for tridiagonal Ax = b problems
END SUBROUTINE


 

SUBROUTINE absorbing_boundaries()
dx = L/(N+1)            !Setting space-step
dt = r*dx**2            !Setting timestep for euler schemes
Nt =nint(t/dt)   !Deciding number of iterations to do forward in time
K = D*dt/(dx*dx)        !Declaring the reduced diffusion coefficient
uVector = 0._wp         !Setting all elements of the concentration vector to zero
uVector(int((xZero/L)*(N+1))) = uTilde/dx        !Implementing the starting concentration as a diracs delta-function (approximation, obviously)
print*, "Nt = ",Nt
print*, "dt = ",dt

open (unit=out_unit1, file="startConcentration.dat")            !Writing inital concentration til file, so that it can be plotted using gnuplot
write(out_unit1,*) 0._wp
do i = 1,N
        write(out_unit1,*) uVector(i)
end do
open (unit=out_unit, file="testConcentration.dat")              !Opening to prepare file for writing 
write(out_unit1,*) 0._wp
write(out_unit,*) 0._wp, 0._wp
do i = 1,Nt                                                             !Iterating forward in time
        call optimized_theta_scheme_absorbing_boundaries(uVector)       !Using the theta-scheme implementet above to do the actual calculations
end do
do j = 1,N                                                              !Writing the final concenctration to file
        write(out_unit,*) uVector(j), j*dx               
end do
dummy = Nt*dt                                                           !Calculating what the time is at the end of the iterations, so that the file with the final plot makes sense
print*, "t = ", dummy 
write(out_unit,*) 0._wp,1._wp                                          
END SUBROUTINE

SUBROUTINE optimized_theta_scheme_reflecting_boundaries(vector_RB)
real(kind=wp), dimension(N), intent(inout) :: vector_RB
real(kind=wp), dimension(N) :: tmp_RB
integer :: i_RB
integer :: info_RB

do i_RB=1,N                                       !Initializing the diagonals
        ADiag(i_RB) = 1._wp+2._wp*theta*K
        BDiag(i_RB) = 1._wp-2._wp*(1._wp-theta)*K
end do

do i_RB =1,N-1                                    !Initializing the upper and lower diagonals
        AUpDiag(i_RB) = -theta*k
        ALowDiag(i_RB) = -theta*k
        BUpDiag(i_RB) = (1._wp-theta)*K
        BLowDiag(i_RB) = (1._wp-theta)*K
end do

AUpDiag(1) =    2*AUpDiag(1)                     !Fixing the reflective boundaries, according to u(i+1) = u(i-1)
ALowDiag(N-1) = 2*ALowDiag(N-1)
BUpDiag(1) =    2*BUpDiag(1)
BLowDiag(N-1) = 2*BLowDiag(N-1)

do i_RB=1,N                                       !This section is a matrix multiplication for a tridiagonal matrix with a vector, which takes less computating time than matmul
        if (i_RB .eq. 1) then
                tmp_RB(i_RB) = BDiag(i_RB)*vector_RB(i_RB) + BUpDiag(i_RB)*vector_RB(i_RB+1)
        else if (i_RB .eq. N) then
                tmp_RB(i_RB) = BLowDiag(i_RB-1)*vector_RB(i_RB-1) + BDiag(i_RB)*vector_RB(i_RB)
        else 
                tmp_RB(i_RB) = BLowDiag(i_RB-1)*vector_RB(i_RB-1) + BDiag(i_RB)*vector_RB(i_RB) + BUpDiag(i_RB)*vector_RB(i_RB+1)
        end if
end do
vector_RB=tmp_RB
call dgtsv(N,1,ALowDiag,ADiag,AUpDiag,vector_RB,N,info_RB)          !Using the lapack routine for tridiagonal Ax = b problems
END SUBROUTINE


SUBROUTINE reflecting_boundaries()
dx = L/(N-1)            !Setting space-step
dt = r*dx**2            !Setting timestep for euler schemes
Nt = nint(t/dt)         !Deciding number of iterations to do forward in time
K = D*dt/(dx*dx)        !Declaring the reduced diffusion coefficient
uVector = 0._wp         !Setting all elements of the concentration vector to zero
uVector(int((xZero/L)*(N+1))) = uTilde/dx        !Implementing the starting concentration as a diracs delta-function (approximation, obviously)

open (unit=out_unit3, file="refBound.dat")            !Writing inital concentration til file, so that it can be plotted using gnuplot

do i = 1,Nt                                                             !Iterating forward in time
        call optimized_theta_scheme_reflecting_boundaries(uVector)       !Using the theta-scheme implementet above to do the actual calculations
end do
do j = 1,N                                                              !Writing the final concenctration to file
        if (j .eq. 1) then
                write(out_unit3,*) uVector(j),j*dx                           !Making sure that the end points are included. This is done in accourdance do d/dx(u) = 0 at the bordes
                write(out_unit3,*) uVector(j),j*dx               
        else if (j .eq. N) then
                write(out_unit3,*) uVector(j),j*dx                           !Including the second end point.
                write(out_unit3,*) uVector(j),j*dx
        else
                write(out_unit3,*) uVector(j),j*dx                           !Writing the outher points normally
        end if
end do
END SUBROUTINE



SUBROUTINE vnAbsorbing(x_,n_,retVar_)
real(kind=wp),intent(out) :: retVar_
real(kind=wp),intent(in)  :: x_
integer,intent(in) :: n_
if (n_ .eq. 0) then
        retVar_ = 0._wp
else
        retVar_ = dsqrt(2._wp/L)*dsin(n_*pi*x_/L)
end if
END SUBROUTINE

SUBROUTINE vnReflecting(x_RB,n_RB,retVar_RB)
real(kind=wp), intent(out)      :: retVar_RB
real(kind=wp), intent(in)       :: x_RB
integer, intent(in)             :: n_RB

if(n_RB .eq. 0) then
        retVar_RB = dsqrt(1._wp/L)
else
        retVar_RB = dsqrt(2._wp/L)*dcos(n_RB*pi*x_RB/L)       
end if
END SUBROUTINE


SUBROUTINE analytical_solution_absorbing_boundaries()
real(kind=wp) :: sum_= 0._wp
real(kind=wp) :: value_, value_xZero 
real(kind=wp) :: xVal = 0._wp
open (unit=out_unit2, file="AAB.dat")
do i = 1,xGridsize
        sum_= 0._wp
        do j= 0,nAnalyticalMax
                call vnAbsorbing(xVal,j,value_) 
                call vnAbsorbing(xZero,j,value_xZero)
                sum_ = sum_ + dexp(-(j*pi/L)**2._wp*D*t)*value_xZero*value_ 
        end do
        xVal = xVal + 1._wp/xGridSize
        write(out_unit2,*) uTilde*sum_, xVal
end do
END SUBROUTINE

SUBROUTINE analytical_solution_reflecting_boundaries()
real(kind=wp)   :: sum_RB = 0._wp
real(kind=wp)   :: value_RB
real(kind=wp)   :: value_xZeroRB
real(kind=wp)   :: xValRB = 0._wp
open (unit=out_unit4, file="ARB.dat")
do i = 1,xGridSize
        sum_RB = 0._wp
        do j = 0,nAnalyticalMax
                call vnReflecting(xValRB,j,value_RB)
                call vnReflecting(xZero,j,value_xZeroRB)
                sum_RB = sum_RB + dexp(-(j*pi/L)**2*D*t)*value_xZeroRB*value_RB
        end do
        xValRB = xValRB + 1._wp/xGridSize
        write(out_unit4,*) uTilde*sum_RB,xValRB
end do
END SUBROUTINE

SUBROUTINE diffusion_profile_heavyside(x_DP,retVal1_DP,retVal2_DP)
real(kind=wp), intent(in)               :: x_DP
real(kind=wp), intent(inout)              :: retVal1_DP,retVal2_DP

if ((nint(x_DP*rcc) .gt. nint((xZero-disp)*rcc)) .or. (nint(x_DP*rcc) .eq.nint(( xZero-disp)*rcc))) then
        retVal1_DP = DPluss  
else
        retVal1_DP = DMinus
end if

!if (nint(x_DP*rcc) .eq. nint((xZero-disp)*rcc))then
!        retVal2_DP = DPluss-DMinus
!else
        retVal2_DP = 0._wp
!end if
END SUBROUTINE

SUBROUTINE diffusion_profile_constant(retVal1_DP,retVal2_DP)
real(kind=wp), intent(out)              :: retVal1_DP,retVal2_DP
retVal1_DP = D
retVal2_DP = 0._wp
END SUBROUTINE

SUBROUTINE diffusion_profile_linear(x_DP,retVal1_DP,retVal2_DP)
real(kind=wp), intent(in)               :: x_DP
real(kind=wp), intent(out)              :: retVal1_DP,retVal2_DP
retVal1_DP = x_DP +0.01
retVal2_DP = 1._wp

END SUBROUTINE

SUBROUTINE CN_scheme_diffusivity(vector_AB)
real(kind=wp), dimension(N), intent(inout) :: vector_AB
real(kind=wp), dimension(N) :: tmp_AB
real(kind=wp), dimension(N) :: diffusivity, beta
real(kind=wp) :: xVal_AB = 0._wp
integer :: i_AB
integer :: info_AB


do i_AB = 1,N
        call diffusion_profile_heavyside(xVal_AB,diffusivity(i_AB),beta(i_AB))
        !call diffusion_profile_constant(diffusivity(i_AB),beta(i_AB))
        !call diffusion_profile_linear(xVal_AB,diffusivity(i_AB),beta(i_AB))
        xVal_AB = xVal_AB + 1._wp/N
        
end do
xVal_AB = 0._wp
 
do i_AB=1,N                                       !Initializing the diagonals
        ADiag(i_AB) = 2._wp + 4._wp*diffusivity(i_AB)*r*theta
        BDiag(i_AB) = -4._wp*diffusivity(i_AB)*(1._wp-theta)*r + 2._wp 
end do

do i_AB =1,N-1                                    !Initializing the upper and lower diagonals
        BUpDiag(i_AB) = r*(1._wp-theta)*(dx*beta(i_AB) + 2._wp*diffusivity(i_AB))
        ALowDiag(i_AB) = r*theta*(dx*beta(i_AB+1) - 2._wp*diffusivity(i_AB+1))
        AUpDiag(i_AB) = -r*theta*(dx*beta(i_AB) + 2._wp*diffusivity(i_AB))
        BLowDiag(i_AB) = r*(1._wp-theta)*(-dx*beta(i_AB+1) + 2._wp*diffusivity(i_AB+1))
end do


do i_AB=1,N                                       !This section is a matrix multiplication for a tridiagonal matrix with a vector, which takes less computating time than matmul
        if (i_AB .eq. 1) then
                tmp_AB(i_AB) = BDiag(i_AB)*vector_AB(i_AB) + BUpDiag(i_AB)*vector_AB(i_AB+1)
        else if (i_AB .eq. N) then
                tmp_AB(i_AB) = BLowDiag(i_AB-1)*vector_AB(i_AB-1) + BDiag(i_AB)*vector_AB(i_AB)
        else 
                tmp_AB(i_AB) = BLowDiag(i_AB-1)*vector_AB(i_AB-1) + BDiag(i_AB)*vector_AB(i_AB) + BUpDiag(i_AB)*vector_AB(i_AB+1)
        end if
end do
vector_AB=tmp_AB
call dgtsv(N,1,ALowDiag,ADiag,AUpDiag,vector_AB,N,info_AB)          !Using the lapack routine for tridiagonal Ax = b problems
END SUBROUTINE


SUBROUTINE absorbing_boundaries_diffusivity()
dx = L/(N+1)            !Setting space-step
dt = r*dx**2            !Setting timestep for euler schemes
Nt = nint(t/dt)         !Deciding number of iterations to do forward in time
KDiff = dt/(2*dx*dx)        !Declaring the reduced diffusion coefficient
uVector = 0._wp         !Setting all elements of the concentration vector to zero
uVector(int((xZero/L)*(N+1))) = uTilde/dx        !Implementing the starting concentration as a diracs delta-function (approximation, obviously)

open (unit=out_unit5, file="absBoundDiffusivity.dat")            !Writing inital concentration til file, so that it can be plotted using gnuplot
write(out_unit5,*) 0._wp, 0._wp
do i = 1,Nt                                                             !Iterating forward in time
        call CN_scheme_diffusivity(uVector)       !Using the theta-scheme implementet above to do the actual calculations
end do

do j = 1,N                                                              !Writing the final concenctration to file
        write(out_unit5,*) uVector(j), j*dx                
end do
write(out_unit5,*) 0._wp, 1._wp
END SUBROUTINE



SUBROUTINE CN_scheme_diffusivity_ref(vector_RB)
real(kind=wp), dimension(N), intent(inout) :: vector_RB
real(kind=wp), dimension(N) :: tmp_RB
real(kind=wp), dimension(N) :: diffusivity, beta
real(kind=wp) :: xVal_RB = 0._wp
integer :: i_RB
integer :: info_RB


do i_RB = 1,N
        call diffusion_profile_heavyside(xVal_RB,diffusivity(i_RB),beta(i_RB))
        !call diffusion_profile_constant(diffusivity(i_RB),beta(i_RB))
        !call diffusion_profile_linear(xVal_RB,diffusivity(i_RB),beta(i_RB))
        xVal_RB = xVal_RB + 1._wp/N
        
end do
xVal_RB = 0._wp
 
do i_RB=1,N                                       !Initializing the diagonals
        ADiag(i_RB) = 2._wp + 4._wp*diffusivity(i_RB)*r*theta
        BDiag(i_RB) = -4._wp*diffusivity(i_RB)*(1._wp-theta)*r + 2._wp 
end do

do i_RB =1,N-1                                    !Initializing the upper and lower diagonals
        BUpDiag(i_RB) = r*(1._wp-theta)*(dx*beta(i_RB) + 2._wp*diffusivity(i_RB))
        ALowDiag(i_RB) = r*theta*(dx*beta(i_RB+1) - 2._wp*diffusivity(i_RB+1))
        AUpDiag(i_RB) = -r*theta*(dx*beta(i_RB) + 2._wp*diffusivity(i_RB))
        BLowDiag(i_RB) = r*(1._wp-theta)*(-dx*beta(i_RB+1) + 2._wp*diffusivity(i_RB+1))
end do

AUpDiag(1) =    -4._wp*r*theta*diffusivity(1)                    !Fixing the reflective boundaries, according to u(i+1) = u(i-1)
ALowDiag(N-1) = -4._wp*r*theta*diffusivity(N)
BUpDiag(1) =     4._wp*r*diffusivity(1)*(1._wp-theta)
BLowDiag(N-1) =  4._wp*r*diffusivity(N)*(1._wp-theta)

do i_RB=1,N                                       !This section is a matrix multiplication for a tridiagonal matrix with a vector, which takes less computating time than matmul
        if (i_RB .eq. 1) then
                tmp_RB(i_RB) = BDiag(i_RB)*vector_RB(i_RB) + BUpDiag(i_RB)*vector_RB(i_RB+1)
        else if (i_RB .eq. N) then
                tmp_RB(i_RB) = BLowDiag(i_RB-1)*vector_RB(i_RB-1) + BDiag(i_RB)*vector_RB(i_RB)
        else 
                tmp_RB(i_RB) = BLowDiag(i_RB-1)*vector_RB(i_RB-1) + BDiag(i_RB)*vector_RB(i_RB) + BUpDiag(i_RB)*vector_RB(i_RB+1)
        end if
end do
vector_RB=tmp_RB
call dgtsv(N,1,ALowDiag,ADiag,AUpDiag,vector_RB,N,info_RB)          !Using the lapack routine for tridiagonal Ax = b problems
END SUBROUTINE




SUBROUTINE reflecting_boundaries_diffusivity()
dx = L/(N-1)            !Setting space-step
dt = r*dx**2            !Setting timestep for euler schemes
Nt = nint(t/dt)         !Deciding number of iterations to do forward in time
KDiff = dt/(2*dx*dx)        !Declaring the reduced diffusion coefficient
uVector = 0._wp         !Setting all elements of the concentration vector to zero
uVector(int((xZero/L)*(N+1))) = uTilde/dx        !Implementing the starting concentration as a diracs delta-function (approximation, obviously)

open (unit=out_unit7, file="refBoundDiffusivity.dat")            !Writing inital concentration til file, so that it can be plotted using gnuplot
!write(out_unit7,*) 0._wp
do i = 1,Nt                                                             !Iterating forward in time
        call CN_scheme_diffusivity_ref(uVector)       !Using the theta-scheme implementet above to do the actual calculations
end do

do j = 1,N                                                              !Writing the final concenctration to file
         if (j .eq. 1) then
                write(out_unit7,*) uVector(j),j*dx                           !Making sure that the end points are included. This is done in accourdance do d/dx(u) = 0 at the bordes
               ! write(out_unit7,*) uVector(j),j*dx               
        else if (j .eq. N) then
                !write(out_unit7,*) uVector(j),j*dx                           !Including the second end point.
                write(out_unit7,*) uVector(j),j*dx
        else
                write(out_unit7,*) uVector(j),j*dx                           !Writing the outher points normally
        end if
end do
!write(out_unit5,*) 0._wp
!print*, uVector

END SUBROUTINE


END MODULE
