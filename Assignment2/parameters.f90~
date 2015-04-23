MODULE parameters
IMPLICIT NONE
integer, parameter              :: wp=kind(0.d0)        !Set working precission
integer, parameter              :: N = 99              !Setting number of iterations for the program to run
integer                         :: Nt = 0               !Setting how long the program will iterate forward in time
real(kind=wp), parameter        :: r = 0.25_wp          !For the stability of euler
integer                         :: i = 0                !iterator
integer                         :: j = 0                !iterator
real(kind=wp), parameter        :: t = 5.0_wp              !Setting simulation time
real(kind=wp)                   :: dt = 0               !Time iteration step
real(kind=wp)                   :: dx = 0               !Spatial iteration step
real(kind=wp)                   :: D = 0.01_wp           !Diffusion constant
real(kind=wp), parameter        :: L = 1.0_wp           !L = b-a, i.e the spatial length of the system
real(kind=wp), dimension(N,N)   :: A                    !Initializing the matrix A
real(kind=wp), dimension(N,N)   :: B                    !Initializing the matrix B
real(kind=wp), dimension(N)     :: ADiag,BDiag          !Diagonal
real(kind=wp), dimension(N-1)   :: AUpDiag,BUpDiag      !Upper diagonal
real(kind=wp), dimension(N-1)   :: ALowDiag,BLowDiag    !Lower diagonal
real(kind=wp), dimension(N)     :: uVector              !Initial vector for the concentration
real(kind=wp), parameter        :: uTilde = 1._wp       !Characteristic or typical concentration
real(kind=wp), parameter        :: xZero = 0.5_wp       !Location of the starting concentration
real(kind=wp)                   :: K = 0                !Initializing the reduced diffusion constant
real(kind=wp), parameter        :: theta = 0.5_wp       !Making sure we use the Crank-Nicolson scheme. Set to 1 for Euler.          
real(kind=wp)                   :: KDiff = 0._wp        !Declaring the constant used in the diffusive matrices
real(kind=wp), parameter        :: DPluss = 0.1         !Declaring the piecewise constant diffusion coefficients, Heavyside-profile 
real(kind=wp), parameter        :: DMinus = 0.1         
integer, parameter              :: rcc=10**8            !Setting a number for maintaining precission when comparing real(kind=wp) numbers
real(kind=wp)                   :: xValues
real(kind=wp),parameter         :: disp = 0             !Setting displacement for the step         

real(kind=wp) :: dummy = 0

!real(kind=wp)                   :: plotTime = 0.1         !Determining at what time you want to plot the concentration
integer, parameter              :: out_unit = 29        !To save data to file
integer, parameter              :: out_unit1 = 28       !the same
integer, parameter              :: out_unit2 = 27       !For saving the analytical absorbing boundaries
integer, parameter              :: out_unit3 = 26       !For saving the reflective boundary solution
integer, parameter              :: out_unit4 = 25       !For saving the analytical reflecting boundaries
integer, parameter              :: out_unit5 = 24       !For saving the absorbing boundaries problem with spatially dependent diffusivity
!integer, parameter              :: out_unit6 = 23       !For saving x-steps
integer, parameter              :: out_unit7 = 22       !For saving reflective boundaries, spatially dependent diffusion coefficient

!-----------------------Section for analytical solutions -------------------------------------------
integer, parameter                :: nAnalyticalMax = 10**4
integer                           :: nAnalytical = 0
integer, parameter                :: xGridSize = 100
real(kind=wp), dimension(xGridSize)     :: xAnalytical
real(kind=wp), parameter                :: pi = 3.141592653589793238_wp



END MODULE
