PROGRAM main
use parameters
use diffusion
use waveparameters
use wave
IMPLICIT NONE

!----------------------Diffusion equation section------------

!call absorbing_boundaries()
!call analytical_solution_absorbing_boundaries()

!call reflecting_boundaries()
!call analytical_solution_reflecting_boundaries()

!call absorbing_boundaries_diffusivity()
!call reflecting_boundaries_diffusivity()

!call plot_gnuplot()
!----------------------Wave equation section-----------------

call wave_equation()


END PROGRAM
