PROGRAM MAIN 
  USE nbp 
  IMPLICIT NONE
  REAL(8),DIMENSION(totaldim):: x0
  REAL(8),DIMENSION(nbody)::mass
  REAL(8)::PER
  INCLUDE 'initial_condition.h'
  CALL print_initial_solution(x0,mass,1.D0,PER,1000)
  RETURN
END PROGRAM MAIN

