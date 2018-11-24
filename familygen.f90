
       INCLUDE 'nbp_potential.f90'
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   nbp :    N- Body Problem continuation with general scheme for periodic orbits
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
! 
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 
!  
      USE nbp
      IMPLICIT NONE
      REAL(8), DIMENSION(NDIM),INTENT(in)::U
      REAL(8), DIMENSION(NDIM),INTENT(out)::F
      REAL(8), DIMENSION(*),INTENT(in):: PAR
      INTEGER, INTENT(in)::NDIM,IJAC
      REAL(8), DIMENSION(*)::DFDU,DFDP
      INTEGER,DIMENSION(*)::ICP
      REAL(8), DIMENSION(3)::mass
      !      The power is par6
      mass=(/PAR(1),PAR(2),PAR(3)/)
      CALL nbpeq(U,mass,F,PAR(6))
      !       PAR(10), PAR(9),PAR(8),PAR(7) are perturbing parameters
      ! First to be perturbated is the vector field 
      CALL perturbate_Hamiltonian(F,U,mass,PAR(6),PAR(7))
      CALL perturbate_linear_momentum(F,1,PAR(10))
      CALL perturbate_linear_momentum(F,2,PAR(9))
      CALL perturbate_angular_momentum(F,U,1,2,PAR(8))
       F=PAR(11)*F

! 
      RETURN 
       END 
!---------------------------------------------------------------------- 
! 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 
! 
      IMPLICIT NONE
      INTEGER,INTENT(in)::NDIM
      REAL(8),DIMENSION (NDIM),INTENT(in out)::U
      REAL(8),DIMENSION (*),INTENT(in out)::PAR
      REAL(8),INTENT(in)::T 
!
      PAR(1)=1.D0  ! masses
      PAR(2)=1.D0
      PAR(3)=1.D0  
      PAR(6)=1.D0 ! power
      ! perturbation parameters
      PAR(7)=0.D0
      PAR(8)=0.D0
      PAR(9)=0.D0
      PAR(10)=0.D0
      ! Period
       PAR(11)=3.14159265358979323846D0*2.D0
!
! 
      RETURN 
      END 
!---------------------------------------------------------------------- 
! 
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 
! 
      IMPLICIT NONE
      INTEGER, INTENT(in):: NDIM,NBC,IJAC
      REAL(8),DIMENSION(*),INTENT(in)::PAR
      REAL(8),DIMENSION(NDIM),INTENT(in)::U0,U1
      REAL(8),DIMENSION(NBC),INTENT(out)::FB
      INTEGER,DIMENSION(*)::ICP  
      REAL(8),DIMENSION(*)::DBC 
      ! Auxiliary variable for loop 
      INTEGER::I
      !
      FB(1:NDIM)=(/(U1(I)-U0(I),I=1,NDIM)/)
! 
      RETURN 
      END 
!---------------------------------------------------------------------- 
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
        USE nbp
        IMPLICIT NONE
        INTEGER,INTENT(in)::NDIM,NINT,IJAC !Dimension of the ODE system 
        REAL(8),DIMENSION(*),INTENT(in)::PAR,DINT
        INTEGER,DIMENSION(*),INTENT(in):: ICP
        REAL(8),DIMENSION(NDIM),INTENT(in)::U,UOLD,UPOLD,UDOT
        REAL(8),DIMENSION(NINT),INTENT(out)::FI
        FI(1)=reduce_linear_momentum(U-UOLD,1)
        FI(2)=reduce_linear_momentum(U-UOLD,2)
        FI(3)=reduce_angular_momentum(U,U-UOLD,1,2)
        FI(4)=reduce_phase(U,U-UOLD,(/PAR(1),PAR(2),PAR(3)/),PAR(6))
      RETURN 
      END 
!
      SUBROUTINE FOPT 
      RETURN 
      END 
! 
      SUBROUTINE PVLS(NDIM,U,PAR)
        USE nbp
      IMPLICIT NONE
      INTEGER,INTENT(in)::NDIM
      REAL(8),DIMENSION(*),INTENT(in)::U
      REAL(8),DIMENSION(*),INTENT(in out)::PAR
      REAL(8)::GETP
      REAL(8),DIMENSION(totaldim)::ic
      INTEGER::I
      
      ic=(/(GETP('BV0',I,U) ,I=1,totaldim)/)
      PAR(4)=nbp_Hamiltonian(ic,(/PAR(1),PAR(2),PAR(3)/),PAR(6))
      PAR(5)=getp('BIF',1,U)
      RETURN 
      END 
!---------------------------------------------------------------------- 
