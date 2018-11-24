
       INCLUDE 'nbp_potential.f90'
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   nbp :    N- Body Problem
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
!      The power is par6
      CALL nbpeq(U,(/PAR(1),PAR(2),PAR(3)/),F,PAR(6))
!       PAR(10) es un parametro a perturbar 
       CALL perturbate_angular_momentum(F,U,1,2,PAR(10))
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
       PAR(10)=0.D0 ! perturbation parameter 
       PAR(11)=3.14159265358979323846D0
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
!      
       FB(1)=U0(1)
       FB(2)=U0(2)
       FB(3)=U0(3)+U0(5) 
       FB(4)=U0(4)+U0(6)
       FB(5)=U0(9)-U0(11)
       FB(6)=U0(10)-U0(12)
       FB(7)=U1(1)
       FB(8)=U1(2)
       FB(9)=U1(3)+U1(5) 
       FB(10)=U1(4)+U1(6)
       FB(11)=U1(9)-U1(11)
       FB(12)=U1(10)-U1(12)
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
        FI(1)=reduce_angular_momentum(U,U-UOLD,1,2)
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
        REAL(8)::getp
        REAL(8),DIMENSION(totaldim)::ic
        INTEGER::I
      
        ic=(/(getp("BV0",I,U) ,I=1,totaldim)/)
        PAR(4)=nbp_Hamiltonian(ic,(/PAR(1),PAR(2),PAR(3)/),PAR(6))
        PAR(5)=getp("BIF",1,U)
        RETURN 
      END 
!---------------------------------------------------------------------- 
