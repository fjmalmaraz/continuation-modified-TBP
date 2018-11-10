MODULE nbp
  IMPLICIT NONE 
  PRIVATE::nbpaux 
  PUBLIC::nbpeq,grad,print_initial_solution,nbp_Hamiltonian
  PUBLIC::reduce_linear_momentum,reduce_angular_momentum,reduce_phase  
  INCLUDE 'nbp_constants_potential.h'
  INTEGER,PARAMETER::nfg=nconf*nbody
  INTEGER,PARAMETER::totaldim=2*nfg
  CONTAINS 
    SUBROUTINE nbpaux(q1,q2,f1,f2,m1,m2,alpha)
          REAL(8), DIMENSION(nconf),INTENT(in)::q1,q2
          REAL(8), DIMENSION(nconf),INTENT(in out)::f1,f2
          REAL(8), INTENT(in)::m1,m2
          REAL(8), INTENT(in)::alpha
          REAL(8)::aux
          REAL(8), DIMENSION(nconf)::q12
          aux=0.
          q12=q1-q2
          aux=alpha*DOT_PRODUCT(q12,q12)**(-(alpha+2)*0.5D0)
          f1=f1-m2*aux*q12
          f2=f2+m1*aux*q12
          RETURN 
     END SUBROUTINE nbpaux
     !
     SUBROUTINE nbpeq(U,mass,F,alpha)
          REAL(8), TARGET,DIMENSION(totaldim),INTENT(in)::U
          REAL(8), TARGET,DIMENSION(totaldim),INTENT(out)::F
          REAL(8), DIMENSION(nbody),INTENT(in)::mass
          REAL(8), INTENT(in)::alpha
          REAL(8), POINTER,DIMENSION(:)::u1,u2,f1,f2
          INTEGER::i,j,ic,jc
          F((nfg+1):totaldim)=(/(0.D0,i=1,nfg)/)
          F(1:nfg)=U((nfg+1):totaldim)
          DO i=1,nbody
            ic=(i-1)*nconf
            u1=>U(ic+1:ic+nconf)
            f1=>F(nfg+ic+1:nfg+ic+nconf)
              DO j=i+1,nbody
               jc=(j-1)*nconf
               u2=>U((jc+1):(jc+nconf))
               f2=>F((nfg+jc+1):(nfg+jc+nconf))
               CALL nbpaux(u1,u2,f1,f2,mass(i),mass(j),alpha)
              END DO
          END DO
          RETURN
    END SUBROUTINE nbpeq
    !
    SUBROUTINE extract_position(U,I,q)
      REAL(8),DIMENSION(totaldim),INTENT(in)::U
      REAL(8),DIMENSION(nconf),INTENT(out)::q
      INTEGER::I
      q=U(((I-1)*nconf+1):I*nconf)
      RETURN
    END SUBROUTINE
    !
    SUBROUTINE extract_momentum(U,I,p)
      REAL(8),DIMENSION(totaldim),INTENT(in)::U
      REAL(8),DIMENSION(nconf),INTENT(out)::p
      INTEGER,INTENT(in)::I
      p=U((nfg+(I-1)*nconf+1):nfg+I*nconf)
      RETURN
    END SUBROUTINE
    !
    REAL(8) FUNCTION nbp_Hamiltonian(ic,mass,alpha)
         REAL(8),DIMENSION(totaldim),INTENT(in)::ic
         REAL(8),DIMENSION(nbody),INTENT(in)::mass
         REAL(8),INTENT(in)::alpha
         REAL(8),DIMENSION(nconf)::q1,q2,q12,p
         INTEGER::I,J
         REAL(8)::H
         DO I=1,nbody
           CALL extract_momentum(ic,I,p)
           H=H+0.5D0*DOT_PRODUCT(p,p)/mass(I)
         ENDDO
         DO I=1,(nbody-1)
           CALL extract_position(ic,I,q1)  !ic(((I-1)*nconf+1):I*nconf)
           DO J=(I+1),nbody 
             CALL extract_position(ic,J,q2) !q2=ic(((J-1)*nconf+1):J*nconf)
             q12=q2-q1
             H=H-mass(I)*mass(J)/(DOT_PRODUCT(q12,q12)**(-alpha*0.5D0))
           ENDDO
         ENDDO
         nbp_Hamiltonian=H
         RETURN 
    END FUNCTION nbp_Hamiltonian
    ! 
    REAL(8) FUNCTION SKEW_PRODUCT(X,Y)
          REAL(8), DIMENSION(totaldim),INTENT(in)::X,Y
          SKEW_PRODUCT=DOT_PRODUCT(X(1:nfg),Y(nfg+1:2*nfg))
          SKEW_PRODUCT=SKEW_PRODUCT-DOT_PRODUCT(X(nfg+1:2*nfg),Y(1:nfg))
          RETURN
    END FUNCTION SKEW_PRODUCT
    !
    REAL(8) FUNCTION reduce_linear_momentum(DU,I)
           REAL(8),DIMENSION(totaldim),INTENT(in)::DU 
           INTEGER::I
           reduce_linear_momentum=SUM(DU(I:nfg:nconf)) 
           RETURN
    END FUNCTION reduce_linear_momentum
    !
    REAL(8) FUNCTION reduce_angular_momentum(U,DU,I,J)
           REAL(8),DIMENSION(totaldim),INTENT(in)::U,DU 
           INTEGER,INTENT(in)::I,J
           REAL(8),DIMENSION(totaldim)::GR
           CALL grad(U,I,J,GR)
           reduce_angular_momentum=SKEW_PRODUCT(GR,DU) 
           RETURN
    END FUNCTION reduce_angular_momentum
    !
    REAL(8) FUNCTION reduce_phase(U,DU,mass,alpha)
           REAL(8),DIMENSION(totaldim),INTENT(in)::U,DU
           REAL(8),DIMENSION(nbody),INTENT(in)::mass
           REAL(8),INTENT(in)::alpha
           REAL(8),DIMENSION(totaldim)::F
           CALL nbpeq(U,mass,F,alpha)
           reduce_phase=DOT_PRODUCT(F,DU) 
           RETURN
    END FUNCTION reduce_phase
    !
    SUBROUTINE grad(U,I1,I2,GR)
          REAL(8), DIMENSION(2*nfg),INTENT(in)::U
          INTEGER, INTENT(in) ::I1,I2
          REAL(8), DIMENSION(2*nfg),INTENT(out)::GR
          INTEGER::I
          GR=(/(0,I=1,2*nfg)/)
          GR(I1:nfg:nconf)=U((nfg+I2):2*nfg:nconf)
          GR(I2:nfg:nconf)=-U((nfg+I1):2*nfg:nconf)
          GR((nfg+I1):(2*nfg):nconf)=-U(I2:nfg:nconf)
          GR((nfg+I2):(2*nfg):nconf)=U(I1:nfg:nconf)
          RETURN
    END SUBROUTINE grad
    !
    SUBROUTINE perturbate_Hamiltonian(F,U,mass,alpha,thisparameter)
      REAL(8),INTENT(in out),DIMENSION(totaldim):: F
      REAL(8),INTENT(in),DIMENSION(totaldim):: U
      REAL(8),INTENT(in),DIMENSION(nbody):: mass
      REAL(8),INTENT(in)::alpha
      REAL(8),INTENT(in)::thisparameter
      REAL(8),DIMENSION(totaldim)::FOLD
      !CALL nbpeq(U,(/PAR(1),1.D0,1.D0/),F)
      CALL nbpeq(U,mass,F,alpha)
      FOLD=F
      F(1:nfg)=F(1:nfg)-thisparameter*FOLD(nfg+1:totaldim)
      F(nfg+1:totaldim)=F(nfg+1:totaldim)+thisparameter*FOLD(1:nfg)
      RETURN
    END SUBROUTINE perturbate_Hamiltonian
    !
    SUBROUTINE perturbate_linear_momentum(F,component,theparameter)
      INTEGER,INTENT(in):: component
      REAL(8),INTENT(in out),DIMENSION(totaldim):: F
      REAL(8) ,INTENT(in):: theparameter
      INTEGER::inicio,I
      inicio=nfg+component
      F(inicio:totaldim:nconf)=(/(F(I)+theparameter,I=inicio,totaldim,nconf)/)
      RETURN
    END SUBROUTINE perturbate_linear_momentum
    !
    SUBROUTINE perturbate_angular_momentum(F,U,compI,compJ,theparameter)
      INTEGER,INTENT(in):: compI,compJ
      REAL(8),INTENT(in out),DIMENSION(totaldim):: F
      REAL(8),INTENT(in),DIMENSION(totaldim):: U
      REAL(8) ,INTENT(in):: theparameter
      INTEGER::inicio,I
      REAL(8),DIMENSION(totaldim)::GR
      CALL grad(U,compI,compJ,GR)
      F=F+theparameter*GR
      RETURN
    END SUBROUTINE perturbate_angular_momentum
    !
    SUBROUTINE print_initial_solution(x0,mass,alpha,PER,DIV)
      REAL(8),DIMENSION(totaldim),INTENT(in)::x0
      REAL(8),DIMENSION(nbody),INTENT(in)::mass
      REAL(8),INTENT(in)::alpha
      REAL(8),INTENT(in)::PER
      INTEGER,INTENT(in)::DIV
      REAL(8),DIMENSION(totaldim):: y,F,yaux,ynew
      REAL(8)::h
      INTEGER::I
      h=PER/DIV
      y=x0
      OPEN(UNIT=10,FILE="nbp.dat",STATUS="unknown")
      write(10,*) 0.D0,y
      DO I=1,DIV
        CALL nbpeq(y,mass,F,alpha)
        yaux=y+h*F*0.5
        ynew=F/6.0D0;
        CALL nbpeq(yaux,mass,F,alpha)
        yaux=y+h*F*0.5
        ynew=ynew+F/3.0D0
        CALL nbpeq(yaux,mass,F,alpha)
        yaux=y+h*F
        ynew=ynew+F/3.D0
        CALL nbpeq(yaux,mass,F,alpha)
        ynew=(ynew+F/6.D0)*h
        y=y+ynew
        IF (MOD(I,10)==0) THEN
          WRITE(10,*) I*h,y
        END IF
      END DO
      CLOSE(10)
      RETURN
    END SUBROUTINE print_initial_solution
END MODULE nbp 

