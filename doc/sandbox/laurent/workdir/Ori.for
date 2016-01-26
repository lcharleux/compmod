          SUBROUTINE  Ori(nblock, N, Nel, Q)
C
         IMPLICIT NONE
         INTEGER I, J, nblock, N, Nel, nset 
         REAL*8 S(100,143), Q(nblock,3,3), Ang(100,3), phi1, phi, phi2
C
         call Sets(S)           !Sets delements
         DO I=1,100
              DO J=1,143
                    IF (Nel==S(I,J)) THEN 
                        nset=I      
                        EXIT
                    ENDIF
              END DO
         END DO
C
         call Angles(Ang)
C
              phi1 = Ang(nset,1) 
              phi = Ang(nset,2) 
              phi2 = Ang(nset,3) 
C
         call Korientinit(phi1,phi,phi2,Q,nblock,N)
C
         RETURN
         END
