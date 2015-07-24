      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
C
C     Single value decomposition.
C     Adaptation from numerical recipes, 2nd ed..
C          K. H. FHI, September 2003
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'param.h'
      PARAMETER (TINY=1.D-14)
C
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NAUXFN)
C
      ZERO=0.D0
      ONE=1.D0
      ONET=ONE+TINY
      ITMAX=60
C
      G=ZERO
      SCALE=ZERO
      ANORM=ZERO
      DO I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M) THEN
           DO K=I,M
             SCALE=SCALE+DABS(A(K,I))
           ENDDO
           IF (SCALE.NE.ZERO) THEN
              SCALINV=ONE/SCALE
              DO K=I,M
                A(K,I)=A(K,I)*SCALINV
                S=S+A(K,I)**2
              ENDDO
              F=A(I,I)
              G=-SIGN(DSQRT(S),F)
              H=F*G-S
              A(I,I)=F-G
              NUM =M-I+1
              HINV=ONE/H
              DO J=L,N
                S=ZERO
                IF (NUM.GT.0) THEN
                   S=DDOT(NUM,A(I,I),1,A(I,J),1)
                   F=S*HINV
                   CALL DAXPY(NUM,F,A(I,I),1,A(I,J),1)
                ENDIF
              ENDDO
              DO  K=I,M
                A(K,I)=SCALE*A(K,I)
              ENDDO
           ENDIF
        ENDIF
C
        W(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M.AND.I.NE.N) THEN
           DO K=L,N
             SCALE=SCALE+DABS(A(I,K))
           ENDDO
           IF (SCALE.NE.ZERO) THEN
              SCALINV=ONE/SCALE
              DO K=L,N
                A(I,K)=A(I,K)*SCALINV
                S=S+A(I,K)**2
              ENDDO
              F=A(I,L)
              G=-SIGN(DSQRT(S),F)
              H=F*G-S
              A(I,L)=F-G
              HINV=ONE/H
              DO K=L,N
                RV1(K)=A(I,K)*HINV
              ENDDO
              NUM=N-L+1
              DO J=L,M
                S=ZERO
                IF (NUM.GT.0) THEN
                   S=DDOT(NUM,A(J,L),MP,A(I,L),MP)
                   CALL DAXPY(NUM,S,RV1(L),1,A(J,L),MP)
                ENDIF
              ENDDO
              DO K=L,N
                A(I,K)=A(I,K)*SCALE
              ENDDO
           ENDIF
        ENDIF
        ANORM=DMAX1(ANORM,(DABS(W(I))+DABS(RV1(I))))
        ANORMM=ANORM*ONET
      ENDDO
C
      DO I=N,1,-1
        IF (I.LT.N) THEN
           IF (G.NE.ZERO) THEN
              DO J=L,N
                V(J,I)=(A(I,J)/A(I,L))/G
              ENDDO
              NUM=N-L+1
              DO J=L,N
                S=ZERO
                IF (NUM.GT.0) THEN
                   S=DDOT(NUM,A(I,L),MP,V(L,J),1)
                   CALL DAXPY(NUM,S,V(L,I),1,V(L,J),1)
                ENDIF
              ENDDO
           ENDIF
           DO J=L,N
             V(I,J)=ZERO
             V(J,I)=ZERO
           ENDDO
        ENDIF
        V(I,I)=ONE
        G=RV1(I)
        L=I
      ENDDO
C
      DO I=MIN0(M,N),1,-1
         L=I+1
         G=W(I)
         DO J=L,N
            A(I,J)=ZERO
         ENDDO
         IF (G.NE.ZERO) THEN
            G=ONE/G
            NUM1=M-L+1
            NUM2=M-I+1
            DO J=L,N
              S=ZERO
              IF (NUM1.GT.0) S=DDOT(NUM1,A(L,I),1,A(L,J),1)
              F=S/A(I,I)*G
              IF (NUM2.GT.0) CALL DAXPY(NUM2,F,A(I,I),1,A(I,J),1)
            ENDDO
            DO J=I,M
              A(J,I)=A(J,I)*G
            ENDDO
         ELSE
            DO J=I,M
              A(J,I)=ZERO
            ENDDO
         ENDIF
         A(I,I)=A(I,I)+ONE
      ENDDO
C
      DO K=N,1,-1
        DO ITS=1,ITMAX
          DO L=K,1,-1
            NM=L-1
CKH            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GOTO 2
CKH            IF ((DABS(W(NM)) +ANORM).EQ.ANORM)  GOTO 1
            IF ((DABS(RV1(L))+ANORM).LE.ANORMM)  GOTO 2
            IF ((DABS(W(NM)) +ANORM).LE.ANORMM)  GOTO 1
          ENDDO
C
    1     C=ZERO
          S=ONE
          DO I=L,K
            F=S*RV1(I)
            RV1(I)=C*RV1(I)
CKH            IF (DABS(F)+ANORM.NE.ANORM) GOTO 2
            IF ((DABS(F)+ANORM).GT.ANORMM) GOTO 2
            G=W(I)
            H=PYTHAG(F,G)
            W(I)=H
            H=ONE/H
            C=G*H
            S=-F*H
            DO J=1,M
              Y=A(J,NM)
              Z=A(J,I)
              A(J,NM)=Y*C+Z*S
              A(J,I)=-Y*S+Z*C
            ENDDO
          ENDDO
C
    2     Z=W(K)
          IF (L.EQ.K) THEN
             IF (Z.LT.ZERO) THEN
                W(K)=-Z
                DO J=1,N
                  V(J,K)=-V(J,K)
                ENDDO
             ENDIF
             GOTO 3
          ENDIF
          IF (ITS.EQ.ITMAX) THEN
             WRITE(6,'(//A/8X,A,I4,A)')
     +           ' WARNING: failure in singular value decomposition,',
     +           'no convergence after',ITMAX,
     +           ' iterations (subroutine svdcmp)'
C             STOP
          ENDIF
C
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0D0*H*Y)
          G=PYTHAG(F,ONE)
          F =((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
C
          C=ONE
          S=ONE
          DO J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=PYTHAG(F,H)
            RV1(J)=Z
            ZINV =ONE/Z
            C=F*ZINV
            S=H*ZINV
            F= X*C+G*S
            G=-X*S+G*C
            H=Y*S
            Y=Y*C
            DO II=1,N
              X=V(II,J)
              Z=V(II,I)
              V(II,J)= X*C+Z*S
              V(II,I)=-X*S+Z*C
            ENDDO
C
            Z=PYTHAG(F,H)
            W(J)=Z
            IF (Z.NE.ZERO) THEN
               Z=ONE/Z
               C=F*Z
               S=H*Z
            ENDIF
            F= C*G+S*Y
            X=-S*G+C*Y
            DO II=1,M
              Y=A(II,J)
              Z=A(II,I)
              A(II,J)= Y*C+Z*S
              A(II,I)=-Y*S+Z*C
            ENDDO
          ENDDO
          RV1(L)=ZERO
          RV1(K)=F
          W(K)=X
        ENDDO
    3   CONTINUE
      ENDDO
      RETURN
      END

      REAL*8 FUNCTION PYTHAG(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      ABSA=DABS(A)
      ABSB=DABS(B)
      IF (ABSA.GT.ABSB) THEN
         PYTHAG=ABSA*DSQRT(1.D0+(ABSB/ABSA)**2)
      ELSE
         IF(ABSB.EQ.0.D0) THEN
            PYTHAG=0.D0
         ELSE
            PYTHAG=ABSB*DSQRT(1.D0+(ABSA/ABSB)**2)
         ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE SVDCMPO(A,M,N,MP,NP,W,V)
C
C     Single value decomposition. Adaptation from numerical recipes.
C          K. H. FHI, August 2003
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'param.h'
C
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NAUXFN)
C
      ZERO=0.D0
      ONE=1.D0
      ITMAX=MIN0(M,N)
C
      G=ZERO
      SCALE=ZERO
      ANORM=ZERO
      DO I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M) THEN
           DO K=I,M
             SCALE=SCALE+DABS(A(K,I))
           ENDDO
           IF (SCALE.NE.ZERO) THEN
              SCALINV=ONE/SCALE
              DO K=I,M
                A(K,I)=A(K,I)*SCALINV
                S=S+A(K,I)**2
              ENDDO
              F=A(I,I)
              G=-SIGN(DSQRT(S),F)
              H=F*G-S
              A(I,I)=F-G
              IF (I.NE.N) THEN
                 NUM =M-I+1
                 HINV=ONE/H
                 DO J=L,N
                   S=DDOT(NUM,A(I,I),1,A(I,J),1)
                   F=S*HINV
                   CALL DAXPY(NUM,F,A(I,I),1,A(I,J),1)
                 ENDDO
              ENDIF
              DO  K=I,M
                A(K,I)=SCALE*A(K,I)
              ENDDO
           ENDIF
        ENDIF
C
        W(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M.AND.I.NE.N) THEN
           DO K=L,N
             SCALE=SCALE+DABS(A(I,K))
           ENDDO
           IF (SCALE.NE.ZERO) THEN
              SCALINV=ONE/SCALE
              DO K=L,N
                A(I,K)=A(I,K)*SCALINV
                S=S+A(I,K)**2
              ENDDO
              F=A(I,L)
              G=-SIGN(DSQRT(S),F)
              H=F*G-S
              A(I,L)=F-G
              HINV=ONE/H
              DO K=L,N
                RV1(K)=A(I,K)*HINV
              ENDDO
              IF (I.NE.M) THEN
                 NUM=N-L+1
                 DO J=L,M
                   S=ZERO
                   IF (NUM.GT.0) S=DDOT(NUM,A(J,L),MP,A(I,L),MP)
                   CALL DAXPY(NUM,S,RV1(L),1,A(J,L),MP)
                 ENDDO
              ENDIF
              DO K=L,N
                A(I,K)=A(I,K)*SCALE
              ENDDO
           ENDIF
        ENDIF
        ANORM=DMAX1(ANORM,(DABS(W(I))+DABS(RV1(I))))
      ENDDO
C
      DO I=N,1,-1
        IF (I.LT.N) THEN
           IF (G.NE.ZERO) THEN
              GINV=ONE/(A(I,L)*G)
              DO J=L,N
                V(J,I)=A(I,J)*GINV
              ENDDO
              NUM=N-L+1
              DO J=L,N
                S=ZERO
                IF (NUM.GT.0) S=DDOT(NUM,A(I,L),MP,V(L,J),1)
                CALL DAXPY(NUM,S,V(L,I),1,V(L,J),1)
              ENDDO
           ENDIF
           DO J=L,N
             V(I,J)=ZERO
             V(J,I)=ZERO
           ENDDO
        ENDIF
        V(I,I)=ONE
        G=RV1(I)
        L=I
      ENDDO
C
      DO I=N,1,-1
         L=I+1
         G=W(I)
         IF (I.LT.N) THEN
            DO J=L,N
               A(I,J)=ZERO
            ENDDO
         ENDIF
         IF (G.NE.ZERO) THEN
            G=ONE/G
            IF (I.NE.N) THEN
               AINV=ONE/A(I,I)
               NUM1=M-L+1
               NUM2=M-I+1
               DO J=L,N
                 S=DDOT(NUM1,A(L,I),1,A(L,J),1)
                 F=S*AINV*G
                 CALL DAXPY(NUM2,F,A(I,I),1,A(I,J),1)
               ENDDO
            ENDIF
            DO J=I,M
              A(J,I)=A(J,I)*G
            ENDDO
         ELSE
            DO J=I,M
              A(J,I)=ZERO
            ENDDO
         ENDIF
         A(I,I)=A(I,I)+ONE
      ENDDO
C
      DO K=N,1,-1
        DO ITS=1,ITMAX
          DO L=K,1,-1
            NM=L-1
            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GOTO 2
            IF ((DABS(W(NM)) +ANORM).EQ.ANORM)  GOTO 1
          ENDDO
C
    1     C=ZERO
          S=ONE
          DO I=L,K
            F=S*RV1(I)
            IF (DABS(F)+ANORM.NE.ANORM) THEN
               G=W(I)
               H=DSQRT(F*F+G*G)
               W(I)=H
               H=ONE/H
               C=G*H
               S=-F*H
               DO J=1,M
                 Y=A(J,NM)
                 Z=A(J,I)
                 A(J,NM)=Y*C+Z*S
                 A(J,I)=-Y*S+Z*C
               ENDDO
            ENDIF
          ENDDO
C
    2     Z=W(K)
          IF (L.EQ.K) THEN
             IF (Z.LT.ZERO) THEN
                W(K)=-Z
                DO J=1,N
                  V(J,K)=-V(J,K)
                ENDDO
             ENDIF
             GOTO 3
          ENDIF
          IF (ITS.EQ.ITMAX) THEN
             WRITE(6,'(//A/8X,A,I4,A)')
     +           ' ERROR: failure in singular value decomposition,',
     +           'no convergence after',ITMAX,
     +           ' iterations (subroutine svdcmp)'
             STOP
          ENDIF
CKH ===============================
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0D0*H*Y)
          G=DSQRT(F*F+ONE)
          F =((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
C
          C=ONE
          S=ONE
          DO J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=DSQRT(F*F+H*H)
            RV1(J)=Z
            ZINV =ONE/Z
            C=F*ZINV
            S=H*ZINV
            F=X*C+G*S
            G=-X*S+G*C
            H=Y*S
            Y=Y*C
            DO II=1,N
              X=V(II,J)
              Z=V(II,I)
              V(II,J)=X*C+Z*S
              V(II,I)=-X*S+Z*C
            ENDDO
C
            Z=DSQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.ZERO) THEN
               Z=ONE/Z
               C=F*Z
               S=H*Z
            ENDIF
            F=C*G+S*Y
            X=-S*G+C*Y
            DO II=1,M
              Y=A(II,J)
              Z=A(II,I)
              A(II,J)=Y*C+Z*S
              A(II,I)=-Y*S+Z*C
            ENDDO
          ENDDO
          RV1(L)=ZERO
          RV1(K)=F
          W(K)=X
        ENDDO
    3   CONTINUE
      ENDDO
      RETURN
      END

      SUBROUTINE SVDCMPMC(A,M,N,MP,NP,W,V)
CKH
CKH   Single value decomposition. From numerical recipes.
CKH
      implicit real*8(a-h,o-z)
      parameter (itmax = 80)

      INCLUDE 'param.h'

      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NAUXFN)

      ZERO  = 0.D0
      ONE   = 1.D0

      G     = ZERO
      SCALE = ZERO
      ANORM = ZERO

      DO 25 I = 1,N
         L      = I + 1
         RV1(I) = SCALE*G
         G      = ZERO
         S      = ZERO
         SCALE  = ZERO

         IF (I .LE. M) THEN

            DO K = I,M
               SCALE = SCALE + DABS(A(K,I))
            ENDDO

            IF (SCALE .NE. ZERO) THEN

               SCALINV = ONE/SCALE
               DO K = I,M
                  A(K,I) = A(K,I)*SCALINV
                  S = S + A(K,I)*A(K,I)
               ENDDO
               F      = A(I,I)
               G      =-SIGN(DSQRT(S),F)
               H      = F*G - S
               A(I,I) = F-G
               IF (I .NE. N) THEN

                  NUM  = M - I + 1
                  HINV = ONE/H
                  DO J = L,N
                     S = DDOT(NUM,A(I,I),1,A(I,J),1)
                     F = S*HINV
                     CALL DAXPY(NUM,F,A(I,I),1,A(I,J),1)
                  ENDDO

               ENDIF
               DO  K = I,M
                  A(K,I) = SCALE*A(K,I)
               ENDDO

            ENDIF

         ENDIF
         W(I) = SCALE *G
         G    = ZERO
         S    = ZERO
         SCALE= ZERO
         IF ((I .LE. M) .AND. (I .NE. N)) THEN

            DO K = L,N
               SCALE = SCALE + DABS(A(I,K))
            ENDDO
            IF (SCALE .NE. ZERO) THEN

               SCALINV = ONE/SCALE
               DO K = L,N
                  A(I,K) = A(I,K)*SCALINV
                  S      = S + A(I,K)*A(I,K)
               ENDDO
               F = A(I,L)
               G =-SIGN(DSQRT(S),F)
               H = F*G - S
               A(I,L) = F - G
               HINV = ONE/H
               DO K = L,N
                  RV1(K) = A(I,K)*HINV
               ENDDO
               IF (I .NE. M) THEN

                  NUM = N - L + 1
                  DO J = L,M
                     S = ZERO
                     DO K = L,N
                        S = S + A(J,K)*A(I,K)
                     ENDDO
                     CALL DAXPY(NUM,S,RV1(L),1,A(J,L),MP)
                  ENDDO

               ENDIF
               DO K = L,N
                  A(I,K) = SCALE*A(I,K)
               ENDDO

            ENDIF

         ENDIF
         ANORM = MAX(ANORM,(DABS(W(I))+DABS(RV1(I))))
 25   CONTINUE

      DO I = N,1,-1
         IF (I .LT. N) THEN
            IF (G .NE. ZERO) THEN
               GINV = ONE/(A(I,L)*G)
               DO J = L,N
                  V(J,I) = A(I,J)*GINV
               ENDDO
               NUM = N - L + 1
               DO J = L,N
                  S = ZERO
                  DO K = L,N
                     S = S + A(I,K)*V(K,J)
                  ENDDO
                  CALL DAXPY(NUM,S,V(L,I),1,V(L,J),1)
               ENDDO
            ENDIF
            DO J = L,N
               V(I,J) = ZERO
               V(J,I) = ZERO
            ENDDO
         ENDIF
         V(I,I) = ONE
         G      = RV1(I)
         L      = I
      ENDDO

      DO I = N,1,-1
         L = I + 1
         G = W(I)
         IF (I .LT. N) THEN
            DO J = L,N
               A(I,J) = ZERO
            ENDDO
         ENDIF
         IF (G .NE. ZERO) THEN
            G = ONE/G
            IF (I .NE. N) THEN
               AINV = ONE/A(I,I)
               NUM1 = M - L + 1
               NUM2 = M - I + 1
               DO J = L,N
                  S = DDOT(NUM1,A(L,I),1,A(L,J),1)
                  F = S*AINV*G
                  CALL DAXPY(NUM2,F,A(I,I),1,A(I,J),1)
               ENDDO
            ENDIF
            DO J = I,M
               A(J,I) = A(J,I)*G
            ENDDO
         ELSE
            DO J = I,M
               A(J,I) = ZERO
            ENDDO
         ENDIF
         A(I,I) = A(I,I) + ONE
      ENDDO

c--- Checked till here

      DO 49 K = N,1,-1
         DO 48 ITS = 1,30
            DO L = K,1,-1
               NM = L - 1
               IF ((DABS(RV1(L))+ANORM) .EQ. ANORM)  GO TO 2
               IF ((DABS(W(NM)) +ANORM) .EQ. ANORM)  GO TO 1
            ENDDO
 1          C = ZERO
            S = ONE
            DO I = L,K
               F = S*RV1(I)
               IF ((DABS(F)+ANORM) .NE. ANORM) THEN
                  G    = W(I)
                  H    = DSQRT(F*F+G*G)
                  W(I) = H
                  H    = ONE/H
                  C    = (G*H)
                  S    =-(F*H)
                  DO J = 1,M
                     Y       = A(J,NM)
                     Z       = A(J,I)
                     A(J,NM) = (Y*C) + (Z*S)
                     A(J,I)  =-(Y*S)+(Z*C)
                  ENDDO
               ENDIF
            ENDDO
 2          Z = W(K)
            IF (L .EQ. K) THEN
               IF (Z .LT. ZERO) THEN
                  W(K) = -Z
                  DO J = 1,N
                     V(J,K) = -V(J,K)
                  ENDDO
               ENDIF
               GO TO 3
            ENDIF
CKH            IF (ITS.EQ.ITMAX) THEN
            IF (ITS.EQ.30) THEN
               WRITE(6,*) ' ***************************************'
               WRITE(6,*) ' FAILURE IN SINGULAR VALUE DECOMPOSITION'
               WRITE(6,*) ' (Subroutine svdcmp)'
               WRITE(6,*) ' No convergence in 30 iterations'
               WRITE(6,*) ' PROGRAM STOPS '
               WRITE(6,*) ' ***************************************'
               STOP
            ENDIF
            X  = W(L)
            NM = K - 1
            Y  = W(NM)
            G  = RV1(NM)
            H  = RV1(K)
            F  = ((Y-Z)*(Y+Z) + (G-H)*(G+H))/(2.0D0*H*Y)
            G  = DSQRT(F*F+ONE)
            F  =((X-Z)*(X+Z) + H*((Y/(F+SIGN(G,F)))-H))/X
            C  = ONE
            S  = ONE
            DO J = L,NM
               I     = J + 1
               G     = RV1(I)
               Y     = W(I)
               H     = S*G
               G     = C*G
               Z     = DSQRT(F*F+H*H)
               RV1(J)= Z
               ZINV  = ONE/Z
               C     = F*ZINV
               S     = H*ZINV
               F     = (X*C)+(G*S)
               G     =-(X*S)+(G*C)
               H     = Y*S
               Y     = Y*C
               DO II = 1,N
                  X = V(II,J)
                  Z = V(II,I)
                  V(II,J) = (X*C)+(Z*S)
                  V(II,I) =-(X*S)+(Z*C)
               ENDDO
               Z    = DSQRT(F*F+H*H)
               W(J) = Z
               IF (Z .NE. ZERO) THEN
                  Z = ONE/Z
                  C = F*Z
                  S = H*Z
               ENDIF
               F = (C*G)+(S*Y)
               X =-(S*G)+(C*Y)
               DO II = 1,M
                  Y = A(II,J)
                  Z = A(II,I)
                  A(II,J) = (Y*C)+(Z*S)
                  A(II,I) =-(Y*S)+(Z*C)
               ENDDO
            ENDDO
            RV1(L) = ZERO
            RV1(K) = F
            W(K)   = X
 48      CONTINUE
 3       CONTINUE
 49   CONTINUE

      RETURN
      END
