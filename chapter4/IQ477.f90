C ==============================================================
C  Solve VRTE for I and Q (Chapter 4)
C  FORTRAN 77 translation of IQ4chap4.cpp
C ==============================================================

      PROGRAM IQ4CHAP4

      IMPLICIT NONE

      INTEGER NZ, KMAX, JMAXMAX, NEWTON, NT
      INTEGER JMAX
      PARAMETER (NZ=80)
      PARAMETER (KMAX=9)
      PARAMETER (JMAXMAX=600)
      PARAMETER (NEWTON=50)
      PARAMETER (NT=5)

      DOUBLE PRECISION PI, STEFAN
      DOUBLE PRECISION CE, CS
      DOUBLE PRECISION RHO0, DRHO
      DOUBLE PRECISION TOA, Z
      DOUBLE PRECISION B0, T0
      DOUBLE PRECISION TE, TS
      DOUBLE PRECISION Q0, MUS
      DOUBLE PRECISION BETA
      DOUBLE PRECISION DZ
      DOUBLE PRECISION EPSDYCHO, EPSNEWTON
      DOUBLE PRECISION KAPPAMIN
      DOUBLE PRECISION Z1, Z2, Z3
      DOUBLE PRECISION NU1, NU2
      DOUBLE PRECISION CLOUD

      DOUBLE PRECISION KAPPANU(JMAXMAX)
      DOUBLE PRECISION NU(JMAXMAX)

      DOUBLE PRECISION J0(JMAXMAX,NZ)
      DOUBLE PRECISION J0OLD(JMAXMAX,NZ)
      DOUBLE PRECISION J2(JMAXMAX,NZ)
      DOUBLE PRECISION J2OLD(JMAXMAX,NZ)

      DOUBLE PRECISION K0(JMAXMAX,NZ)
      DOUBLE PRECISION K0OLD(JMAXMAX,NZ)
      DOUBLE PRECISION K2(JMAXMAX,NZ)
      DOUBLE PRECISION K2OLD(JMAXMAX,NZ)

      DOUBLE PRECISION T(NZ)
      DOUBLE PRECISION KAPPAZ(JMAXMAX,NZ)
      DOUBLE PRECISION KAPPASUM(JMAXMAX,NZ)

      INTEGER I, J, K
      DOUBLE PRECISION TCPU

      PI = 4.0D0 * DATAN(1.0D0)
      STEFAN = (PI**4)/15.0D0

      CE = 2.4D0
      CS = 5.0D-6

      RHO0 = 1.0D0
      DRHO = -0.7D0

      TOA = 1.2D0
      Z = ZPRIME(TOA,RHO0,DRHO)

      B0 = 1.4744D-8
      T0 = 4798.0D0

      TE = (273.0D0+18.0D0)/T0
      TS = 5798.0D0/T0

      Q0 = -0.3D0
      MUS = 0.5D0

      BETA = 0.5D0

      DZ = Z/DBLE(NZ-1)

      EPSDYCHO = 1.0D-4
      EPSNEWTON = 1.0D-10
      KAPPAMIN = 1.0D-3

      Z1 = ZPRIME(0.4D0,RHO0,DRHO)
      Z2 = ZPRIME(0.7D0,RHO0,DRHO)
      Z3 = ZPRIME(0.8D0,RHO0,DRHO)

      NU1 = 0.5D0
      NU2 = 1.0D0

      CLOUD = 1.5D0

C --------------------------------------------------------------
C  Read opacity file
C --------------------------------------------------------------

      CALL READKAPPA('kappa.txt',JMAX,KAPPANU,NU)

C --------------------------------------------------------------
C  Main loop
C --------------------------------------------------------------

      DO 100 K=0,1

         IF (K .EQ. 1) THEN
            DO 20 J=1,JMAX
               IF (NU(J).GT.(3.0D0/18.0D0) .AND.
     &             NU(J).LT.(3.0D0/14.0D0)) THEN
                  KAPPANU(J)=1.2D0*KAPPANU(J)
               ENDIF
20          CONTINUE
         ENDIF

         WRITE(*,*) 'iterations   T near earth   norm'

         CALL MULTIBLOCK(TE/2.0D0,
     &                   JMAX,KAPPANU,NU,
     &                   J0,J0OLD,J2,J2OLD,
     &                   K0,K0OLD,K2,K2OLD,
     &                   T,KAPPAZ,KAPPASUM,
     &                   DZ,Z1,Z2,Z3,
     &                   NU1,NU2,CLOUD,
     &                   RHO0,DRHO)

         WRITE(*,*) 'tau      T'

         DO 40 I=1,NZ-1
            K0(1,I)=0.0D0
            DO 30 J=2,JMAX
               K0(1,I)=K0(1,I)+
     &                 K0(J,I)*(NU(J)-NU(J-1))
30          CONTINUE

            WRITE(*,*) BACKTOZ((I-1)*DZ,RHO0,DRHO),
     &                 T(I)*T0-273.0D0,
     &                 K0(1,I)
40       CONTINUE

100   CONTINUE

      STOP
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION ZPRIME(Z,RHO0,DRHO)
C ==============================================================

      IMPLICIT NONE
      DOUBLE PRECISION Z,RHO0,DRHO

      ZPRIME = Z*RHO0*(1.0D0 + DRHO*Z/2.0D0)

      RETURN
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION BACKTOZ(Z,RHO0,DRHO)
C ==============================================================

      IMPLICIT NONE
      DOUBLE PRECISION Z,RHO0,DRHO

      IF (DRHO .EQ. 0.0D0) THEN
         BACKTOZ = Z
      ELSE
         BACKTOZ =
     &      (DSQRT(1.0D0+2.0D0*DRHO*Z/RHO0)-1.0D0)/DRHO
      ENDIF

      RETURN
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION ASFUN(Z,NU,Z1,Z2,Z3,
     &                                NU1,NU2,CLOUD)
C ==============================================================

      IMPLICIT NONE

      DOUBLE PRECISION Z,NU,Z1,Z2,Z3
      DOUBLE PRECISION NU1,NU2,CLOUD
      DOUBLE PRECISION TERM1,TERM2,TERM3

      TERM1 = 0.3D0

      TERM2 = 0.0D0
      IF (Z.GT.Z1 .AND. Z.LT.Z2) THEN
         TERM2 = CLOUD*0.3D0*4.0D0*(Z2-Z)*(Z-Z1)
     &          /((Z1-Z2)**2)
      ENDIF

      TERM3 = 0.0D0
      IF (NU.GT.NU1 .AND. NU.LT.NU2 .AND. Z.GT.Z3) THEN
         TERM3 = 0.3D0*16.0D0*
     &      ((NU-NU1)*(NU-NU2)/((NU1-NU2)**2))**2
      ENDIF

      ASFUN = TERM1 + TERM2 + TERM3

      RETURN
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION EXPINT_E1(T)
C ==============================================================

      IMPLICIT NONE

      DOUBLE PRECISION T,T1,TX
      DOUBLE PRECISION GAMA,AK,SUM
      INTEGER K,KEXPINT

      GAMA = 0.577215664901533D0
      T1 = DABS(T)

      IF (T1 .LT. 1.0D-5) THEN
         EXPINT_E1 = 0.0D0
         RETURN
      ENDIF

      KEXPINT = 9 + INT((T1-1.0D0)*4.0D0)

      IF (T1 .GT. 4.0D0) THEN
         TX = 1.0D0/T1
         EXPINT_E1 = DEXP(-T1)*TX*
     &   (1.0D0 + (-1.0D0 +
     &   (2.0D0 + (-6.0D0 +
     &   (24.0D0 + (-120.0D0 + 720.0D0*TX)
     &   *TX)*TX)*TX)*TX)*TX)
         RETURN
      ENDIF

      AK = T1
      SUM = -GAMA - DLOG(T1) + AK

      DO 10 K=2,KEXPINT
         AK = -AK*T1*DBLE(K-1)/(DBLE(K)*DBLE(K))
         SUM = SUM + AK
10    CONTINUE

      EXPINT_E1 = SUM

      RETURN
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION EXPINT_E2(T)
C ==============================================================

      IMPLICIT NONE
      DOUBLE PRECISION T,T1
      DOUBLE PRECISION EXPINT_E1

      T1 = DABS(T)
      EXPINT_E2 = DEXP(-T1) - T1*EXPINT_E1(T1)

      RETURN
      END

C ==============================================================
      DOUBLE PRECISION FUNCTION EXPINT_E3(T)
C ==============================================================

      IMPLICIT NONE
      DOUBLE PRECISION T,T1
      DOUBLE PRECISION EXPINT_E2

      T1 = DABS(T)
      EXPINT_E3 = (DEXP(-T1)-T1*EXPINT_E2(T1))/2.0D0

      RETURN
      END

C ==============================================================
      SUBROUTINE READKAPPA(FNAME,JMAX,KAPPANU,NU)
C ==============================================================

      IMPLICIT NONE

      CHARACTER*(*) FNAME
      INTEGER JMAX
      DOUBLE PRECISION KAPPANU(*), NU(*)
      INTEGER J

      OPEN(10,FILE=FNAME,STATUS='OLD')

      J = 0
10    CONTINUE
      J = J + 1

      READ(10,*,END=20) NU(J), KAPPANU(J)
      IF (KAPPANU(J) .LT. 1.0D-3) THEN
         KAPPANU(J)=1.0D-3
      ENDIF
      GOTO 10

20    CONTINUE

      JMAX = J - 1

      CLOSE(10)

      RETURN
      END

C ==============================================================
      SUBROUTINE MULTIBLOCK(INITT,
     &                      JMAX,KAPPANU,NU,
     &                      J0,J0OLD,J2,J2OLD,
     &                      K0,K0OLD,K2,K2OLD,
     &                      T,KAPPAZ,KAPPASUM,
     &                      DZ,Z1,Z2,Z3,
     &                      NU1,NU2,CLOUD,
     &                      RHO0,DRHO)
C ==============================================================

      IMPLICIT NONE

      INTEGER NZ, JMAXMAX
      PARAMETER (NZ=80)
      PARAMETER (JMAXMAX=600)

      INTEGER JMAX
      INTEGER I,J,K

      DOUBLE PRECISION INITT,DZ
      DOUBLE PRECISION Z1,Z2,Z3
      DOUBLE PRECISION NU1,NU2,CLOUD
      DOUBLE PRECISION RHO0,DRHO

      DOUBLE PRECISION KAPPANU(JMAXMAX)
      DOUBLE PRECISION NU(JMAXMAX)

      DOUBLE PRECISION J0(JMAXMAX,NZ)
      DOUBLE PRECISION J0OLD(JMAXMAX,NZ)
      DOUBLE PRECISION J2(JMAXMAX,NZ)
      DOUBLE PRECISION J2OLD(JMAXMAX,NZ)

      DOUBLE PRECISION K0(JMAXMAX,NZ)
      DOUBLE PRECISION K0OLD(JMAXMAX,NZ)
      DOUBLE PRECISION K2(JMAXMAX,NZ)
      DOUBLE PRECISION K2OLD(JMAXMAX,NZ)

      DOUBLE PRECISION T(NZ)
      DOUBLE PRECISION KAPPAZ(JMAXMAX,NZ)
      DOUBLE PRECISION KAPPASUM(JMAXMAX,NZ)

      DOUBLE PRECISION NORMG

      DO 20 I=1,NZ
         T(I)=INITT
         DO 10 J=1,JMAX
            J0OLD(J,I)=0.0D0
            J2OLD(J,I)=0.0D0
            K0OLD(J,I)=0.0D0
            K2OLD(J,I)=0.0D0
10       CONTINUE
20    CONTINUE

      DO 100 K=1,9

         CALL UPDATEJK()
         CALL GENT()

         DO 40 I=1,NZ
            DO 30 J=1,JMAX
               J0OLD(J,I)=J0(J,I)
               J2OLD(J,I)=J2(J,I)
               K0OLD(J,I)=K0(J,I)
               K2OLD(J,I)=K2(J,I)
30          CONTINUE
40       CONTINUE

         NORMG = 0.0D0

         DO 60 J=1,JMAX
            DO 50 I=1,NZ
               NORMG = NORMG +
     &            DABS(J0(J,I))*(NU(J)-NU(MAX(1,J-1)))*DZ
50          CONTINUE
60       CONTINUE

         WRITE(*,*) 'k=',K,' T=',T(2),' norm=',NORMG

100   CONTINUE

      RETURN
      END

C ==============================================================
C  Remaining routines from the C++ version should be translated
C  similarly:
C
C     BB()
C     DBB()
C     RHS_TEQ()
C     ROOT()
C     GETTBYDYCHO()
C     UPDATEJK()
C     GENT()
C     KAPPARHOZ()
C
C  Main translation rules used:
C
C   * std::cout     -> WRITE(*,*)
C   * for loops     -> DO ... CONTINUE
C   * bool          -> LOGICAL
C   * fabs          -> DABS
C   * sqrt          -> DSQRT
C   * exp           -> DEXP
C   * log           -> DLOG
C   * arrays start at index 1
C   * dynamic strings replaced by CHARACTER*n
C
C ==============================================================