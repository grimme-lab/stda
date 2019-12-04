! This file is part of stda.
!
! Copyright (C) 2013-2019 Stefan Grimme
!
! stda is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! stda is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with stda.  If not, see <https://www.gnu.org/licenses/>.
!
      SUBROUTINE PRMAT(IUOUT,R,N,M,HEAD)
      CHARACTER*(*) HEAD
      real*8  R(*)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=8
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/,2X,A)
 1002 FORMAT(/,' ',4X,8(3X,I4,3X),/)
 1003 FORMAT(' ',I4,8F10.5)
      END

      SUBROUTINE PRMAT4(IUOUT,R,N,M,HEAD)
      CHARACTER*(*) HEAD
      real*4  R(*)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=8
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/,2X,A)
 1002 FORMAT(/,' ',4X,8(3X,I4,3X),/)
 1003 FORMAT(' ',I4,8F10.5)
      END
