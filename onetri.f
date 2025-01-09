! This file is part of std2.
!
! Copyright (C) 2013-2025 Stefan Grimme and Marc de Wergifosse
!
! std2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! std2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with std2.  If not, see <https://www.gnu.org/licenses/>.
!
!! ------------------------------------------------------------------------
!
      SUBROUTINE ONETRI(ITY,S,S1,S2,ARRAY,N,IVAL)
      IMPLICIT REAL*8 (A-H,O-Z)
C     ******DESIGNED FOR ABELIAN GROUPS ONLY******
C
C     CALLING SEQUENCE:
C     ITY       =1 SYM AO
C               =-1 ANTI SYM AO
C     S         INPUT PROPERTY MATRIX OVER AO'S
C     S1        TRANSFORMED INTEGRALS ON OUTPUT
C     S2        SCRATCH ARRAYS LARGE ENOUGH TO HOLD A SQUARE MATRIX
C     ARRAY     MO MATRIX OVER SO'S
C     N         LINEAR DIMENSION OF ARRAYS
C     BERND HESS, UNIVERSITY OF BONN, JANUARY 1991
      DIMENSION S(*),S1(*),S2(*),ARRAY(N,IVAL)

C
C     DETERMINE IF WE HAVE AN ANTISYMMETRIC INTEGRAL
      IF (ITY.EQ.-1) GOTO 99
C
C     BLOW UP SYMMETRIC MATRIX S
      CALL BLOWSY(ITY,S,S1,N)

C
C     TRANSFORMATION OF S
      CALL DSYMM('L','L',N,IVAL,1.D0,S1,N,ARRAY,N,0.D0,S2,N)
      CALL DGEMM('T','N',IVAL,IVAL,N,1.D0,ARRAY,N,S2,N,0.D0,S1,IVAL)
      RETURN
99    CONTINUE
C
C     BLOW UP ANTI-SYMMETRIC MATRIX S
      CALL BLOWSY(ITY,S,S1,N)
C
C     TRANSFORMATION OF S
      CALL DGEMM('N','N',N,IVAL,N,1.D0,S1,N,ARRAY,N,0.D0,S2,N)
      CALL DGEMM('T','N',IVAL,IVAL,N,1.D0,ARRAY,N,S2,N,0.D0,S1,IVAL)
      END
      SUBROUTINE BLOWSY(ITY,A,B,N)
C
C     BLOW UP SYMMETRIC OR ANTISYMMETRIC MATRIX TO FULL SIZE
      REAL*8 A(*),B(N,N)
C
C     DETERMINE IF WE HAVE AN ANTISYMMETRIC INTEGRAL

      IF (ITY.EQ.-1) GOTO 99
      IJ=0
      DO 1 I=1,N
      DO 2 J=1,I-1
      IJ=IJ+1
      B(J,I)=A(IJ)
2     B(I,J)=A(IJ)
      IJ=IJ+1
      B(I,I)=A(IJ)
1     CONTINUE
      RETURN
99    IJ=0
      DO 11 I=1,N
      DO 12 J=1,I-1
      IJ=IJ+1
      B(J,I)=-A(IJ)
12    B(I,J)=A(IJ)
      IJ=IJ+1
      B(I,I)=0.D0
11    CONTINUE
      RETURN
      END
