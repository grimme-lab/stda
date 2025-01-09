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
      SUBROUTINE pckao(NPR,NAO,A,B)
      use stdacommon
      IMPLICIT REAL*8(A-H,O-Z)
      integer*8 ij,kl,k,iaa,iii

      dimension a(*),b(*)

      if(nao.eq.0) then
      k=0
      do i=1,npr
         do j=1,i
            k=k+1
            b(k)=a(k)
         enddo
      enddo
      return
      endif

      ij=nao
      ij=ij*(ij+1)/2
      b(1:ij)=0.0d0

      kl=0
      do i=1,npr
         iai=ipao(i)
         c1=cxip(i)
         do j=1,i-1
            kl=kl+1
            c2=cxip(j)
            iaj=ipao(j)
            iaa=max(iaj,iai)
            iii=min(iaj,iai)
            ij=iii+iaa*(iaa-1)/2
            b(ij)=b(ij)+a(kl)*c1*c2*2.0d0
         enddo
         kl=kl+1
         ij=iai
         ij=ij+ij*(ij-1)/2
         b(ij)=b(ij)+a(kl)*c1*c1
      enddo

      ij=0
      do i=1,nao
         do j=1,i-1
            ij=ij+1
            b(ij)=b(ij)*0.5
         enddo
         ij=ij+1
      enddo

      return
      end

      SUBROUTINE pckao3(NPR,NAO,A1,A2,A3,B1,B2,B3)
      use stdacommon
      IMPLICIT REAL*8(A-H,O-Z)
      integer*8 ij,kl,k

      dimension a1(*),b1(*)
      dimension a2(*),b2(*)
      dimension a3(*),b3(*)

      if(nao.eq.0) then
      k=0
      do i=1,npr
         do j=1,i
            k=k+1
            b1(k)=a1(k)
            b2(k)=a2(k)
            b3(k)=a3(k)
         enddo
      enddo
      return
      endif

      ij=nao
      ij=ij*(ij+1)/2
      b1(1:ij)=0.0d0
      b2(1:ij)=0.0d0
      b3(1:ij)=0.0d0

      kl=0
      do i=1,npr
         iai=ipao(i)
         c1=cxip(i)
         do j=1,i-1
            kl=kl+1
            c2=cxip(j)
            iaj=ipao(j)
            iaa=max(iaj,iai)
            iii=min(iaj,iai)
            ij=iii+iaa*(iaa-1)/2
            ccf=c1*c2*2.0d0
            b1(ij)=b1(ij)+a1(kl)*ccf
            b2(ij)=b2(ij)+a2(kl)*ccf
            b3(ij)=b3(ij)+a3(kl)*ccf
         enddo
         kl=kl+1
         ij=iai
         ij=ij+ij*(ij-1)/2
         b1(ij)=b1(ij)+a1(kl)*c1*c1
         b2(ij)=b2(ij)+a2(kl)*c1*c1
         b3(ij)=b3(ij)+a3(kl)*c1*c1
      enddo

      ij=0
      do i=1,nao
         do j=1,i-1
            ij=ij+1
            b1(ij)=b1(ij)*0.50d0
            b2(ij)=b2(ij)*0.50d0
            b3(ij)=b3(ij)*0.50d0
         enddo
         ij=ij+1
      enddo

      return
      end
