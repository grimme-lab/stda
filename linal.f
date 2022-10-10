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
      subroutine blow(nbf,F,X)
      implicit none
      real*8  F(nbf*(nbf+1)/2),X(nbf,nbf)
      integer nbf,i,j,k

c blow it up
      k=0
      do i=1,nbf
         do j=1,i
            k=k+1
            X(i,j)=F(k)
            X(j,i)=F(k)
         enddo
      enddo

      end

