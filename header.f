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
      subroutine header(aarg,iarg)
      IMPLICIT REAL*8(A-H,O-Z)
      character*(*) aarg

      write(*,110)
      if(iarg.ne.0) then
         write(*,120)aarg,iarg
      else
         write(*,121)aarg
      endif
      write(*,110)
110   format(70('='))
130   format(/)
120   format(20x,a,5x,i4)
121   format(20x,a)

      return
      end
