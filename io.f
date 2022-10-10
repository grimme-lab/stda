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
***********************************************************************

      subroutine mwrite(n,iwo,v,irec)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension v(n)
      write(iwo,rec=irec) v
      return
      end


      subroutine mread(n,iwo,v,irec)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension v(n)
      read(iwo,rec=irec) v
      return
      end
      
