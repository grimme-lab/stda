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
ccccccccccccccccccccccccccccccccccc
! optional: Normalization of AOs  c
ccccccccccccccccccccccccccccccccccc
! This normalizes the contraction coefficients for each contraction
      subroutine normalize(cartbas,nprims,ipao,ipty,exip,cxip)
      implicit none
      logical, intent ( in ) :: cartbas
      integer, intent ( in ) :: ipao(nprims),ipty(nprims)
      integer, intent ( in ) :: nprims
      real*8, intent( in ) :: exip(nprims)
      real*8, intent( inout ) :: cxip(nprims)
      integer i,j,k,l,iprimao,iprimtyp,jprimao,jprimtyp
      integer ifac,lang,lx,my,nz
      real*8 fnorm,summe,dzaehl,dnenn,expon,xlinf

      open(unit=11,file='fnorm')
      write(*,'(A)',advance='no') 'normalizing...'
      j=1
      do i=1,nprims+1
       if(i.le.nprims) then
         iprimao=ipao(i)
         iprimtyp=ipty(i)
       endif
       if(i.gt.1)then
        if(i.eq.nprims+1) iprimao=0
        if(iprimao.ne.jprimao) then
         call deflmna(jprimtyp,lx,my,nz,lang,xlinf)
! we override xlinf here, if a Cartesian basis is present (i.e. selv overlap is 1.0)
         if(cartbas)  xlinf=1.0d0
         fnorm=3.14159265358979323846**1.50d0
         fnorm=fnorm/dble(2**lang)
         ifac=0
         call dblfac(lx,ifac)
         fnorm=fnorm*dble(ifac)
         call dblfac(my,ifac)
         fnorm=fnorm*dble(ifac)
         call dblfac(nz,ifac)
         fnorm=fnorm*dble(ifac)
         fnorm=dsqrt(fnorm)
         summe=0.0d0
         do k=j,i-1
          do l=j,i-1
           dzaehl=cxip(k)*cxip(l)
           dnenn=exip(k)+exip(l)
           expon=dble(lang)+1.50d0
           dnenn=dnenn**expon
           summe=summe+dzaehl/dnenn
          enddo
         enddo
         summe=dsqrt(summe)
         fnorm=fnorm*summe
         fnorm=xlinf/fnorm
         do k=j,i-1
!          write(*,*) k,ipty(k)
!          write(*,*) exip(k),cxip(k),fnorm*cxip(k)
          cxip(k)=fnorm*cxip(k)
          write(11,*)fnorm
!          write(*,*)k,jprimao,jprimtyp,cxip(k),cxip(k)**2
         enddo
         j=i
        endif
       endif
       jprimao=iprimao
       jprimtyp=iprimtyp
      enddo
      close(11)
!!!!!!!!!!!!!!!!!!!!!!!
!! check normalization
!      j=1
!      do i=1,nprims+1
!       iprimao=ipao(i)
!       iprimtyp=ipty(i)
!       if(i.gt.1)then
!        if(i.eq.nprims+1) iprimao=0
!        if(iprimao.ne.jprimao) then
!         call deflmna(jprimtyp,lx,my,nz,lang,xlinf)
!         fnorm=3.14159265358979323846**1.50d0
!         fnorm=fnorm/dble(2**lang)
!         ifac=0
!         call dblfac(lx,ifac)
!         fnorm=fnorm*ifac
!         call dblfac(my,ifac)
!         fnorm=fnorm*ifac
!         call dblfac(nz,ifac)
!         fnorm=fnorm*ifac
!         summe=0.0d0
!         do k=j,i-1
!          do l=j,i-1
!           dzaehl=cxip(k)*cxip(l)
!           dnenn=exip(k)+exip(l)
!           expon=dble(lang)+1.50d0
!           dnenn=dnenn**expon
!           summe=summe+dzaehl/dnenn
!          enddo
!         enddo
!         summe=fnorm*summe
!         write(*,*) 'self ovlp:', summe
!         j=i
!        endif
!       endif
!       jprimao=iprimao
!       jprimtyp=iprimtyp
!      enddo
!!!
      end subroutine normalize


      subroutine deflmna(iprtyp,lx,my,nz,lang,xlinf)
      implicit none
      integer lx,my,nz,lang,iprtyp
      real*8 xlinf
      xlinf=1.0d0
      lx=0
      my=0
      nz=0
      lang=0
c=======================================================================
c cartesian gaussian functions (6d,10f...)
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c
c assign for each ipty, the angular momentum (L), and the expinents l,m,n of x,y,z
c linf is the factor to multiply normalized functions that are linearly dependent (e.g. dx**2,dy**2,dz**2)
c=======================================================================
      select case(iprtyp)
       case(1)
        lx=0
        my=0
        nz=0
        lang=0
        xlinf=1.0d0
       case(2)
        lx=1
        my=0
        nz=0
        lang=1
        xlinf=1.0d0
       case(3)
        lx=0
        my=1
        nz=0
        lang=1
        xlinf=1.0d0
       case(4)
        lx=0
        my=0
        nz=1
        lang=1
        xlinf=1.0d0
       case(5)
        lx=2
        my=0
        nz=0
        lang=2
        xlinf=dsqrt(3.0d0)
       case(6)
        lx=0
        my=2
        nz=0
        lang=2
        xlinf=dsqrt(3.0d0)
       case(7)
        lx=0
        my=0
        nz=2
        lang=2
        xlinf=dsqrt(3.0d0)
       case(8)
        lx=1
        my=1
        nz=0
        lang=2
        xlinf=dsqrt(3.0d0)
c        xlinf=1.0d0
       case(9)
        lx=1
        my=0
        nz=1
        lang=2
        xlinf=dsqrt(3.0d0)
c        xlinf=1.0d0
       case(10)
        lx=0
        my=1
        nz=1
        lang=2
        xlinf=dsqrt(3.0d0)
c        xlinf=1.0d0
       case(11)
        lx=3
        my=0
        nz=0
        lang=3
        xlinf=dsqrt(15.0d0)
       case(12)
        lx=0
        my=3
        nz=0
        lang=3
        xlinf=dsqrt(15.0d0)
       case(13)
        lx=0
        my=0
        nz=3
        lang=3
        xlinf=dsqrt(15.0d0)
       case(14)
        lx=2
        my=1
        nz=0
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(15)
        lx=2
        my=0
        nz=1
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(16)
        lx=1
        my=2
        nz=0
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(17)
        lx=0
        my=2
        nz=1
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(18)
        lx=1
        my=0
        nz=2
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(19)
        lx=0
        my=1
        nz=2
        lang=3
!        xlinf=dsqrt(3.0d0)
        xlinf=dsqrt(15.0d0)
       case(20)
        lx=1
        my=1
        nz=1
        lang=3
!        xlinf=1.0d0
        xlinf=dsqrt(15.0d0)
       case default
        write(*,*)'unrecognized cartesian function/ang. momentum'
        write(0,*)'unrecognized cartesian function/ang. momentum'
        stop 'normalization impossible! Exiting...'
       end select
! for checking - is the factor important??
!       xlinf=1.0d0
      return
      end

      subroutine dblfac(iin,iout)
      implicit none
      integer i,iin,iout,jdem,jnum,jfac
cccccccccccccccccccccccccccccccccccccccccccccc
c
c this calculates the double faculty (2n-1)!!
c i.e. expression in brackets is always odd
c
cccccccccccccccccccccccccccccccccccccccccccccc
      jdem=1
      do i=1,iin
       jdem=jdem*i
      enddo
      jnum=jdem
      do i=iin+1,2*iin
       jdem=jdem*i
      enddo
      jfac=2**iin
      jnum=jnum*jfac
      iout=jdem/jnum
      return
      end
