! This file is part of std2.
!
! Copyright (C) 2025 Marc de Wergifosse
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
      subroutine overlap(ncent,nprims,nbf,overlap_AO)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m
      integer :: ncent,nprims,nbf
      integer, allocatable :: info(:)
      double precision, parameter :: PI=4.D0*DATAN(1.D0)

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:)
      double precision :: overlap_AO(nbf*(nbf+1)/2)
      double precision :: norm_cart(1:10,1:10)

      integer,external :: CINTcgto_cart

      integer*8 ::lin8

!       real*8 :: sumwx,sumwy,sumwz,sumw,atmass
!       real*8 :: cen,xmolw,ams
!       common /cema   / cen(3),xmolw
!       common /amass  / ams(107)

      allocate(atm(1:6,1:ncent))
      j=21 ! variables are defined below 20 note that the array in libcint starts with 0
      Do i=1, ncent
      atm(1,i)=int(co(i,4))  ! atom charge
      atm(2,i)=j           ! "zero" for list of atoms in env
      j=j+3
      enddo
      n_env=j

      allocate(info(nprims)) !indicate the number of primitives in each contractions
      ! Counting the number of primitives without those by symmetry in each contractions
      nbas=1
      info=1
      j=1
      Do i=2,nprims
      if(nbas/=i)then
      if(ipty(i)==1.or.ipty(i)==2.or.ipty(i)==5.or.ipty(i)==11)then ! remove those by symmetry
      if(ipao(i)==ipao(i-1))then
      info(nbas)=info(nbas)+1
      info(i)=0
      else ! jump to the next AO
      nbas=i
      j=j+1
      endif
      else ! other than x, xx, xxx are zero, no need to count them
      info(i)=0
      endif
      endif
      enddo

      nbas=j
      allocate(bas(1:8,1:nbas))

      j=0
      n_env=n_env+1
      Do i=1, nprims
      if(info(i)/=0)then
      bas(1,j+1)=ipat(i)-1        ! atom id starts with zero
      if(ipty(i)==1)bas(2,j+1)=0  ! orbital quantum number l
      if(ipty(i)==2)bas(2,j+1)=1
      if(ipty(i)==5)bas(2,j+1)=2
      if(ipty(i)==11)bas(2,j+1)=3
      bas(3,j+1)=info(i)         ! number of primitives in each contraction
      bas(4,j+1)=1               ! one contraction, no more
      bas(6,j+1)=n_env-1         ! "zero" for list of exponents in env
      bas(7,j+1)=n_env+info(i)-1 ! "zero" for list of contraction coefficients in env
      n_env=n_env+info(i)*2
      j=j+1
      endif
      enddo

      allocate(env(1:n_env)) ! start with 20 because variables are passed with values <=
      !env(1) threshold on integrals
      !env(1)=7.0d0
c center of nuclear charge and molar mass
!       sumwx=0.d0
!       sumwy=0.d0
!       sumwz=0.d0
!       sumw=0.0d0
!       xmolw=0.0d0
!       do i=1,ncent
!          atmass=co(i,4)
!          sumw=sumw+atmass
!          sumwx=sumwx+atmass*co(i,1)
!          sumwy=sumwy+atmass*co(i,2)
!          sumwz=sumwz+atmass*co(i,3)
!          xmolw=xmolw+ams(idint(atmass))
!       enddo
!       cen(1)=sumwx/sumw
!       cen(2)=sumwy/sumw
!       cen(3)=sumwz/sumw
      env(2:4)=0.0 ! env(PTR_COMMON_ORIG:PTR_COMMON_ORIG+2)


      ! copy all data in env vector
      j=21
      Do i=1,ncent
      env(j+1)=co(i,1)
      env(j+2)=co(i,2)
      env(j+3)=co(i,3)
      j=j+3
      enddo

      j=j+1
      n=1
      Do i=1,nprims
      if(info(i)/=0)then
      Do k=0,info(i)-1
      env(j)=exip(i+k)
      j=j+1
      enddo
      Do k=0,info(i)-1
      If(ipty(i)==1)then
      env(j)=cxip(i+k)*dsqrt(PI*4.0d0)    ! different normalization in libcint
      endif
      If(ipty(i)==2)then
      env(j)=cxip(i+k)*dsqrt(PI*4.0d0/3.0d0)
      endif
      If(ipty(i)==5)then
      env(j)=cxip(i+k)
      endif
      If(ipty(i)==11)then
      env(j)=cxip(i+k)
      endif
      j=j+1
      enddo
      n=n+1
      endif
      enddo

      allocate(di_all(nbas),sum_di(nbas))
      sum_di=0
      Do i=1, nbas
      di_all(i)=CINTcgto_cart(i-1,bas)
      if(i<nbas)sum_di(i+1)=sum_di(i)+di_all(i)
      enddo

      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
!$omp parallel private(i,j,k,l,shls,di,dj,buf)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj))
      buf=0.0
      call cint1e_ovlp_cart(buf,shls,atm,ncent,bas,nbas,env)
!       write(*,*)i,j
!       write(*,*)buf(:,1:dj)
      Do k=1,di
      Do l=1,dj
      overlap_AO(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l))
     .)=buf(k,l)*norm_cart(di,k)*norm_cart(dj,l)
!       write(*,*)k,l,buf(k,l)*norm_cart(di,k)*norm_cart(dj,l),
!      .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l))
      enddo
      enddo
      deallocate(buf)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di))
      call cint1e_ovlp_cart(buf,shls,atm,ncent,bas,nbas,env)
!       write(*,*)i,i
      Do k=1,di
      Do l=1,k-1
      overlap_AO(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l))
     .)=buf(k,l)*norm_cart(di,k)*norm_cart(di,l)
!       write(*,*)k,l,buf(k,l)*norm_cart(di,k)*norm_cart(di,l),
!      .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l))
      enddo
      overlap_AO(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k))
     .)=buf(k,k)*norm_cart(di,k)*norm_cart(di,k)
!       write(*,*)k,k,buf(k,k)*norm_cart(di,k)*norm_cart(di,k),
!      .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k))
      enddo
      deallocate(buf)
      enddo
!$omp end do
!$omp end parallel
      deallocate(info)
      end subroutine overlap

      subroutine dipole_moment(ncent,nprims,nbf,mu)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m,ii,jj
      integer :: ncent,nprims,nbf

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:)
      double precision :: mu(nbf*(nbf+1)/2,1:3)
      double precision :: norm_cart(1:10,1:10)

      integer*8 ::lin8


      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
!$omp parallel private(i,j,k,l,n,shls,di,dj,buf)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj*3))
      buf=0.0
      call cint1e_r_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do n=1,3
      Do k=1,di
      Do l=1,dj
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),n
     .)=-buf(k,l+(n-1)*dj)*norm_cart(di,k)*norm_cart(dj,l)
      enddo
      enddo
      enddo
      deallocate(buf)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di*3))
      call cint1e_r_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do n=1,3
      Do k=1,di
      Do l=1,k-1
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=-buf(k,l+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      enddo
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),n
     .)=-buf(k,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,k)
      enddo
      enddo
      deallocate(buf)
      enddo
!$omp end do
!$omp end parallel
      end subroutine dipole_moment

      subroutine magnetic_moment(ncent,nprims,nbf,mu)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m
      integer :: ncent,nprims,nbf

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:)
      double precision :: mu(nbf*(nbf+1)/2,1:3)
      double precision :: norm_cart(1:10,1:10)

      integer*8 ::lin8


      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
      mu=0.0d0
!$omp parallel private(i,j,k,l,n,shls,di,dj,buf)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj*3))
      buf=0.0
      call cint1e_cg_irxp_cart(buf,shls,atm,ncent,bas,nbas,env)
      !write(*,*)i,di,j,dj,buf
      Do n=1,3
      Do k=1,di
      Do l=1,dj
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),n
     .)=buf(k,l+(n-1)*dj)*norm_cart(di,k)*norm_cart(dj,l)
      enddo
      enddo
      enddo
      deallocate(buf)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di*3))
      call cint1e_cg_irxp_cart(buf,shls,atm,ncent,bas,nbas,env)

      Do n=1,3
      Do k=1,di
      Do l=1,k-1
      if(di==6.and.l==2.and.k==4)then !incorrect signed element bug with skew-symmetric matrix in libcint
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==6.and.l==3.and.k==6)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==6.and.l==5.and.k==6)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==2.and.k==10)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==3.and.k==7)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==4.and.k==7)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==4.and.k==10)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==6.and.k==7)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==6.and.k==10)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==6.and.k==8)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==5.and.k==8)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==5.and.k==6)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==5.and.k==9)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else if(di==10.and.l==9.and.k==10)then
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(l,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      else
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(k,l+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      endif
      enddo
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),n
     .)=buf(k,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,k)
      enddo
      enddo
      deallocate(buf)
      enddo
!$omp end do
!$omp end parallel
      end subroutine magnetic_moment

      subroutine magnetic_moment_other(ncent,nprims,nbf,mu)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m
      integer :: ncent,nprims,nbf

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:),buf1(:,:)
      double precision :: mu(nbf*(nbf+1)/2,1:3)
      double precision :: norm_cart(1:10,1:10)

      integer*8 ::lin8


      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
      mu=0.0d0
!$omp parallel private(i,j,k,l,n,shls,di,dj,buf,buf1)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj*9),buf1(di,dj*3))
      call cint1e_irp_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do k=1,di
      Do l=1,dj
      buf1(k,l+(0)*dj)=buf(k,l+(5)*dj)-buf(k,l+(7)*dj)
      buf1(k,l+(1)*dj)=-buf(k,l+(2)*dj)+buf(k,l+(6)*dj)
      buf1(k,l+(2)*dj)=buf(k,l+(1)*dj)-buf(k,l+(3)*dj)
      enddo
      enddo

      Do n=1,3
      Do k=1,di
      Do l=1,dj
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),n
     .)=buf1(k,l+(n-1)*dj)*norm_cart(di,k)*norm_cart(dj,l)
      enddo
      enddo
      enddo
      deallocate(buf,buf1)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di*9),buf1(di,di*3))
      call cint1e_irp_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do k=1,di
      Do l=1,di
      buf1(k,l+(0)*di)=buf(k,l+(5)*di)-buf(k,l+(7)*di)
      buf1(k,l+(1)*di)=-buf(k,l+(2)*di)+buf(k,l+(6)*di)
      buf1(k,l+(2)*di)=buf(k,l+(1)*di)-buf(k,l+(3)*di)
      enddo
      enddo
      if(di==6)then
      Do n=0,2
      Do k=1,di
      write(*,'(6f10.3)')buf1(k,1+n*di),buf1(k,2+n*di),buf1(k,3+n*di),
     .buf1(k,4+n*di),buf1(k,5+n*di),buf1(k,6+n*di)
      enddo
            write(*,*)'**',n
      enddo

      Do n=0,8
      Do k=1,di
      write(*,'(6f10.3)')buf(k,1+n*di),buf(k,2+n*di),buf(k,3+n*di),
     .buf(k,4+n*di),buf(k,5+n*di),buf(k,6+n*di)
      enddo
            write(*,*)'**',n
      enddo

      endif

      Do n=1,3
      Do k=1,di
      Do l=1,k-1
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf1(k,l+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      enddo
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),n
     .)=buf1(k,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,k)
      enddo
      enddo
      deallocate(buf,buf1)
      enddo
!$omp end do
!$omp end parallel
      end subroutine magnetic_moment_other

      subroutine velo_moment(ncent,nprims,nbf,mu)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m
      integer :: ncent,nprims,nbf

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:)
      double precision :: mu(nbf*(nbf+1)/2,1:3)
      double precision :: norm_cart(1:10,1:10)

      double precision,external :: CINTgto_norm
      integer,external :: CINTcgto_cart

      integer*8 ::lin8


      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
!$omp parallel private(i,j,k,l,n,shls,di,dj,buf)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj*3))
      buf=0.0
      call cint1e_p_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do n=1,3
      Do k=1,di
      Do l=1,dj
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),n
     .)=buf(k,l+(n-1)*dj)*norm_cart(di,k)*norm_cart(dj,l)
      enddo
      enddo
      enddo
      deallocate(buf)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di*3))
      call cint1e_p_cart(buf,shls,atm,ncent,bas,nbas,env)
      Do n=1,3
      Do k=1,di
      Do l=1,k-1
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),n
     .)=buf(k,l+(n-1)*di)*norm_cart(di,k)*norm_cart(di,l)
      enddo
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),n
     .)=buf(k,k+(n-1)*di)*norm_cart(di,k)*norm_cart(di,k)
      enddo
      enddo
      deallocate(buf)
      enddo
!$omp end do
!$omp end parallel
      end subroutine velo_moment

      subroutine quadrupole_moment(ncent,nprims,nbf,mu)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,n,m
      integer :: ncent,nprims,nbf

      integer :: shls(2), di,dj,dk,dl, orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:)
      double precision :: mu(nbf*(nbf+1)/2,1:6)! XX,YY,ZZ,XY,XZ,YZ
      double precision :: norm_cart(1:10,1:10)

      double precision,external :: CINTgto_norm
      integer,external :: CINTcgto_cart

      integer*8 ::lin8

      ! To get the same order as in libcint
      orb_cart(1,1)=1
      orb_cart(3,1)=1
      orb_cart(3,2)=2
      orb_cart(3,3)=3
      orb_cart(6,1)=1
      orb_cart(6,2)=4
      orb_cart(6,3)=5
      orb_cart(6,4)=2
      orb_cart(6,5)=6
      orb_cart(6,6)=3
      orb_cart(10,1)=1
      orb_cart(10,2)=4
      orb_cart(10,3)=5
      orb_cart(10,4)=6
      orb_cart(10,5)=10
      orb_cart(10,6)=8
      orb_cart(10,7)=2
      orb_cart(10,8)=7
      orb_cart(10,9)=9
      orb_cart(10,10)=3
      ! For cartesian basis libcint uses d and f not normalized
      if(spherical)then
      norm_cart=1.0
      else
        norm_cart=1.0d0
        norm_cart(6,1)=1.0d0
        norm_cart(6,2)=dsqrt(3.0d0)
        norm_cart(6,3)=dsqrt(3.0d0)
        norm_cart(6,4)=1.0d0
        norm_cart(6,5)=dsqrt(3.0d0)
        norm_cart(6,6)=1.0d0
        norm_cart(10,1)=1.0d0
        norm_cart(10,2)=dsqrt(5.0d0)
        norm_cart(10,3)=dsqrt(5.0d0)
        norm_cart(10,4)=dsqrt(5.0d0)
        norm_cart(10,5)=dsqrt(15.0d0)
        norm_cart(10,6)=dsqrt(5.0d0)
        norm_cart(10,7)=1.0d0
        norm_cart(10,8)=dsqrt(5.0d0)
        norm_cart(10,9)=dsqrt(5.0d0)
        norm_cart(10,10)=1.0d0
      endif
!$omp parallel private(i,j,k,l,shls,di,dj,buf)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)

      Do j=1,i-1
      shls(2)=j-1 ; dj=di_all(j)
      allocate(buf(di,dj*9))!XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
      call cint1e_rr_cart(buf,shls,atm,ncent,bas,nbas,env)

      Do k=1,di
      Do l=1,dj
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),1
     .)=buf(k,l+(1-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !XX
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),2
     .)=buf(k,l+(5-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !YY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),3
     .)=buf(k,l+(9-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !ZZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),4
     .)=buf(k,l+(2-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !XY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),5
     .)=buf(k,l+(3-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !XZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(j)+orb_cart(dj,l)),6
     .)=buf(k,l+(6-1)*dj)*norm_cart(di,k)*norm_cart(dj,l) !YZ
      enddo
      enddo

      deallocate(buf)
      enddo

      shls(2)=i-1 ; dj=di
      allocate(buf(di,di*9))
      call cint1e_rr_cart(buf,shls,atm,ncent,bas,nbas,env)

      Do k=1,di
      Do l=1,k-1
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),1
     .)=buf(k,l+(1-1)*di)*norm_cart(di,k)*norm_cart(di,l) !XX
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),2
     .)=buf(k,l+(5-1)*di)*norm_cart(di,k)*norm_cart(di,l) !YY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),3
     .)=buf(k,l+(9-1)*di)*norm_cart(di,k)*norm_cart(di,l) !ZZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),4
     .)=buf(k,l+(2-1)*di)*norm_cart(di,k)*norm_cart(di,l) !XY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),5
     .)=buf(k,l+(3-1)*di)*norm_cart(di,k)*norm_cart(di,l) !XZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,l)),6
     .)=buf(k,l+(6-1)*di)*norm_cart(di,k)*norm_cart(di,l) !YZ
      enddo
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),1
     .)=buf(k,k+(1-1)*di)*norm_cart(di,k)*norm_cart(di,k) !XX
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),2
     .)=buf(k,k+(5-1)*di)*norm_cart(di,k)*norm_cart(di,k) !YY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),3
     .)=buf(k,k+(9-1)*di)*norm_cart(di,k)*norm_cart(di,k) !ZZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),4
     .)=buf(k,k+(2-1)*di)*norm_cart(di,k)*norm_cart(di,k) !XY
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),5
     .)=buf(k,k+(3-1)*di)*norm_cart(di,k)*norm_cart(di,k) !XZ
      mu(
     .lin8(sum_di(i)+orb_cart(di,k),sum_di(i)+orb_cart(di,k)),6
     .)=buf(k,k+(6-1)*di)*norm_cart(di,k)*norm_cart(di,k) !YZ
      enddo

      deallocate(buf)
      enddo
!$omp end do
!$omp end parallel
      end subroutine quadrupole_moment
