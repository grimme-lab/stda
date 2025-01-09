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
! Correct TDA eigenvector for Rv  c
!      (A+0.5*B)/omega * X        c
ccccccccccccccccccccccccccccccccccc
      subroutine apbtrafo(n,nroot,x,e,xl,yl,zl,xv,yv,zv,xm,ym,zm,xmss
     .                  ,maxconf,iconf,ak,rvpout)
      use commonlogicals
      implicit none
      integer, intent ( in ) :: maxconf,n,nroot,iconf(maxconf,2)
      real*8, intent( in ) :: xl(*),yl(*),zl(*)
      real*8, intent( in ) :: xv(*),yv(*),zv(*)
      real*8, intent( in ) :: xm(*),ym(*),zm(*)
      real*4, intent( in ) :: x(n,n),e(n)
      real*8, intent( in ) :: xmss
      real*8 , intent( out ) :: rvpout(nroot)

      integer i,j,k,l,io,iv
      integer*8 lin8,ij
      real*8 de,ef,ak,xp,xlp,ylp,zlp,xvp,yvp,zvp,xmp,ymp,zmp
      real*8 xvp2,yvp2,zvp2,xmp2,ymp2,zmp2,efmod
      real*8 flp,fvp,rlp,rvp,p23,fact,unew,uold,hilf !,enew(nroot)
      real*8 flx,fly,flz,fvx,fvy,fvz,rlx,rly,rlz,rvx,rvy,rvz ! resolved along the three orientations
      real*4, allocatable :: bmat(:),upper(:,:), xnew(:,:)
      write(*,'(A)',advance='yes') '  perform velo correction for X...'
      ij=n
      ij=ij*(ij+1)/2
      allocate(bmat(ij))

************** read B matrix (packed) **************************
      open(unit=52,file='bmat',form='unformatted',status='old')
      read(52)bmat
      close(52,status='delete')
******************************************************************

************** blow up 0.5*B and compute (0.5*B)*X *******************
      allocate(upper(n,n))
      call spack2tri(n,bmat,upper) ! blow up matrix, but we only need upper triangle
      deallocate(bmat)
      allocate(xnew(n,nroot))
      call ssymm('l','u',n,nroot,1.0e0,upper,n,x,n,0.e0,xnew,n)
      deallocate(upper)
******************************************************************
!
!******************** scale with omega_TDA ************************
!   compute divide by omega_TDA to yield  X_new
      do i=1,nroot
         ef=1.0d0/(dble(e(i))+1.0d-8)
         do j=1,n
            xnew(j,i)=ef*xnew(j,i)
         enddo
      enddo
!******************************************************************
!
!!! OLD ORTHOGONALIZATION PART - NOT USED ANYMORE
!
************** compute overlap: S = (X_new)**T * (X_new) *********
!      allocate(upper(nroot,nroot))
!      upper=0.0
!      ! compute overlap between xpy and xnew
!      call sgemm('T','n',nroot,nroot,n,1.0,xnew,n,xnew,n,0.0,
!     .           upper,nroot)
!
!!      enew=0.0d0
!!      do i=1,nroot
!!         enew(i)=dble(e(i))/dble(upper(i,i))
!!      enddo
!******************************************************************
!!
!************** do Cholesky - get upper triangular: S=U**T * U ****
!      call spotrf('u', nroot, upper, nroot, l )
!      if(l.ne.0) stop 'error in Cholesky decomposition'
!******************************************************************
!!
!************** orthogonalize vector: X' = X_new * U**-1 **********
!      call strsm('R','U','N','N',n, nroot,1.0,upper,nroot,xnew,n)
!*****************************************************************
!
!*********** check orthogonality *****************
!
!         upper=0.0
!         call sgemm('T','n',nroot,nroot,n,1.0,xnew,n,xnew,n,0.0,
!     .              upper,nroot)
!
!
!         do i=1,min(12,nroot)
!          write(*,'(12f10.6)') (upper(j,i),j=1,min(12,nroot))
!         enddo
!         write(*,*)
!         do i=max(1,nroot-11),nroot
!          write(*,'(12f10.6)') (upper(j,i),j=max(1,nroot-11),nroot)
!         enddo
!************************************************
!!!!
!      deallocate(upper)

      write(*,*)' writing trafoed spectral data to tda.dat ...'
      open(unit=28,file='tda.dat',status='replace')
      write(28,*)'NM'
      write(28,*)'VELO'
      write(28,*)'MMASS'
      write(28,*)xmss
      write(28,*)'LFAKTOR'
      write(28,*)' 0.5'
      write(28,*)'RFAKTOR'
      write(28,*)' 1.0'
      write(28,*)'WIDTH'
      write(28,*)' 0.20'
      write(28,*)'SHIFT'
      write(28,*)' 0.00'
      write(28,*)'DATXY'
      if(aniso) then
        write(*,*)'writing anisotropic data ...'
        write(*,*)'...tdax.dat,tday.dat,tdaz.dat'
        open(unit=25,file='tdax.dat')
        write(25,*)'NM'
        write(25,*)'VELO'
        write(25,*)'MMASS'
        write(25,*)xmss
        write(25,*)'LFAKTOR'
        write(25,*)' 0.5'
        write(25,*)'RFAKTOR'
        write(25,*)' 1.0'
        write(25,*)'WIDTH'
        write(25,*)' 0.20'
        write(25,*)'SHIFT'
        write(25,*)' 0.00'
        write(25,*)'DATXY'
        open(unit=26,file='tday.dat')
        write(26,*)'NM'
        write(26,*)'VELO'
        write(26,*)'MMASS'
        write(26,*)xmss
        write(26,*)'LFAKTOR'
        write(26,*)' 0.5'
        write(26,*)'RFAKTOR'
        write(26,*)' 1.0'
        write(26,*)'WIDTH'
        write(26,*)' 0.20'
        write(26,*)'SHIFT'
        write(26,*)' 0.00'
        write(26,*)'DATXY'
        open(unit=27,file='tdaz.dat')
        write(27,*)'NM'
        write(27,*)'VELO'
        write(27,*)'MMASS'
        write(27,*)xmss
        write(27,*)'LFAKTOR'
        write(27,*)' 0.5'
        write(27,*)'RFAKTOR'
        write(27,*)' 1.0'
        write(27,*)'WIDTH'
        write(27,*)' 0.20'
        write(27,*)'SHIFT'
        write(27,*)' 0.00'
        write(27,*)'DATXY'
      endif

      p23=ak* 2.0d0/3.0d0
      do i=1,nroot
! A+B transformed stuff
         de=dble(e(i))
         ef=1.0d0/(de+1.0d-8)
! damp exponent from -1 to 0 for small energies
         hilf=1.0d0-exp(-150.0d0*de*de)
         efmod=1.0d0/((de**hilf)+1.0d-8)
         xlp=0.0d0
         ylp=0.0d0
         zlp=0.0d0
         xmp=0.0d0
         ymp=0.0d0
         zmp=0.0d0
         xvp=0.0d0
         yvp=0.0d0
         zvp=0.0d0
         rlp=0.0d0
         rvp=0.0d0
         flp=0.0d0
         fvp=0.0d0
         flx=0.0d0
         fly=0.0d0
         flz=0.0d0
         fvx=0.0d0
         fvy=0.0d0
         fvz=0.0d0
         rlx=0.0d0
         rly=0.0d0
         rlz=0.0d0
         rvx=0.0d0
         rvy=0.0d0
         rvz=0.0d0
         xmp2=0.0d0
         ymp2=0.0d0
         zmp2=0.0d0
         xvp2=0.0d0
         yvp2=0.0d0
         zvp2=0.0d0
         do j=1,n
           io=iconf(j,1)
           iv=iconf(j,2)
           unew=dble(xnew(j,i))
           uold=dble(x(j,i))
           ij=lin8(io,iv)
! A+B transformed stuff
           xlp=xlp+xl(ij)*uold
           ylp=ylp+yl(ij)*uold
           zlp=zlp+zl(ij)*uold
           xvp=xvp+xv(ij)*uold
           yvp=yvp+yv(ij)*uold
           zvp=zvp+zv(ij)*uold
           xmp=xmp+xm(ij)*uold
           ymp=ymp+ym(ij)*uold
           zmp=zmp+zm(ij)*uold
           xvp2=xvp2+xv(ij)*unew
           yvp2=yvp2+yv(ij)*unew
           zvp2=zvp2+zv(ij)*unew
           xmp2=xmp2+xm(ij)*unew
           ymp2=ymp2+ym(ij)*unew
           zmp2=zmp2+zm(ij)*unew
         enddo
! A+0.5*B transformed stuff
         ! fL stays the same
           flx=xlp*xlp
           fly=ylp*ylp
           flz=zlp*zlp
           xp=flx+fly+flz
           flp=xp*p23*de
         ! fV stays the same
           fvx=xvp*xvp
           fvy=yvp*yvp
           fvz=zvp*zvp
           xp=fvx+fvy+fvz
           fvp=p23*ef*xp
! RL is the same
           rlx=xlp*xmp
           rly=ylp*ymp
           rlz=zlp*zmp
           xp=rlx+rly+rlz
           rlp=-235.7220d0*xp*ak
! Rv gets (symmetric) corrections from B (*vp2 and *mp2 terms)
           rvx=xvp*xmp+xvp2*xmp+xvp*xmp2
           rvy=yvp*ymp+yvp2*ymp+yvp*ymp2
           rvz=zvp*zmp+zvp2*zmp+zvp*zmp2
           xp=rvx+rvy+rvz
           rvp=-235.7220d0*xp*ak*efmod
! print out transformed stuff
!           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           rvpout(i)=rvp
           if(aniso)then
             hilf=de*2.0d0*ak
             flx=flx*hilf
             fly=fly*hilf
             flz=flz*hilf
             hilf=ef*2.0d0*ak
             fvx=fvx*hilf
             fvy=fvy*hilf
             fvz=fvz*hilf
             hilf=-235.7220d0*ak
             rlx=rlx*hilf
             rly=rly*hilf
             rlz=rlz*hilf
             hilf=-235.7220d0*ak*efmod
             rvx=rvx*hilf
             rvy=rvy*hilf
             rvz=rvz*hilf
             write(25,'(i4,F10.4,4f13.6)')i,de*27.21139,flx,fvx,rlx,rvx
             write(26,'(i4,F10.4,4f13.6)')i,de*27.21139,fly,fvy,rly,rvy
             write(27,'(i4,F10.4,4f13.6)')i,de*27.21139,flz,fvz,rlz,rvz
           endif
      enddo
      close(28)
      if(aniso)then
        close(25)
        close(26)
        close(27)
        write(*,*)'data given such that f = (f_x + f_y + f_z)/3'
        write(*,*)'                 and R = R_x + R_y + R_z'
      endif
      deallocate(xnew)
      end subroutine apbtrafo



      subroutine spack2tri(n,a,b)
c     blow up symmetric matrix to its lower (upper in Lapack) triangular form
      implicit real*4 (a-h,o-z)
      real*4 a(n*(n+1)/2),b(n,n)
      integer*8 ij
      ij=0
      do i=1,n
         do j=1,i
            ij=ij+1
            b(j,i)=a(ij)
!            b(i,j)=a(ij)
         enddo
      enddo
      return
      end

      subroutine apbtrafo_uks(n,na,nb,nroot,x,e,xla,yla,zla,xva,yva,zva,
     . xma,yma,zma,xlb,ylb,zlb,xvb,yvb,zvb,xmb,ymb,zmb,xmss,maxconfa,
     . maxconfb,iconfa,iconfb,rvpout)
      use commonlogicals
      implicit none
      integer, intent ( in ) :: maxconfa,maxconfb,n,na,nb,nroot
      integer, intent ( in ) :: iconfa(maxconfa,2),iconfb(maxconfb,2)
      real*8, intent(in) :: xla(*),yla(*),zla(*)
      real*8, intent(in) :: xva(*),yva(*),zva(*)
      real*8, intent(in) :: xma(*),yma(*),zma(*)
      real*8, intent(in) :: xlb(*),ylb(*),zlb(*)
      real*8, intent(in) :: xvb(*),yvb(*),zvb(*)
      real*8, intent(in) :: xmb(*),ymb(*),zmb(*)
      real*4, intent(in) :: x(n,n),e(n)
      real*8, intent(in) :: xmss
      real*8, intent( out ) :: rvpout(nroot)

      integer i,j,k,l,io,iv
      integer*8 ij,lin8
      real*8 de,ef,xp,xlp,ylp,zlp,xvp,yvp,zvp,xmp,ymp,zmp
      real*8 xvp2,yvp2,zvp2,xmp2,ymp2,zmp2,hilf,efmod
      real*8 flp,fvp,rlp,rvp,p23,fact,unew,uold !,enew(nroot)
      real*8 flx,fly,flz,fvx,fvy,fvz,rlx,rly,rlz,rvx,rvy,rvz ! resolved along the three orientations
      real*4, allocatable :: bmat(:),upper(:,:), xnew(:,:)
      write(*,'(A)',advance='yes') '  perform velo correction for X...'
      rvpout=0.0d0
      ij=n
      ij=ij*(ij+1)/2
      allocate(bmat(ij))

************** read 0.5*B matrix (packed) ***********************
      open(unit=52,file='bmat',form='unformatted',status='old')
      read(52)bmat
      close(52,status='delete')
******************************************************************

************** blow up 0.5*B and compute (0.5*B)*X ***************
      allocate(upper(n,n))
      call spack2tri(n,bmat,upper) ! blow up matrix, but we only need upper triangle
      deallocate(bmat)
      allocate(xnew(n,nroot))
      call ssymm('l','u',n,nroot,1.e0,upper,n,x,n,0.e0,xnew,n)
      deallocate(upper)
******************************************************************
!
******************** scale with omega_TDA ************************
!   divide by omega_TDA to yield  X_new
      do i=1,nroot
         ef=1.0d0/(dble(e(i))+1.0d-8)
         do j=1,n
            xnew(j,i)=ef*xnew(j,i)
         enddo
      enddo
******************************************************************
!
! testwise: orthogonalize the vectors
!
!************** compute overlap: S = (X_new)**T * (X_new) *********
!      allocate(upper(nroot,nroot))
!      upper=0.0
!      ! compute overlap between xpy and xnew
!      call sgemm('T','n',nroot,nroot,n,1.0,xnew,n,xnew,n,0.0,
!     .           upper,nroot)
!
!******************************************************************
!!
!************** do Cholesky - get upper triangular: S=U**T * U ****
!      call spotrf('u', nroot, upper, nroot, l )
!      if(l.ne.0) stop 'error in Cholesky decomposition'
!******************************************************************
!!
!************** orthogonalize vector: X' = X_new * U**-1 **********
!      call strsm('R','U','N','N',n, nroot,1.0,upper,nroot,xnew,n)
!******************************************************************
!!
!*********** check orthogonality *****************
!!
!!         upper=0.0
!!         call sgemm('T','n',nroot,nroot,n,1.0,xnew,n,xnew,n,0.0,
!!     .              upper,nroot)
!!
!!
!!         do i=1,min(12,nroot)
!!          write(*,'(12f10.6)') (upper(j,i),j=1,min(12,nroot))
!!         enddo
!!         write(*,*)
!!         do i=max(1,nroot-11),nroot
!!          write(*,'(12f10.6)') (upper(j,i),j=max(1,nroot-11),nroot)
!!         enddo
!************************************************
!!!!
!      deallocate(upper)

      write(*,*)' writing trafoed spectral data to tda.dat ...'
      open(unit=28,file='tda.dat',status='replace')
      write(28,*)'NM'
      write(28,*)'VELO'
      write(28,*)'MMASS'
      write(28,*)xmss
      write(28,*)'LFAKTOR'
      write(28,*)' 0.5'
      write(28,*)'RFAKTOR'
      write(28,*)' 1.0'
      write(28,*)'WIDTH'
      write(28,*)' 0.20'
      write(28,*)'SHIFT'
      write(28,*)' 0.00'
      write(28,*)'DATXY'

      p23=2.0d0/3.0d0
      do i=1,nroot
! A+B transformed stuff
         de=dble(e(i))
         ef=1.0d0/(de+1.0d-8)
! damp exponent from -1 to 0 for small energies
         hilf=1.0d0-exp(-150.0d0*de*de)
         efmod=1.0d0/((de**hilf)+1.0d-8)
         xlp=0.0d0
         ylp=0.0d0
         zlp=0.0d0
         xmp=0.0d0
         ymp=0.0d0
         zmp=0.0d0
         xvp=0.0d0
         yvp=0.0d0
         zvp=0.0d0
         rlp=0.0d0
         rvp=0.0d0
         flp=0.0d0
         fvp=0.0d0
         flx=0.0d0
         fly=0.0d0
         flz=0.0d0
         fvx=0.0d0
         fvy=0.0d0
         fvz=0.0d0
         rlx=0.0d0
         rly=0.0d0
         rlz=0.0d0
         rvx=0.0d0
         rvy=0.0d0
         rvz=0.0d0
         xmp2=0.0d0
         ymp2=0.0d0
         zmp2=0.0d0
         xvp2=0.0d0
         yvp2=0.0d0
         zvp2=0.0d0
!alpha
         do j=1,na
           io=iconfa(j,1)
           iv=iconfa(j,2)
           unew=dble(xnew(j,i))
           uold=dble(x(j,i))
           ij=lin8(io,iv)
! A+B transformed stuff
           xlp=xlp+xla(ij)*uold
           ylp=ylp+yla(ij)*uold
           zlp=zlp+zla(ij)*uold
           xvp=xvp+xva(ij)*uold
           yvp=yvp+yva(ij)*uold
           zvp=zvp+zva(ij)*uold
           xmp=xmp+xma(ij)*uold
           ymp=ymp+yma(ij)*uold
           zmp=zmp+zma(ij)*uold
           xvp2=xvp2+xva(ij)*unew
           yvp2=yvp2+yva(ij)*unew
           zvp2=zvp2+zva(ij)*unew
           xmp2=xmp2+xma(ij)*unew
           ymp2=ymp2+yma(ij)*unew
           zmp2=zmp2+zma(ij)*unew
         enddo
!beta
         do j=1,nb
           k=na+j
           io=iconfb(j,1)
           iv=iconfb(j,2)
           unew=dble(xnew(k,i))
           uold=dble(x(k,i))
           ij=lin8(io,iv)
! A+B transformed stuff
           xlp=xlp+xlb(ij)*uold
           ylp=ylp+ylb(ij)*uold
           zlp=zlp+zlb(ij)*uold
           xvp=xvp+xvb(ij)*uold
           yvp=yvp+yvb(ij)*uold
           zvp=zvp+zvb(ij)*uold
           xmp=xmp+xmb(ij)*uold
           ymp=ymp+ymb(ij)*uold
           zmp=zmp+zmb(ij)*uold
           xvp2=xvp2+xvb(ij)*unew
           yvp2=yvp2+yvb(ij)*unew
           zvp2=zvp2+zvb(ij)*unew
           xmp2=xmp2+xmb(ij)*unew
           ymp2=ymp2+ymb(ij)*unew
           zmp2=zmp2+zmb(ij)*unew
         enddo
! A+B transformed stuff
         ! fL stays the same
           flx=xlp*xlp
           fly=ylp*ylp
           flz=zlp*zlp
           xp=flx+fly+flz
           flp=xp*p23*de
         ! fV stays the same
           fvx=xvp*xvp
           fvy=yvp*yvp
           fvz=zvp*zvp
           xp=fvx+fvy+fvz
           fvp=p23*ef*xp
! RL is the same
           rlx=xlp*xmp
           rly=ylp*ymp
           rlz=zlp*zmp
           xp=rlx+rly+rlz
           rlp=-235.7220d0*xp
! Rv gets (symmetric) corrections from B (*vp2 and *mp2 terms)
           rvx=xvp*xmp+xvp2*xmp+xvp*xmp2
           rvy=yvp*ymp+yvp2*ymp+yvp*ymp2
           rvz=zvp*zmp+zvp2*zmp+zvp*zmp2
           xp=rvx+rvy+rvz
           rvp=-235.7220d0*xp*efmod
! print out transformed stuff
           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           rvpout(i)=rvp
           if(aniso)then
             hilf=de*2.0d0
             flx=flx*hilf
             fly=fly*hilf
             flz=flz*hilf
             hilf=ef*2.0d0
             fvx=fvx*hilf
             fvy=fvy*hilf
             fvz=fvz*hilf
             hilf=-235.7220d0
             rlx=rlx*hilf
             rly=rly*hilf
             rlz=rlz*hilf
             hilf=-235.7220d0*efmod
             rvx=rvx*hilf
             rvy=rvy*hilf
             rvz=rvz*hilf
             write(25,'(i4,F10.4,4f13.6)')i,de*27.21139,flx,fvx,rlx,rvx
             write(26,'(i4,F10.4,4f13.6)')i,de*27.21139,fly,fvy,rly,rvy
             write(29,'(i4,F10.4,4f13.6)')i,de*27.21139,flz,fvz,rlz,rvz
           endif
      enddo
      close(28)
      if(aniso)then
        close(25)
        close(26)
        close(29)
        write(*,*)'data given such that f = (f_x + f_y + f_z)/3'
        write(*,*)'                 and R = R_x + R_y + R_z'
      endif
      deallocate(xnew)
      end subroutine apbtrafo_uks


***********************************************************************
* set up 0.5*B  (packed form) in RKS case
***********************************************************************
      subroutine rtdacorr(nci,ncent,no,nv,mxcnf,iconf,dak,dax
     .                    ,ed,pia,qia,pij,qab)
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*4, intent(in)  :: qia(ncent,mxcnf),pia(ncent,mxcnf)
      real*4, intent(in)  :: pij(ncent,no*(no+1)/2)
      real*4, intent(in)  :: qab(ncent,nv*(nv+1)/2)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, allocatable :: qj(:),qk(:),bmat(:)
      integer i,j,io,iv,jo,jv,ierr,iiv,jjv,iwrk,jwrk
      integer*8 ij,lin8
      real*4 ek,ej,sdot,ak,ax,de,fact
      ij=nci
      ij=ij*(ij+1)/2
      allocate(qj(ncent),qk(ncent),bmat(ij), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj/bmat crashed'
      ak=real(dak)
      ax=real(dax)
! calculate 0.5*B
      bmat=0.0e0
      fact=0.50d0 ! this is the scaling of the B-contribution
      open(unit=52,file='bmat',form='unformatted',status='replace')
      ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
      do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           qk(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin8(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv + jjv
              ek=sdot(ncent,qk,1,qia(1,jwrk),1) ! ek = (ia|bj)
              bmat(ij)=(fact)*ak*ek
              jwrk=(io-1)*nv+jjv
              qj(1:ncent)=pia(1:ncent,jwrk)
              jwrk=(jo-1)*nv+iiv
              ek=sdot(ncent,qj,1,qia(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              bmat(ij)=bmat(ij)-fact*ax*ek ! scaled by ax
           enddo
           ij=lin8(i,i)
           ek=sdot(ncent,qk,1,qia(1,iwrk),1)
           bmat(ij)=fact*(ak*ek-ax*ek) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
      write(52)bmat
      close(52)
      deallocate(bmat,qk,qj)
      return

      end subroutine rtdacorr
***********************************************************************



***********************************************************************
* set up 0.5*B  (packed form) in UKS case !
***********************************************************************
       subroutine utdacorr(nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,
     .                    mxcnfb,iconfa,iconfb,dax,piaa,qiaa,
     .                    piab,qiab,pija,qaba,pijb,qabb)
       use omp_lib
      implicit none
      integer, intent(in) :: nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa
      integer, intent(in) :: mxcnfb,iconfa(mxcnfa,2),iconfb(mxcnfb,2)
      real*4, intent(in)  :: qiaa(ncent,mxcnfa),piaa(ncent,mxcnfa)
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qaba(ncent,nva*(nva+1)/2)
      real*4, intent(in)  :: qiab(ncent,mxcnfb),piab(ncent,mxcnfb)
      real*4, intent(in)  :: pijb(ncent,nob*(nob+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*8, intent(in)  :: dax
      real*4, allocatable :: qj(:),qk(:),bmat(:)
      integer i,j,io,iv,jo,jv,ierr,nex,iiv,jjv,iwrk,jwrk
      integer*8 lin8,ij
      real*4 ek,ej,sdot,ax,de,fact
      nex=nexa+nexb
      ij=nex
      ij=ij*(ij+1)/2
      allocate(qj(ncent),qk(ncent),bmat(ij), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj/bmat crashed'
      ax=real(dax)
! calculate 0.5*B
      fact=0.50e0 ! this is the scaling of the B-contribution
      bmat=0.0e0
      open(unit=52,file='bmat',form='unformatted',status='replace')
      ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
c alpha-alpha block
      do i = 1,nexa
         io=iconfa(i,1)
         iv=iconfa(i,2)
         iiv=iv-noa
         iwrk=(io-1)*nva + iiv
         qk(1:ncent)=piaa(1:ncent,iwrk)
         do j=1,i-1
            ij=lin8(i,j)
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
            bmat(ij)=(fact)*ek
            jwrk=(io-1)*nva+jjv
            qj(1:ncent)=piaa(1:ncent,jwrk)
            jwrk=(jo-1)*nva+iiv
            ek=sdot(ncent,qj,1,qiaa(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
            bmat(ij)=bmat(ij)-fact*ax*ek
         enddo
         ij=lin8(i,i)
         ek=sdot(ncent,qk,1,qiaa(1,iwrk),1)
         bmat(ij)=fact*(ek-ax*ek) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
      ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
! beta...
      do i = nexa+1,nex
         io=iconfb(i-nexa,1)
         iv=iconfb(i-nexa,2)
         iiv=iv-nob
         iwrk=(io-1)*nvb + iiv
         qk(1:ncent)=piab(1:ncent,iwrk)
! ...alpha block
         do j = 1,nexa
            ij=lin8(i,j)
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
            bmat(ij)=(fact)*ek
         enddo
! ...beta block
         do j = nexa+1,i-1
            ij=lin8(i,j)
            jo=iconfb(j-nexa,1)
            jv=iconfb(j-nexa,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
            ek=sdot(ncent,qk,1,qiab(1,jwrk),1)
            bmat(ij)=(fact)*ek
            jwrk=(io-1)*nvb+jjv
            qj(1:ncent)=piab(1:ncent,jwrk)
            jwrk=(jo-1)*nvb+iiv
            ek=sdot(ncent,qj,1,qiab(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
            bmat(ij)=bmat(ij)-fact*ax*ek
         enddo
         ij=lin8(i,i)
         ek=sdot(ncent,qk,1,qiab(1,iwrk),1)
         bmat(ij)=fact*(ek-ax*ek) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
!      call prmat4(6,bmat,nexa+nexb,0,'A+ 0.5 * B')
      write(52)bmat
      close(52)
      deallocate(bmat,qj,qk)
      return

      end subroutine utdacorr
***********************************************************************


cccccccccccccccccccccccccccccccccccccc
! TDA eigenvector with exciton print c
cccccccccccccccccccccccccccccccccccccc
      subroutine apbtrafoexc(n,nroot,x,e,xl,yl,zl,xv,yv,zv,xm,ym,zm,xmss
     .                     ,no,nv,coc,ncent,qia,maxconf,iconf,ak,rvpout)
      use commonlogicals
      implicit none
      integer, intent ( in ) :: maxconf,n,nroot,iconf(maxconf,2)
      integer, intent ( in ) :: ncent,no,nv
      real*8, intent( in ) :: xl(*),yl(*),zl(*)
      real*8, intent( in ) :: xv(*),yv(*),zv(*)
      real*8, intent( in ) :: xm(*),ym(*),zm(*),coc(3)
      real*4, intent( in ) :: x(n,n),e(n),qia(ncent,n)
      real*8, intent( in ) :: xmss
      real*8 , intent( out ) :: rvpout(nroot)

      integer i,j,k,l,ij,io,iv,lin
      real*8 de,ef,ak,xp,xlp,ylp,zlp,xvp,yvp,zvp,xmp,ymp,zmp
      real*8 xvu,yvu,zvu,xmu,ymu,zmu,aksqrt
      real*8 xvp2,yvp2,zvp2,xmp2,ymp2,zmp2,efmod,xms,yms,zms
      real*8 flp,fvp,rlp,rvp,p23,fact,unew,uold,hilf !,enew(nroot)
      real*8 flx,fly,flz,fvx,fvy,fvz,rlx,rly,rlz,rvx,rvy,rvz
      real*4, allocatable :: bmat(:),upper(:,:), xnew(:,:),q1(:)
      write(*,'(A)',advance='yes') '  perform velo correction for X...'
      allocate(bmat(n*(n+1)/2))

************** read B matrix (packed) **************************
      open(unit=52,file='bmat',form='unformatted',status='old')
      read(52)bmat
      close(52,status='delete')
******************************************************************

************** blow up 0.5*B and compute (0.5*B)*X *******************
      allocate(upper(n,n))
      call spack2tri(n,bmat,upper) ! blow up matrix, but we only need upper triangle
      deallocate(bmat)
      allocate(xnew(n,nroot))
      call ssymm('l','u',n,nroot,1.0e0,upper,n,x,n,0.e0,xnew,n)
      deallocate(upper)
******************************************************************
!
!******************** scale with omega_TDA ************************
!   compute divide by omega_TDA to yield  X_new
      do i=1,nroot
         ef=1.0d0/(dble(e(i))+1.0d-8)
         do j=1,n
            xnew(j,i)=ef*xnew(j,i)
         enddo
      enddo
!******************************************************************

      write(*,*)' writing trafoed spectral data to tda.dat ...'
      open(unit=28,file='tda.dat',status='replace')
      write(28,*)'NM'
      write(28,*)'VELO'
      write(28,*)'MMASS'
      write(28,*)xmss
      write(28,*)'LFAKTOR'
      write(28,*)' 0.5'
      write(28,*)'RFAKTOR'
      write(28,*)' 1.0'
      write(28,*)'WIDTH'
      write(28,*)' 0.20'
      write(28,*)'SHIFT'
      write(28,*)' 0.00'
      write(28,*)'DATXY'

      if(aniso) then
        write(*,*)'writing anisotropic data ...'
        write(*,*)'...tdax.dat,tday.dat,tdaz.dat'
        open(unit=25,file='tdax.dat')
        write(25,*)'NM'
        write(25,*)'VELO'
        write(25,*)'MMASS'
        write(25,*)xmss
        write(25,*)'LFAKTOR'
        write(25,*)' 0.5'
        write(25,*)'RFAKTOR'
        write(25,*)' 1.0'
        write(25,*)'WIDTH'
        write(25,*)' 0.20'
        write(25,*)'SHIFT'
        write(25,*)' 0.00'
        write(25,*)'DATXY'
        open(unit=26,file='tday.dat')
        write(26,*)'NM'
        write(26,*)'VELO'
        write(26,*)'MMASS'
        write(26,*)xmss
        write(26,*)'LFAKTOR'
        write(26,*)' 0.5'
        write(26,*)'RFAKTOR'
        write(26,*)' 1.0'
        write(26,*)'WIDTH'
        write(26,*)' 0.20'
        write(26,*)'SHIFT'
        write(26,*)' 0.00'
        write(26,*)'DATXY'
        open(unit=29,file='tdaz.dat')
        write(29,*)'NM'
        write(29,*)'VELO'
        write(29,*)'MMASS'
        write(29,*)xmss
        write(29,*)'LFAKTOR'
        write(29,*)' 0.5'
        write(29,*)'RFAKTOR'
        write(29,*)' 1.0'
        write(29,*)'WIDTH'
        write(29,*)' 0.20'
        write(29,*)'SHIFT'
        write(29,*)' 0.00'
        write(29,*)'DATXY'
      endif

      allocate(q1(ncent))
      q1=0.0e0

      p23=2.0d0/3.0d0
      aksqrt=sqrt(ak)
      do i=1,nroot
! A+B transformed stuff
         de=dble(e(i))
         ef=1.0d0/(de+1.0d-8)
! damp exponent from -1 to 0 for small energies
         hilf=1.0d0-exp(-150.0d0*de*de)
         efmod=1.0d0/((de**hilf)+1.0d-8)
         xlp=0.0d0
         ylp=0.0d0
         zlp=0.0d0
         xmp=0.0d0
         ymp=0.0d0
         zmp=0.0d0
         xvp=0.0d0
         yvp=0.0d0
         zvp=0.0d0
         rlp=0.0d0
         rvp=0.0d0
         flp=0.0d0
         fvp=0.0d0
         xmp2=0.0d0
         ymp2=0.0d0
         zmp2=0.0d0
         xvp2=0.0d0
         yvp2=0.0d0
         zvp2=0.0d0
         xms=0.0d0
         yms=0.0d0
         zms=0.0d0
         xmu=0.0d0
         ymu=0.0d0
         zmu=0.0d0
         xvu=0.0d0
         yvu=0.0d0
         zvu=0.0d0
         flx=0.0d0
         fly=0.0d0
         flz=0.0d0
         fvx=0.0d0
         fvy=0.0d0
         fvz=0.0d0
         rlx=0.0d0
         rly=0.0d0
         rlz=0.0d0
         rvx=0.0d0
         rvy=0.0d0
         rvz=0.0d0
         q1=0.0e0
         do j=1,n
           io=iconf(j,1)
           iv=iconf(j,2)
           unew=dble(xnew(j,i))
           uold=dble(x(j,i))
           ij=lin(io,iv)
! A+B transformed stuff
           xlp=xlp+xl(ij)*uold
           ylp=ylp+yl(ij)*uold
           zlp=zlp+zl(ij)*uold
           xvp=xvp+xv(ij)*uold
           yvp=yvp+yv(ij)*uold
           zvp=zvp+zv(ij)*uold
           xmp=xmp+xm(ij)*uold
           ymp=ymp+ym(ij)*uold
           zmp=zmp+zm(ij)*uold
           xvp2=xvp2+xv(ij)*unew
           yvp2=yvp2+yv(ij)*unew
           zvp2=zvp2+zv(ij)*unew
           xmp2=xmp2+xm(ij)*unew
           ymp2=ymp2+ym(ij)*unew
           zmp2=zmp2+zm(ij)*unew
           l=(io-1)*nv+(iv-no)
           q1(1:ncent)=q1(1:ncent)+qia(1:ncent,l)*real(uold)
         enddo

         ! multiply with factor from spin-integration
         xlp=xlp*aksqrt
         ylp=ylp*aksqrt
         zlp=zlp*aksqrt
         xvp=xvp*aksqrt
         yvp=yvp*aksqrt
         zvp=zvp*aksqrt
         xmp=xmp*aksqrt
         ymp=ymp*aksqrt
         zmp=zmp*aksqrt
         xvp2=xvp2*aksqrt
         yvp2=yvp2*aksqrt
         zvp2=zvp2*aksqrt
         xmp2=xmp2*aksqrt
         ymp2=ymp2*aksqrt
         zmp2=zmp2*aksqrt

!------
! print exciton output
         xvu=xvp+xvp2
         yvu=yvp+yvp2
         zvu=zvp+zvp2
         xmu=xmp+xmp2
         ymu=ymp+ymp2
         zmu=zmp+zmp2
         xms=(yvu*coc(3)-zvu*coc(2))
         yms=(zvu*coc(1)-xvu*coc(3))
         zms=(xvu*coc(2)-yvu*coc(1))
         write(27,'(i5,a,x,F10.4,x,a2)') i,':',de*27.21139,'eV'
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'l:',xlp,ylp,zlp
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'p:',xvu,yvu,zvu
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'m:',xmu-xms,ymu-yms,
     .                                                zmu-zms
         do k=1,ncent
           write(27,'(f12.8)') real(aksqrt)*q1(k)
         enddo
!------
! A+0.5*B transformed stuff
         ! fL stays the same
           flx=xlp*xlp
           fly=ylp*ylp
           flz=zlp*zlp
           xp=flx+fly+flz
           flp=xp*p23*de
         ! fV stays the same
           fvx=xvp*xvp
           fvy=yvp*yvp
           fvz=zvp*zvp
           xp=fvx+fvy+fvz
           fvp=p23*ef*xp
! RL is the same
           rlx=xlp*xmp
           rly=ylp*ymp
           rlz=zlp*zmp
           xp=rlx+rly+rlz
           rlp=-235.7220d0*xp
! Rv gets (symmetric) corrections from B (*vp2 and *mp2 terms)
           rvx=xvp*xmp+xvp2*xmp+xvp*xmp2
           rvy=yvp*ymp+yvp2*ymp+yvp*ymp2
           rvz=zvp*zmp+zvp2*zmp+zvp*zmp2
           xp=rvx+rvy+rvz
           rvp=-235.7220d0*xp*efmod
! print out transformed stuff
!           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           rvpout(i)=rvp
           if(aniso)then
             hilf=de*2.0d0
             flx=flx*hilf
             fly=fly*hilf
             flz=flz*hilf
             hilf=ef*2.0d0
             fvx=fvx*hilf
             fvy=fvy*hilf
             fvz=fvz*hilf
             hilf=-235.7220d0
             rlx=rlx*hilf
             rly=rly*hilf
             rlz=rlz*hilf
             hilf=-235.7220d0*efmod
             rvx=rvx*hilf
             rvy=rvy*hilf
             rvz=rvz*hilf
             write(25,'(i4,F10.4,4f13.6)')i,de*27.21139,flx,fvx,rlx,rvx
             write(26,'(i4,F10.4,4f13.6)')i,de*27.21139,fly,fvy,rly,rvy
             write(29,'(i4,F10.4,4f13.6)')i,de*27.21139,flz,fvz,rlz,rvz
           endif
      enddo
      close(28)
      if(aniso)then
        close(25)
        close(26)
        close(29)
        write(*,*)'data given such that f = (f_x + f_y + f_z)/3'
        write(*,*)'                 and R = R_x + R_y + R_z'
      endif
      deallocate(xnew,q1)

      end subroutine apbtrafoexc

      subroutine apbtrafoexc_uks(n,na,nb,nroot,x,e,xla,yla,zla,xva,yva
     . ,zva,xma,yma,zma,xlb,ylb,zlb,xvb,yvb,zvb,xmb,ymb,zmb,xmss
     . ,noa,nva,nob,nvb,coc,ncent,qiaa,qiab,maxconfa,maxconfb,iconfa
     . ,iconfb,rvpout)
      use commonlogicals
      implicit none
      integer, intent ( in ) :: maxconfa,maxconfb,n,na,nb,nroot
      integer, intent ( in ) :: ncent,noa,nva,nob,nvb
      integer, intent ( in ) :: iconfa(maxconfa,2),iconfb(maxconfb,2)
      real*8, intent(in) :: xla(*),yla(*),zla(*)
      real*8, intent(in) :: xva(*),yva(*),zva(*)
      real*8, intent(in) :: xma(*),yma(*),zma(*)
      real*8, intent(in) :: xlb(*),ylb(*),zlb(*)
      real*8, intent(in) :: xvb(*),yvb(*),zvb(*)
      real*8, intent(in) :: xmb(*),ymb(*),zmb(*)
      real*4, intent(in) :: x(n,n),e(n),qiaa(ncent,na),qiab(ncent,nb)
      real*8, intent(in) :: xmss,coc(3)
      real*8, intent( out ) :: rvpout(nroot)

      integer i,j,k,l,ij,io,iv,lin
      real*8 de,ef,xp,xlp,ylp,zlp,xvp,yvp,zvp,xmp,ymp,zmp
      real*8 xmu,ymu,zmu,xvu,yvu,zvu,xms,yms,zms
      real*8 xvp2,yvp2,zvp2,xmp2,ymp2,zmp2,hilf,efmod
      real*8 flp,fvp,rlp,rvp,p23,fact,unew,uold !,enew(nroot)
      real*8 flx,fly,flz,fvx,fvy,fvz,rlx,rly,rlz,rvx,rvy,rvz ! resolved along the three orientations
      real*4, allocatable :: bmat(:),upper(:,:), xnew(:,:),q1(:)
      write(*,'(A)',advance='yes') '  perform velo correction for X...'
      rvpout=0.0d0

      allocate(bmat(n*(n+1)/2))

************** read 0.5*B matrix (packed) ***********************
      open(unit=52,file='bmat',form='unformatted',status='old')
      read(52)bmat
      close(52,status='delete')
******************************************************************

************** blow up 0.5*B and compute (0.5*B)*X ***************
      allocate(upper(n,n))
      call spack2tri(n,bmat,upper) ! blow up matrix, but we only need upper triangle
      deallocate(bmat)
      allocate(xnew(n,nroot))
      call ssymm('l','u',n,nroot,1.e0,upper,n,x,n,0.e0,xnew,n)
      deallocate(upper)
******************************************************************
!
******************** scale with omega_TDA ************************
!   divide by omega_TDA to yield  X_new
      do i=1,nroot
         ef=1.0d0/(dble(e(i))+1.0d-8)
         do j=1,n
            xnew(j,i)=ef*xnew(j,i)
         enddo
      enddo
******************************************************************

      write(*,*)' writing trafoed spectral data to tda.dat ...'
      open(unit=28,file='tda.dat',status='replace')
      write(28,*)'NM'
      write(28,*)'VELO'
      write(28,*)'MMASS'
      write(28,*)xmss
      write(28,*)'LFAKTOR'
      write(28,*)' 0.5'
      write(28,*)'RFAKTOR'
      write(28,*)' 1.0'
      write(28,*)'WIDTH'
      write(28,*)' 0.20'
      write(28,*)'SHIFT'
      write(28,*)' 0.00'
      write(28,*)'DATXY'

      allocate(q1(ncent))
      q1=0.0e0

      p23=2.0d0/3.0d0
      do i=1,nroot
! A+B transformed stuff
         de=dble(e(i))
         ef=1.0d0/(de+1.0d-8)
! damp exponent from -1 to 0 for small energies
         hilf=1.0d0-exp(-150.0d0*de*de)
         efmod=1.0d0/((de**hilf)+1.0d-8)
         xlp=0.0d0
         ylp=0.0d0
         zlp=0.0d0
         xmp=0.0d0
         ymp=0.0d0
         zmp=0.0d0
         xvp=0.0d0
         yvp=0.0d0
         zvp=0.0d0
         rlp=0.0d0
         rvp=0.0d0
         flp=0.0d0
         fvp=0.0d0
         flx=0.0d0
         fly=0.0d0
         flz=0.0d0
         fvx=0.0d0
         fvy=0.0d0
         fvz=0.0d0
         rlx=0.0d0
         rly=0.0d0
         rlz=0.0d0
         rvx=0.0d0
         rvy=0.0d0
         rvz=0.0d0
         xmp2=0.0d0
         ymp2=0.0d0
         zmp2=0.0d0
         xvp2=0.0d0
         yvp2=0.0d0
         zvp2=0.0d0
         q1=0.0e0
!alpha
         do j=1,na
           io=iconfa(j,1)
           iv=iconfa(j,2)
           unew=dble(xnew(j,i))
           uold=dble(x(j,i))
           ij=lin(io,iv)
! A+B transformed stuff
           xlp=xlp+xla(ij)*uold
           ylp=ylp+yla(ij)*uold
           zlp=zlp+zla(ij)*uold
           xvp=xvp+xva(ij)*uold
           yvp=yvp+yva(ij)*uold
           zvp=zvp+zva(ij)*uold
           xmp=xmp+xma(ij)*uold
           ymp=ymp+yma(ij)*uold
           zmp=zmp+zma(ij)*uold
           xvp2=xvp2+xva(ij)*unew
           yvp2=yvp2+yva(ij)*unew
           zvp2=zvp2+zva(ij)*unew
           xmp2=xmp2+xma(ij)*unew
           ymp2=ymp2+yma(ij)*unew
           zmp2=zmp2+zma(ij)*unew
           l=(io-1)*nva+(iv-noa)
           q1(1:ncent)=q1(1:ncent)+qiaa(1:ncent,l)*real(uold)
         enddo
!beta
         do j=1,nb
           k=na+j
           io=iconfb(j,1)
           iv=iconfb(j,2)
           unew=dble(xnew(k,i))
           uold=dble(x(k,i))
           ij=lin(io,iv)
! A+B transformed stuff
           xlp=xlp+xlb(ij)*uold
           ylp=ylp+ylb(ij)*uold
           zlp=zlp+zlb(ij)*uold
           xvp=xvp+xvb(ij)*uold
           yvp=yvp+yvb(ij)*uold
           zvp=zvp+zvb(ij)*uold
           xmp=xmp+xmb(ij)*uold
           ymp=ymp+ymb(ij)*uold
           zmp=zmp+zmb(ij)*uold
           xvp2=xvp2+xvb(ij)*unew
           yvp2=yvp2+yvb(ij)*unew
           zvp2=zvp2+zvb(ij)*unew
           xmp2=xmp2+xmb(ij)*unew
           ymp2=ymp2+ymb(ij)*unew
           zmp2=zmp2+zmb(ij)*unew
           l=(io-1)*nvb+(iv-nob)
           q1(1:ncent)=q1(1:ncent)+qiab(1:ncent,l)*real(uold)
         enddo
!------
! print exciton output
         xvu=xvp+xvp2
         yvu=yvp+yvp2
         zvu=zvp+zvp2
         xmu=xmp+xmp2
         ymu=ymp+ymp2
         zmu=zmp+zmp2
         xms=(yvu*coc(3)-zvu*coc(2))
         yms=(zvu*coc(1)-xvu*coc(3))
         zms=(xvu*coc(2)-yvu*coc(1))
         write(27,'(i5,a,x,F10.4,x,a2)') i,':',de*27.21139,'eV'
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'l:',xlp,ylp,zlp
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'p:',xvu,yvu,zvu
         write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'m:',xmu-xms,ymu-yms,
     .                                                zmu-zms
         do k=1,ncent
           write(27,'(f12.8)') q1(k)
         enddo
!------

! A+B transformed stuff
         ! fL stays the same
           flx=xlp*xlp
           fly=ylp*ylp
           flz=zlp*zlp
           xp=flx+fly+flz
           flp=xp*p23*de
         ! fV stays the same
           fvx=xvp*xvp
           fvy=yvp*yvp
           fvz=zvp*zvp
           xp=fvx+fvy+fvz
           fvp=p23*ef*xp
! RL is the same
           rlx=xlp*xmp
           rly=ylp*ymp
           rlz=zlp*zmp
           xp=rlx+rly+rlz
           rlp=-235.7220d0*xp
! Rv gets (symmetric) corrections from B (*vp2 and *mp2 terms)
           rvx=xvp*xmp+xvp2*xmp+xvp*xmp2
           rvy=yvp*ymp+yvp2*ymp+yvp*ymp2
           rvz=zvp*zmp+zvp2*zmp+zvp*zmp2
           xp=rvx+rvy+rvz
           rvp=-235.7220d0*xp*efmod
! print out transformed stuff
           write(28,'(i4,F10.4,4f13.6)')i,de*27.21139,flp,fvp,rlp,rvp
           rvpout(i)=rvp
           if(aniso)then
             hilf=de*2.0d0
             flx=flx*hilf
             fly=fly*hilf
             flz=flz*hilf
             hilf=ef*2.0d0
             fvx=fvx*hilf
             fvy=fvy*hilf
             fvz=fvz*hilf
             hilf=-235.7220d0
             rlx=rlx*hilf
             rly=rly*hilf
             rlz=rlz*hilf
             hilf=-235.7220d0*efmod
             rvx=rvx*hilf
             rvy=rvy*hilf
             rvz=rvz*hilf
             write(25,'(i4,F10.4,4f13.6)')i,de*27.21139,flx,fvx,rlx,rvx
             write(26,'(i4,F10.4,4f13.6)')i,de*27.21139,fly,fvy,rly,rvy
             write(29,'(i4,F10.4,4f13.6)')i,de*27.21139,flz,fvz,rlz,rvz
           endif
      enddo
      close(28)
      if(aniso)then
        close(25)
        close(26)
        close(29)
        write(*,*)'data given such that f = (f_x + f_y + f_z)/3'
        write(*,*)'                 and R = R_x + R_y + R_z'
      endif
      deallocate(xnew,q1)
      end subroutine apbtrafoexc_uks
