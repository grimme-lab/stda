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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c sRPA routine                                        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c amb: A - B (packed form)                            c
c apb: A + B (packed form)                            c
c ambsqr: (A - B)**0.5 (packed form                   c
c omsq: omega**2 - eigenvalue belonging to Z          c
c xpy: X + Y, (later xpy = x)                         c
c xmy: X - Y, (later xmy = y)                         c
c n: number of configurations in A and B matrices     c
c nroots: number of roots                             c
c ggavec : print z vector                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine srpapack(n,thr,ambsqr,apb,eci,xpy,xmy,nroots,ggavec)
      implicit none
      integer ierror,n,nroots,k,i,j,ij,m,nro,lin
      integer lwork,liwork,il,iu,info,nfound
c matrices
      real*4 ambsqr(n*(n+1)/2),xpy(n,nroots)
      real*4 apb(n*(n+1)/2),xmy(n,nroots)
      real*4 eci(n)
      real*4 summe,x,y,omsqi,vl,vu
      real*8 thr 

      integer, allocatable ::iwork(:),isuppz(:)
      real*4, allocatable ::u(:,:),v(:,:),w(:,:)
      real*4, allocatable ::z(:,:)
      real*4, allocatable ::work(:)
      real*4, allocatable ::e (:)

      logical ggavec
!      allocate(ambsqr(n*(n+1)/2),
      allocate(e(n),u(n,n),v(n,n),w(n,n),stat=ierror)  
      if(ierror.ne.0) stop 'allocation error (rpasolve)'


      call sblow(n,ambsqr,U) ! blow up sqrt(A - B) from vector to matrix
      call sblow(n,apb   ,V) ! blow up A + B from vector to matrix

!      call prmat4(6,U     ,n,n,'(A-B)^0.5') ! for debugging: print ambsqr
!      call prmat4(6,V     ,n,n,' A+B     ') ! for debugging: print apb


c form product: w = (a+b)*(a-b)^0.5
      write(*,*) ' calculating (A+B)*(A-B)^0.5 ...'
      call ssymm('l','l',n,n,1.e0,V,n,U,n,0.e0,W,n)

!      call prmat4(6,W,n,n,'(a+b)*(a-b)^0.5') ! for debugging: print w

c form product: v = (a-b)^0.5 * w

       call ssymm('l','l',n,n,1.e0,U,n,W,n,0.e0,V,n)       

!      call prmat4(6,V,n,n,'M') ! for debugging: print M matrix (V)

c get rid of matrices which are not needed anymore
      deallocate(u,w,stat=ierror)

c set variables for RPA diagonalization
      lwork =26*n
      liwork=10*n
      vl=0  
      vu=thr**2.0 ! set nroot threshold to Ethr^2
      allocate(z(n,n),work(lwork)
     .        ,iwork(liwork),isuppz(n),stat=ierror)  
      if(ierror.ne.0) stop 'allocation error (rpasolve)' 

      write(*,*)'calculate eigenvalues of (A-B)^0.5*(A+B)*(A-B)^0.5 ...'
      call ssyevr('V','V','U',n,v,n,vl,vu,il,iu,1.e-6,
     .            nfound,e,z,n,isuppz,
     .            work,lwork,iwork,liwork,info)
      nroots=nfound
      if(info.ne.0.or.nroots.lt.1) stop 'RPA diag failed'
      deallocate(v,work,iwork,isuppz,stat=ierror)
      if(ierror.ne.0) stop 'deallocation after RPA diag failed'
      
! testing print Z
!      call prmat4(6,z,n,n,'Z')

      ij=0
      do i=1,n
         ij=ij+i
         ambsqr(ij)=ambsqr(ij)*0.5
         apb(ij)=apb(ij)*0.5
      enddo

      do nro=1,nroots

         eci(nro)=sqrt(e(nro))

c (A-B)^0.5  * Z = X+Y ! results from conversion to Hermitian eigenvalue problem
         do m=1,n
            xpy(m,nro)=0.0
         enddo

         k=0
         do i=1,n
            do j=1,i
               k=k+1
               xpy(i,nro)=xpy(i,nro)+ambsqr(k)*z(j,nro)/sqrt(eci(nro)) ! dividing by sqrt(eci) yields correct norm 
               xpy(j,nro)=xpy(j,nro)+ambsqr(k)*z(i,nro)/sqrt(eci(nro))  
            enddo
         enddo

c (A+B)*|X+Y>  = e * (X-Y) ! first row of TD-DFT equation
         do m=1,n
            xmy(m,nro)=0.0
         enddo
         k=0
         do i=1,n
            do j=1,i
               k=k+1
               xmy(i,nro)=xmy(i,nro)+apb(k)*xpy(j,nro)  
               xmy(j,nro)=xmy(j,nro)+apb(k)*xpy(i,nro)  
            enddo
         enddo

c        write(*,'(''x+y'',10f8.4)')(xpy(i,nro),i=1,n)
c        write(*,'(''x-y'',10f8.4)')(xmy(i,nro),i=1,n)

         summe=0.0

         do i=1,n
            xmy(i,nro)=xmy(i,nro)/eci(nro)
            x=(xmy(i,nro)+xpy(i,nro))*0.5
            y=xpy(i,nro)-x
            xpy(i,nro)=x !xpy is now x
            xmy(i,nro)=y !xmy is now y
!            summe=summe+xpy(i,nro)**2-xmy(i,nro)**2
         enddo

c        write(*,'(''x  '',10f8.4)')(xpy(i,nro),i=1,n)
c        write(*,'(''y  '',10f8.4)')(xmy(i,nro),i=1,n)
c        write(*,*) 'norm',summe
c        write(*,*) 'e   ',eci(nro)             
c        write(*,'(''x+y'',10f8.4)')(xpy(i,nro),i=1,n)
! norm the vectors
!         summe=1.0/sqrt(summe)
!         do i=1,n
!            xpy(i,nro)=xpy(i,nro)*summe
!            xmy(i,nro)=xmy(i,nro)*summe
!         enddo

c         write(38)xpy
c         write(38)xmy

      enddo

c      close(36)
      write(*,*)' rpa vectors ok'

!      call prmat4(6,xpy,n,nro,'X')
!      call prmat4(6,xmy,n,nro,'Y')

! internal check for orthonormality
!       nro=nroots
!       z=0.0
!       call sgemm('T','n',n,nro,n,1.d0,xpy,n,xpy,n,0.d0,z,n)
!       call sgemm('T','N',n,nro,n,-1.0d0,xmy,n,xmy,n,1.0d0,z,n)
!       do i=1,min(12,nro)
!        write(6,'(12f10.6)') (z(j,i),j=1,min(12,nro))
!       enddo
!      write(6,*)
!      do i=max(1,nro-11),nro
!        write(6,'(12f10.6)') (z(j,i),j=max(1,nro-11),nro)
!      enddo

      if (ggavec) then
        call printvectda(ggavec,n,nroots,z,e) 
      endif

               
      deallocate(z,stat=ierror)  
c      deallocate(ambsqr,w5,z,stat=ierror)
      return
      end

c----------------------------------------------------------------------
c subroutine to take the power of a matrix (used here for (A-B)**0.5)
      subroutine smatpow(n,a)
      use omp_lib
      implicit none
      integer n,i,ierror,m,j,k,info
      real*4 a(n*(n+1)/2)
      real*4 summe

      real*4, allocatable ::c(:)
      real*4, allocatable ::e(:)
      real*4, allocatable ::w(:)
c working variables for sspevd diagonalization routine
      integer lwork,liwork,lin
      integer, allocatable ::iwork(:)

      lwork =1 + 6*n + n**2
      liwork=3 + 5*n
      allocate(iwork(liwork),stat=ierror)
      if(ierror.ne.0) stop 'allocation error (iwork in matpow)'

c      allocate(c(n*n),w(n*5),e(n),stat=ierror)  
      allocate(c(n*n),w(lwork),e(n),stat=ierror)
      if(ierror.ne.0) stop 'allocation error (matpow)' 

      
c      call shqrii(a,n,n,w,e,c) ! old routine - not used
c      call sspev('V','U',n,a,e,c,n,w,info) ! alternative LAPACK routine (not used)
c used LAPACK routine using divide-and-conquer algorithm (slightly faster than sspev)
      call sspevd('V','U',n,a,e,c,n,w,lwork,iwork,liwork,info) 
      if(e(1).lt.0) stop 'matrix power impossible'
c take square root of diagonal elements
      do i=1,n
         e(i)=sqrt(e(i))
      enddo
      a=0.0e0
c transform back from diagonal to non-diagonal form
      m=0
!$omp parallel private(i,j,m,k,summe)
!$omp do
      do i=1,n
         do j=1,i
            summe=0.0
            m=lin(i,j)
            do k=1,n
              summe=summe+c(i+(k-1)*n)*e(k)*c(j+(k-1)*n)
            enddo
            a(m)=summe
         enddo
      enddo
!$omp end do
!$omp end parallel

c      deallocate(c,w,e,stat=ierror)  
      deallocate(c,w,e,iwork,stat=ierror)  

      return
      end

      subroutine sblow(n,a,b)
c     blow up symmetric matrix to full size
      implicit none
      real*4 a(n*(n+1)/2),b(n,n)
      integer ij,i,n,j,lin
      ij=0
      do i=1,n
         do j=1,i-1
            ij=ij+1
            b(j,i)=a(ij)
            b(i,j)=a(ij)
         enddo
         ij=ij+1
         b(i,i)=a(ij)
      enddo
      return
      end

      subroutine sblow_fast(n,a,b)
c     blow up symmetric matrix to full size
      implicit none
      real*4 a(n*(n+1)/2),b(n,n)
      integer ij,i,n,j,lin
!$omp parallel private(i,j,ij)
!$omp do
      do i=1,n
         do j=1,i-1
            ij=lin(i,j)
            b(j,i)=a(ij)
            b(i,j)=a(ij)
         enddo
         ij=lin(i,i)
         b(i,i)=a(ij)
      enddo
!$omp end do
!$omp end parallel
      return
      end
      
      subroutine sUnblow_fast(n,a,b)
c     blow up symmetric matrix to full size
      implicit none
      real*4 b(n*(n+1)/2),a(n,n)
      integer ij,i,n,j,lin
!$omp parallel private(i,j,ij)
!$omp do
      do i=1,n
         do j=1,i-1
            ij=lin(i,j)
            b(ij)=a(j,i)
         enddo
         ij=lin(i,i)
         b(ij)=a(i,i)
      enddo
!$omp end do
!$omp end parallel
      return
      end      

      subroutine dblow_fast(n,a,b)
c     blow up symmetric matrix to full size
      implicit none
      real*8 a(n*(n+1)/2),b(n,n)
      integer ij,i,n,j,lin
!$omp parallel private(i,j,ij)
!$omp do
      do i=1,n
         do j=1,i-1
            ij=lin(i,j)
            b(j,i)=a(ij)
            b(i,j)=a(ij)
         enddo
         ij=lin(i,i)
         b(i,i)=a(ij)
      enddo
!$omp end do
!$omp end parallel
      return
      end
      
      subroutine dUnblow_fast(n,a,b)
c     blow up symmetric matrix to full size
      implicit none
      real*8 b(n*(n+1)/2),a(n,n)
      integer ij,i,n,j,lin
!$omp parallel private(i,j,ij)
!$omp do
      do i=1,n
         do j=1,i-1
            ij=lin(i,j)
            b(ij)=a(j,i)
         enddo
         ij=lin(i,i)
         b(ij)=a(i,i)
      enddo
!$omp end do
!$omp end parallel
      return
      end   
