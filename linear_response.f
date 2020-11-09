! This file is part of stda.
!
! Copyright (C) 2013-2020 Stefan Grimme
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

! written by Marc de Wegifosse 2017-2020


      SUBROUTINE lresp(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,no,nv)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)

      real*4 ::mu_x(nci)
      real*4 ::mu_y(nci)
      real*4 ::mu_z(nci)

      real*4 ::apb(nci*(nci+1)/2)
      real*4 ::amb(nci*(nci+1)/2)
      real*4, allocatable ::inv_amb(:)
      real*4, allocatable ::inv_resp(:)
      real*4, allocatable ::XpY(:,:)
      real*4 ::omega
      real*4 ::freq(num_freq+1)
      real*4 ::alpha_xx,alpha_xy,alpha_xz
      real*4 ::alpha_yy,alpha_yz
      real*4 ::alpha_zz
      character*1 ::uplo
      integer ::info
      integer, allocatable ::ipiv(:)
      real*4, allocatable ::work (:)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      open(unit=101,file='wavelength',form='formatted',status='old')
      freq(1)=0.0
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      write(*,*) 'Wavelengths (nm)'
      write(*,*) '  infinity'
      Do i=1, num_freq
      read(101,*)freq(i+1)
      write(*,*)freq(i+1)
      freq(i+1)=45.56335/freq(i+1)
      enddo
      close(101)

      allocate(inv_amb(nci*(nci+1)/2))
      inv_amb=amb
      uplo='U'
      allocate(ipiv(1:nci),work(1:nci))
      call ssptrf(uplo,nci,inv_amb,ipiv,info)
      call ssptri(uplo,nci,inv_amb,ipiv,work,info)
      deallocate(ipiv,work)

      allocate( XpY(nci,3))
      ! the dipole moment matrix mu_ai
      mu_x=0.0
      mu_y=0.0
      mu_z=0.0
!$omp parallel private(j,io,iv,idum1,idum2,ij)
!$omp do
              Do j=1, nci
                  io=iconf(j,1)
                  iv=iconf(j,2)
                  idum1=max(io,iv)
                  idum2=min(io,iv)
                  ij=idum2+idum1*(idum1-1)/2
            mu_x(j)=-2.0*xl(ij)
            mu_y(j)=-2.0*yl(ij)
            mu_z(j)=-2.0*zl(ij)
              enddo
!$omp end do
!$omp end parallel

      Do ii=1, num_freq+1
      omega=freq(ii)
      XpY(:,:)=0.0
      call cpu_time(start_time)
      allocate(inv_resp(nci*(nci+1)/2))
      inv_resp=apb-omega**2.0*inv_amb

      uplo='U'
      XpY(:,1)=mu_x(:)
      XpY(:,2)=mu_y(:)
      XpY(:,3)=mu_z(:)
      allocate(ipiv(1:nci))
      call ssptrf(uplo,nci,inv_resp,ipiv,info)
      call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY,nci,info)
      deallocate(ipiv)

      alpha_xx=0.0
      alpha_xy=0.0
      alpha_xz=0.0
      alpha_yy=0.0
      alpha_yz=0.0
      alpha_zz=0.0
!$omp parallel private(i,io,iv,idum1,idum2,ij)
!$omp&                 reduction(+:alpha_xx,alpha_xy
!$omp&                 ,alpha_xz,alpha_yy,alpha_yz
!$omp&                 ,alpha_zz)
!$omp do
      Do j=1, nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      alpha_xx=alpha_xx-(xl(ij)*XpY(j,1))*2.0
      alpha_xy=alpha_xy-(xl(ij)*XpY(j,2))*2.0
      alpha_xz=alpha_xz-(xl(ij)*XpY(j,3))*2.0
      alpha_yy=alpha_yy-(yl(ij)*XpY(j,2))*2.0
      alpha_yz=alpha_yz-(yl(ij)*XpY(j,3))*2.0
      alpha_zz=alpha_zz-(zl(ij)*XpY(j,3))*2.0
      enddo
!$omp end do
!$omp end parallel
      write(*,*)
      write(*,3333) 'Polarizability (',-45.56335/omega,
     .                             ';',45.56335/omega,')'
      write(*,*)
      write(*,1111)'x','y','z'
      write(*,2222)'x',alpha_xx,alpha_xy,alpha_xz
      write(*,2222)'y',alpha_xy,alpha_yy,alpha_yz
      write(*,2222)'z',alpha_xz,alpha_yz,alpha_zz
      write(*,*)
      write(*,111) 'Mean of alpha',(alpha_xx+alpha_yy+alpha_zz)/3.0
      write(*,*)
      deallocate(inv_resp)
      enddo
      deallocate(XpY)
      deallocate(inv_amb)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3333  format(A16,F7.1,A1,F7.1,A1)
      end subroutine lresp



      SUBROUTINE lresp1(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,no,nv)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)

      real*4 ::mu_x(nci)
      real*4 ::mu_y(nci)
      real*4 ::mu_z(nci)
      real*4 ::XpY_int(nci,3)

      real*4 ::apb(nci*(nci+1)/2)
      real*4 ::amb(nci*(nci+1)/2)
      real*4, allocatable ::inv_amb(:)
      real*4, allocatable ::inv_resp(:)
      real*8, allocatable ::XpY(:,:,:,:)
      real*4 ::omega
      real*4 ::freq(num_freq+1)
      real*8 ::alpha_xx,alpha_xy,alpha_xz
      real*8 ::alpha_yy,alpha_yz
      real*8 ::alpha_zz
      real*8, allocatable ::XmY(:,:,:,:)
      real*8, allocatable ::X(:,:,:,:)
      real*8, allocatable ::Y(:,:,:,:)
      character*1 ::uplo
      integer ::info
      integer, allocatable ::ipiv(:)
      real*4, allocatable ::work (:)

      integer ::ix,iy,iz
      real*8 ::beta(3,3,3),A,B
      real*8 ::beta2_ZZZ, beta2_XZZ, beta_HRS, DR
      real*8 ::betaVEC_X, betaVEC_Y, betaVEC_Z

      real*4 ::start_time,end_time

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)
      open(unit=555,file='beta_HRS',form='formatted',status='replace')
      open(unit=556,file='beta_tensor',form='formatted',
     .                                           status='replace')
      open(unit=101,file='wavelength',form='formatted',status='old')
      freq(1)=0.0
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      write(*,*) 'Wavelengths (nm)'
      write(*,*) '  infinity'
      Do i=1, num_freq
      read(101,*)freq(i+1)
      write(*,*)freq(i+1)
      freq(i+1)=45.56335/freq(i+1)
      enddo
      close(101)

      allocate(inv_amb(nci*(nci+1)/2))
      inv_amb=amb
      uplo='U'
      allocate(ipiv(1:nci),work(1:nci))
      call ssptrf(uplo,nci,inv_amb,ipiv,info)
      call ssptri(uplo,nci,inv_amb,ipiv,work,info)
      deallocate(ipiv,work)

      allocate( XpY(nci,num_freq+1,2,3))
      allocate( XmY(nci,num_freq+1,2,3))
      allocate( X(nci,num_freq+1,2,3),
     .          Y(nci,num_freq+1,2,3))
      XpY(:,:,:,:)=0.0
      XmY(:,:,:,:)=0.0
      X(:,:,:,:)=0.0
      Y(:,:,:,:)=0.0

      ! the dipole moment matrix mu_ai
      mu_x=0.0
      mu_y=0.0
      mu_z=0.0
!$omp parallel private(j,io,iv,idum1,idum2,ij)
!$omp do
              Do j=1, nci
                  io=iconf(j,1)
                  iv=iconf(j,2)
                  idum1=max(io,iv)
                  idum2=min(io,iv)
                  ij=idum2+idum1*(idum1-1)/2
            mu_x(j)=-2.0*xl(ij)
            mu_y(j)=-2.0*yl(ij)
            mu_z(j)=-2.0*zl(ij)
              enddo
!$omp end do
!$omp end parallel

      Do ii=1, num_freq+1
      Do jj=1, 2 !2*omega
      omega=freq(ii)*jj
      if(jj==2) omega=-omega ! from GAMESS implementation, this is what seems to be!
      call cpu_time(start_time)
      allocate(inv_resp(nci*(nci+1)/2))
      inv_resp=apb-omega**2.0*inv_amb


      uplo='U'
      XpY_int(:,1)=mu_x(:)
      XpY_int(:,2)=mu_y(:)
      XpY_int(:,3)=mu_z(:)
      allocate(ipiv(1:nci))
      call ssptrf(uplo,nci,inv_resp,ipiv,info)
      call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY_int,nci,info)
      deallocate(ipiv)

      XpY(:,ii,jj,1)=dble(XpY_int(:,1))
      XpY(:,ii,jj,2)=dble(XpY_int(:,2))
      XpY(:,ii,jj,3)=dble(XpY_int(:,3))

      if(jj==1)then !print alpha only for omega and not -2*omega
      alpha_xx=0.0
      alpha_xy=0.0
      alpha_xz=0.0
      alpha_yy=0.0
      alpha_yz=0.0
      alpha_zz=0.0
!$omp parallel private(i,io,iv,idum1,idum2,ij)
!$omp&                 reduction(+:alpha_xx,alpha_xy
!$omp&                 ,alpha_xz,alpha_yy,alpha_yz
!$omp&                 ,alpha_zz)
!$omp do
      Do j=1, nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      alpha_xx=alpha_xx-(xl(ij)*XpY(j,ii,jj,1))*2.0
      alpha_xy=alpha_xy-(xl(ij)*XpY(j,ii,jj,2))*2.0
      alpha_xz=alpha_xz-(xl(ij)*XpY(j,ii,jj,3))*2.0
      alpha_yy=alpha_yy-(yl(ij)*XpY(j,ii,jj,2))*2.0
      alpha_yz=alpha_yz-(yl(ij)*XpY(j,ii,jj,3))*2.0
      alpha_zz=alpha_zz-(zl(ij)*XpY(j,ii,jj,3))*2.0
      enddo
!$omp end do
!$omp end parallel
      write(*,*)
      write(*,3333) 'Polarizability (',-45.56335/omega,
     .                             ';',45.56335/omega,')'
      write(*,*)
      write(*,1111)'x','y','z'
      write(*,2222)'x',alpha_xx,alpha_xy,alpha_xz
      write(*,2222)'y',alpha_xy,alpha_yy,alpha_yz
      write(*,2222)'z',alpha_xz,alpha_yz,alpha_zz
      write(*,*)
      write(*,111) 'Mean of alpha',(alpha_xx+alpha_yy+alpha_zz)/3.0
      write(*,*)
      endif



c frequency-dependent first hyperpolarizability (beta) SHG only

      ! extract X and Y from XpY
      ! (X-Y)=omega*(A-B)^-1 (X+Y)
      ! X=((X+Y)+(X-Y))/2
      ! Y=(X+Y)-X
!$omp parallel private(i,j,ij,ix)
!$omp&                 reduction (+:XmY)
!$omp do
      Do i=1,nci
        Do j=1,nci
        ij=lin8(i,j)
        Do ix=1,3
        XmY(i,ii,jj,ix)=XmY(i,ii,jj,ix)+ dble(omega)
     .               *dble(inv_amb(ij))*XpY(j,ii,jj,ix)
        enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(ix,j)
!$omp&         reduction(+:X,Y)
!$omp do
      Do ix=1,3
      Do j=1,nci
      X(j,ii,jj,ix)=(XpY(j,ii,jj,ix)+XmY(j,ii,jj,ix))/2.0
      Y(j,ii,jj,ix)=XpY(j,ii,jj,ix)-X(j,ii,jj,ix)
      enddo
      enddo
!$omp end do
!$omp end parallel

      deallocate(inv_resp)

      call cpu_time(end_time)
      print '("alpha Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      enddo   !jjdo

c beta  beta = A - B + C (here C=0 since no derivative of g^XC)
      call cpu_time(start_time)

      if(ii==1)then
      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,no,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      endif
      !
      !Fast version
      !
      beta(:,:,:)=0.0
      ! Let's save time applying symmetry conditions
      if(ii==1)then !static case Kleinman condition applied ijk = kij = jki
      ! only 10 tensor components computed
      !Diagonal
      call beta_resp_fast(1,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,1,1)=A-B
      call beta_resp_fast(2,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,2,2)=A-B
      call beta_resp_fast(3,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,3,3)=A-B
      !Off-diagonal
      !ijj
      call beta_resp_fast(1,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,2,2)=A-B
      beta(2,1,2)=beta(1,2,2)
      beta(2,2,1)=beta(1,2,2)
      call beta_resp_fast(1,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,3,3)=A-B
      beta(3,1,3)=beta(1,3,3)
      beta(3,3,1)=beta(1,3,3)
      call beta_resp_fast(2,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,1,1)=A-B
      beta(1,2,1)=beta(2,1,1)
      beta(1,1,2)=beta(2,1,1)
      call beta_resp_fast(2,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,3,3)=A-B
      beta(3,2,3)=beta(2,3,3)
      beta(3,3,2)=beta(2,3,3)
      call beta_resp_fast(3,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,2,2)=A-B
      beta(2,3,2)=beta(3,2,2)
      beta(2,2,3)=beta(3,2,2)
      call beta_resp_fast(3,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,1,1)=A-B
      beta(1,3,1)=beta(3,1,1)
      beta(1,1,3)=beta(3,1,1)
      !ijk
      call beta_resp_fast(1,2,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,2,3)=A-B
      beta(1,3,2)=beta(1,2,3)
      beta(3,1,2)=beta(1,2,3)
      beta(3,2,1)=beta(1,2,3)
      beta(2,3,1)=beta(1,2,3)
      beta(2,1,3)=beta(1,2,3)

      else ! ijk=ikj because SHG (2 identical frequencies)
      !Diagonal
      call beta_resp_fast(1,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,1,1)=A-B
      call beta_resp_fast(2,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,2,2)=A-B
      call beta_resp_fast(3,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,3,3)=A-B
      !Off-diagonal
      !ijj
      call beta_resp_fast(1,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,2,2)=A-B
      call beta_resp_fast(2,1,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,1,2)=A-B
      beta(2,2,1)=beta(2,1,2)

      call beta_resp_fast(1,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,3,3)=A-B
      call beta_resp_fast(3,1,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,1,3)=A-B
      beta(3,3,1)=beta(3,1,3)
      call beta_resp_fast(2,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,1,1)=A-B
      call beta_resp_fast(1,2,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,2,1)=A-B
      beta(1,1,2)=beta(1,2,1)
      call beta_resp_fast(2,3,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,3,3)=A-B
      call beta_resp_fast(3,2,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,2,3)=A-B
      beta(3,3,2)=beta(3,2,3)
      call beta_resp_fast(3,2,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,2,2)=A-B
      call beta_resp_fast(2,3,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,3,2)=A-B
      beta(2,2,3)=beta(2,3,2)
      call beta_resp_fast(3,1,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,1,1)=A-B
      call beta_resp_fast(1,3,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,3,1)=A-B
      beta(1,1,3)=beta(1,3,1)
      !ijk
      call beta_resp_fast(1,2,3,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(1,2,3)=A-B
      beta(1,3,2)=beta(1,2,3)
      call beta_resp_fast(3,1,2,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(3,1,2)=A-B
      beta(3,2,1)=beta(3,1,2)
      call beta_resp_fast(2,3,1,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)

      beta(2,3,1)=A-B
      beta(2,1,3)=beta(2,3,1)

      endif

      call cpu_time(end_time)
      print '("beta  Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

      betaVec_X=0.0
      betaVec_Y=0.0
      betaVec_Z=0.0
      Do k=1,3
      betaVec_X=betaVec_X+beta(1,k,k)+beta(k,1,k)+beta(k,k,1)
      betaVec_Y=betaVec_Y+beta(2,k,k)+beta(k,2,k)+beta(k,k,2)
      betaVec_Z=betaVec_Z+beta(3,k,k)+beta(k,3,k)+beta(k,k,3)
      enddo
      betaVec_X=betaVec_X/5.0
      betaVec_Y=betaVec_Y/5.0
      betaVec_Z=betaVec_Z/5.0


      write(*,*) '_____________________________________________
     .________'
      write(*,*)
      write(*,5555) 'SHG first hyperpolarizability ('
     .,-45.56335/(freq(ii)*2.0),';',45.56335/freq(ii),','
     .,45.56335/freq(ii),')'
      write(*,*)
      write(556,5555) 'SHG first hyperpolarizability ('
     .,-45.56335/(freq(ii)*2.0),';',45.56335/freq(ii),','
     .,45.56335/freq(ii),')'
      write(*,*)
      call PrintBeta(beta,556)
      write(*,*)
      write(*,6666) 'beta_VEC_X',betaVec_X
      write(*,6666) 'beta_VEC_Y',betaVec_Y
      write(*,6666) 'beta_VEC_Z',betaVec_Z
      write(*,*)
      call HRS(beta,beta2_ZZZ, beta2_XZZ, beta_HRS, DR)
      write(*,9007) 6.0*beta2_ZZZ - 9.0*beta2_XZZ,
     .              sqrt(6.0*beta2_ZZZ - 9.0*beta2_XZZ),
     .              63.0/2.0*beta2_XZZ - 7.0/2.0*beta2_ZZZ,
     .              sqrt(63.0/2.0*beta2_XZZ - 7.0/2.0*beta2_ZZZ),
     .              sqrt(63.0/2.0*beta2_XZZ - 7.0/2.0*beta2_ZZZ)/
     .              sqrt(6.0*beta2_ZZZ - 9.0*beta2_XZZ)
      write(*,4444) 'beta**2_ZZZ',beta2_ZZZ, 'beta**2_XZZ', beta2_XZZ
      write(*,4444) 'beta_HRS',beta_HRS,' DR', DR
      write(*,*) '_____________________________________________
     .________'

      write(555,*)freq(ii)*27.21139,beta_HRS

!       beta(:,:,:)=0.0
!       ! Let's save time applying symmetry conditions
!       if(ii==1)then !static case Kleinman condition applied ijk = kij = jki
!       ! only 10 tensor components computed
!       !Diagonal
!       call beta_resp(1,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,1,1)=A-B
!       call beta_resp(2,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,2,2)=A-B
!       call beta_resp(3,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,3,3)=A-B
!       !Off-diagonal
!       !ijj
!       call beta_resp(1,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,2,2)=A-B
!       beta(2,1,2)=beta(1,2,2)
!       beta(2,2,1)=beta(1,2,2)
!       call beta_resp(1,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,3,3)=A-B
!       beta(3,1,3)=beta(1,3,3)
!       beta(3,3,1)=beta(1,3,3)
!       call beta_resp(2,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,1,1)=A-B
!       beta(1,2,1)=beta(2,1,1)
!       beta(1,1,2)=beta(2,1,1)
!       call beta_resp(2,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,3,3)=A-B
!       beta(3,2,3)=beta(2,3,3)
!       beta(3,3,2)=beta(2,3,3)
!       call beta_resp(3,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,2,2)=A-B
!       beta(2,3,2)=beta(3,2,2)
!       beta(2,2,3)=beta(3,2,2)
!       call beta_resp(3,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,1,1)=A-B
!       beta(1,3,1)=beta(3,1,1)
!       beta(1,1,3)=beta(3,1,1)
!       !ijk
!       call beta_resp(1,2,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,2,3)=A-B
!       beta(1,3,2)=beta(1,2,3)
!       beta(3,1,2)=beta(1,2,3)
!       beta(3,2,1)=beta(1,2,3)
!       beta(2,3,1)=beta(1,2,3)
!       beta(2,1,3)=beta(1,2,3)
!
!       else ! ijk=ikj because SHG (2 identical frequencies)
!       !Diagonal
!       call beta_resp(1,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,1,1)=A-B
!       call beta_resp(2,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,2,2)=A-B
!       call beta_resp(3,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,3,3)=A-B
!       !Off-diagonal
!       !ijj
!       call beta_resp(1,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,2,2)=A-B
!       call beta_resp(2,1,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,1,2)=A-B
!       beta(2,2,1)=beta(2,1,2)
!
!       call beta_resp(1,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,3,3)=A-B
!       call beta_resp(3,1,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,1,3)=A-B
!       beta(3,3,1)=beta(3,1,3)
!       call beta_resp(2,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,1,1)=A-B
!       call beta_resp(1,2,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,2,1)=A-B
!       beta(1,1,2)=beta(1,2,1)
!       call beta_resp(2,3,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,3,3)=A-B
!       call beta_resp(3,2,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,2,3)=A-B
!       beta(3,3,2)=beta(3,2,3)
!       call beta_resp(3,2,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,2,2)=A-B
!       call beta_resp(2,3,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,3,2)=A-B
!       beta(2,2,3)=beta(2,3,2)
!       call beta_resp(3,1,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,1,1)=A-B
!       call beta_resp(1,3,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,3,1)=A-B
!       beta(1,1,3)=beta(1,3,1)
!       !ijk
!       call beta_resp(1,2,3,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(1,2,3)=A-B
!       beta(1,3,2)=beta(1,2,3)
!       call beta_resp(3,1,2,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(3,1,2)=A-B
!       beta(3,2,1)=beta(3,1,2)
!       call beta_resp(2,3,1,X,Y,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!
!       beta(2,3,1)=A-B
!       beta(2,1,3)=beta(2,3,1)
!
!       endif
!
!       call cpu_time(end_time)
!       print '("beta  Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0
!
!       betaVec_X=0.0
!       betaVec_Y=0.0
!       betaVec_Z=0.0
!       Do k=1,3
!       betaVec_X=betaVec_X+beta(1,k,k)+beta(k,1,k)+beta(k,k,1)
!       betaVec_Y=betaVec_Y+beta(2,k,k)+beta(k,2,k)+beta(k,k,2)
!       betaVec_Z=betaVec_Z+beta(3,k,k)+beta(k,3,k)+beta(k,k,3)
!       enddo
!       betaVec_X=betaVec_X/5.0
!       betaVec_Y=betaVec_Y/5.0
!       betaVec_Z=betaVec_Z/5.0
!
!
!       write(*,*) '_____________________________________________
!      .________'
!       write(*,*)
!       write(*,5555) 'SHG first hyperpolarizability ('
!      .,-45.56335/(freq(ii)*2.0),';',45.56335/freq(ii),','
!      .,45.56335/freq(ii),')'
!       write(*,*)
!       call PrintBeta(beta)
!       write(*,*)
!       write(*,6666) 'beta_VEC_X',betaVec_X
!       write(*,6666) 'beta_VEC_Y',betaVec_Y
!       write(*,6666) 'beta_VEC_Z',betaVec_Z
!       write(*,*)
!       call HRS(beta,beta2_ZZZ, beta2_XZZ, beta_HRS, DR)
!       write(*,4444) 'beta**2_ZZZ',beta2_ZZZ, 'beta**2_XZZ', beta2_XZZ
!       write(*,4444) 'beta_HRS',beta_HRS,' DR', DR
!       write(*,*) '_____________________________________________
!      .________'

      enddo
      deallocate(XpY)
      deallocate(XmY)
      deallocate(X,Y)
      deallocate(inv_amb)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3333  format(A16,F7.1,A1,F7.1,A1)
4444  format(A12,F20.3,A12,F20.3)
5555  format(A31,F7.1,A1,F7.1,A1,F7.1,A1)
6666  format(A12,F20.3)
9007  format(2X,'|Bj=1|^2 :', F15.4,' |Bj=1| :', F15.4,
     .     /,2X,'|Bj=3|^2 :', F15.4,' |Bj=3| :', F15.4,
     .     /,2X,'rho=|Bj=3|/|Bj=1| :', F15.4,/)
      end subroutine lresp1


      SUBROUTINE beta_resp_fast(ix,iy,iz,X,Y,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::xx,yy,zz
      integer ::ix,iy,iz,wx,wy,wz,no,nv,nci,ii,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      real*8 ::A,B
      real*8 ::A1,A2,A3,A4,A5,A6
      real*8 ::B1,B2,B3,B4,B5,B6
      integer ::counter_A,counter_B
      integer :: A_list(1:counter_A,1:3)
      integer :: B_list(1:counter_B,1:3)


      A=0.0
      B=0.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      B1=0.0
      B2=0.0
      B3=0.0
      B4=0.0
      B5=0.0
      B6=0.0

c A1 ix iy iz
      xx=ix
      wx=2
      yy=iy
      wy=1
      zz=iz
      wz=1
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A1,A_list,counter_A)

c A2 ix iz iy
      xx=ix
      wx=2
      yy=iz
      wy=1
      zz=iy
      wz=1
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A2,A_list,counter_A)

c A3 iy iz ix
      xx=iy
      wx=1
      yy=iz
      wy=1
      zz=ix
      wz=2
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A3,A_list,counter_A)

c A4 iy ix iz
      xx=iy
      wx=1
      yy=ix
      wy=2
      zz=iz
      wz=1
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A4,A_list,counter_A)

c A5 iz ix iy
      xx=iz
      wx=1
      yy=ix
      wy=2
      zz=iy
      wz=1
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A5,A_list,counter_A)

c A6 iz iy ix
      xx=iz
      wx=1
      yy=iy
      wy=1
      zz=ix
      wz=2
      call A_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,A6,A_list,counter_A)

      A=A1+A2+A3+A4+A5+A6

      !write(*,*) '______________'
      !write(*,*) 'A',ix,iy,iz
      !write(*,*) A1,A2,A3
      !write(*,*) A4,A5,A6
      !write(*,*) A

c B1 ix iy iz
      xx=ix
      wx=2
      yy=iy
      wy=1
      zz=iz
      wz=1
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B1,B_list,counter_B)
c B2 ix iz iy
      xx=ix
      wx=2
      yy=iz
      wy=1
      zz=iy
      wz=1
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B2,B_list,counter_B)


c B3 iy iz ix
      xx=iy
      wx=1
      yy=iz
      wy=1
      zz=ix
      wz=2
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B3,B_list,counter_B)


c B4 iy ix iz
      xx=iy
      wx=1
      yy=ix
      wy=2
      zz=iz
      wz=1
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B4,B_list,counter_B)


c B5 iz ix iy
      xx=iz
      wx=1
      yy=ix
      wy=2
      zz=iy
      wz=1
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B5,B_list,counter_B)


c B6 iz iy ix
      xx=iz
      wx=1
      yy=iy
      wy=1
      zz=ix
      wz=2
      call B_beta_fast(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ii,B6,B_list,counter_B)

      B=B1+B2+B3+B4+B5+B6


      !write(*,*) 'B',ix,iy,iz
      !write(*,*) B1,B2,B3
      !write(*,*) B4,B5,B6
      !write(*,*) B
      !write(*,*) '______________'

      end subroutine beta_resp_fast

      Subroutine List_A(maxconf,no,nci,iconf,A_list,counter_A)
      use commonresp
      use omp_lib
      implicit none

      integer ::no,nci,maxconf
      integer ::iconf(maxconf,2)
      integer ::counter_A,counter
      integer ::A_list(1:counter_A,1:3)
      integer ::i,j,ij,kk,io1,io2,idum1,idum2
      integer ::jwrk

      counter=1
!$omp parallel private(j,i,kk,jwrk,io1,io2,idum1,idum2,ij)
!$omp do ordered
      Do i=1,nci
      !$OMP ORDERED
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then

            jwrk=kk
            A_list(counter,1)=i
            io1=iconf(i,1)
            io2=j
            idum1=max(io1,io2)
            idum2=min(io1,io2)
            ij=idum2+idum1*(idum1-1)/2
            A_list(counter,2)=ij
            A_list(counter,3)=jwrk
            counter=counter+1

        endif
        enddo
      enddo
      !$OMP END ORDERED
      enddo
!$omp end do
!$omp end parallel
      end subroutine List_A

      Subroutine List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      use commonresp
      use omp_lib
      implicit none

      integer ::no,nv,nci,maxconf
      integer ::iconf(maxconf,2)
      integer ::counter_B,counter
      integer ::B_list(1:counter_B,1:3)
      integer ::i,j,ab,kk,iv1,iv2,idum1,idum2
      integer ::jwrk

      counter=1
!$omp parallel private(j,i,kk,jwrk,iv1,iv2,idum1,idum2,ab)
!$omp do ordered
      Do i=1,nci
      !$OMP ORDERED
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            jwrk=kk
            B_list(counter,1)=i
            iv1=iconf(i,2)
            iv2=j+no
            idum1=max(iv1,iv2)
            idum2=min(iv1,iv2)
            ab=idum2+idum1*(idum1-1)/2
            B_list(counter,2)=ab
            B_list(counter,3)=jwrk
            counter=counter+1
        endif
        enddo
      enddo
      !$OMP END ORDERED
      enddo
!$omp end do
!$omp end parallel
      end subroutine List_B


      Subroutine A_beta_fast(ix,iy,iz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ifreq,A,A_list,counter_A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,iz,wx,wy,wz,nci,ifreq,moci
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      integer ::i,ii

      integer ::counter_A
      integer ::A_list(1:counter_A,1:3)

      ii=ifreq
      A=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+X(A_list(i,1),ii,wx,ix)*(-mu(A_list(i,2),iy))
     .                         *Y(A_list(i,3),ii,wz,iz)
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_beta_fast

      Subroutine B_beta_fast(ix,iy,iz,wx,wy,wz,X,Y,
     .           mu,nci,moci,ifreq,B,B_list,counter_B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,iz,wx,wy,wz,nci,ifreq,moci
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      integer ::i,ii

      integer ::counter_B
      integer ::B_list(1:counter_B,1:3)

      ii=ifreq
      B=0.0
!$omp parallel private(i)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+X(B_list(i,1),ii,wx,ix)*(-mu(B_list(i,2),iy))
     .                         *Y(B_list(i,3),ii,wz,iz)
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_beta_fast

      SUBROUTINE beta_resp(ix,iy,iz,X,Y,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::xx,yy,zz
      integer ::ix,iy,iz,wx,wy,wz,no,nv,nci,ii,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      real*8 ::A,B
      real*8 ::A1,A2,A3,A4,A5,A6
      real*8 ::B1,B2,B3,B4,B5,B6

      A=0.0
      B=0.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      B1=0.0
      B2=0.0
      B3=0.0
      B4=0.0
      B5=0.0
      B6=0.0

c A1 ix iy iz
      xx=ix
      wx=2
      yy=iy
      wy=1
      zz=iz
      wz=1
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A1)

c A2 ix iz iy
      xx=ix
      wx=2
      yy=iz
      wy=1
      zz=iy
      wz=1
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A2)

c A3 iy iz ix
      xx=iy
      wx=1
      yy=iz
      wy=1
      zz=ix
      wz=2
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A3)
c A4 iy ix iz
      xx=iy
      wx=1
      yy=ix
      wy=2
      zz=iz
      wz=1
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A4)
c A5 iz ix iy
      xx=iz
      wx=1
      yy=ix
      wy=2
      zz=iy
      wz=1
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A5)
c A6 iz iy ix
      xx=iz
      wx=1
      yy=iy
      wy=1
      zz=ix
      wz=2
      call A_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A6)

      A=A1+A2+A3+A4+A5+A6

      !write(*,*) '______________'
      !write(*,*) 'A',ix,iy,iz
      !write(*,*) A1,A2,A3
      !write(*,*) A4,A5,A6
      !write(*,*) A

c B1 ix iy iz
      xx=ix
      wx=2
      yy=iy
      wy=1
      zz=iz
      wz=1
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B1)

c B2 ix iz iy
      xx=ix
      wx=2
      yy=iz
      wy=1
      zz=iy
      wz=1
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B2)


c B3 iy iz ix
      xx=iy
      wx=1
      yy=iz
      wy=1
      zz=ix
      wz=2
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B3)


c B4 iy ix iz
      xx=iy
      wx=1
      yy=ix
      wy=2
      zz=iz
      wz=1
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B4)


c B5 iz ix iy
      xx=iz
      wx=1
      yy=ix
      wy=2
      zz=iy
      wz=1
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B5)


c B6 iz iy ix
      xx=iz
      wx=1
      yy=iy
      wy=1
      zz=ix
      wz=2
      call B_beta1(xx,yy,zz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B6)

      B=B1+B2+B3+B4+B5+B6


      !write(*,*) 'B',ix,iy,iz
      !write(*,*) B1,B2,B3
      !write(*,*) B4,B5,B6
      !write(*,*) B
      !write(*,*) '______________'

      end subroutine beta_resp

      Subroutine A_beta1(ix,iy,iz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,iz,wx,wy,wz,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      integer ::i,j,k,ij,ii,kk,io1,io2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      A=0.0
!$omp parallel private(j,i,kk,jwrk,io1,io2,idum1,idum2,ij)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,nci
      Do j=1,no
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      io1=iconf(i,1)
      io2=j
      idum1=max(io1,io2)
      idum2=min(io1,io2)
      ij=idum2+idum1*(idum1-1)/2
      A=A+X(i,ii,wx,ix)*(-mu(ij,iy))*Y(jwrk,ii,wz,iz)
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_beta1

      Subroutine B_beta1(ix,iy,iz,wx,wy,wz,X,Y,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,iz,wx,wy,wz,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,num_freq+1,2,3)
      real*8 ::Y(nci,num_freq+1,2,3)
      integer ::i,j,k,ab,ii,kk,iv1,iv2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      B=0.0
!$omp parallel private(j,i,kk,jwrk,iv1,iv2,idum1,idum2,ab)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      iv1=iconf(i,2)
      iv2=j+no
      idum1=max(iv1,iv2)
      idum2=min(iv1,iv2)
      ab=idum2+idum1*(idum1-1)/2
      B=B+X(i,ii,wx,ix)*(-mu(ab,iy))*Y(jwrk,ii,wz,iz)
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_beta1

      SUBROUTINE HRS(beta,beta2_ZZZ, beta2_XZZ, beta_HRS, DR)
      implicit none
      INTEGER::i,j,k
      REAL*8::beta(1:3,1:3,1:3)
      REAL*8:: beta2_ZZZ, beta2_XZZ, beta2_iii, beta2_iij, beta2_jii,
     .         beta2_ijk, beta_iii_beta_ijj, beta_jii_beta_iij,
     .         beta_iii_beta_jji,beta_iij_beta_jkk, beta_jii_beta_jkk,
     .         beta_iij_beta_kkj, beta_ijk_beta_jik, beta2_ijj,
     .         beta_iij_beta_jii, beta_ijj_beta_ikk, beta_iik_beta_jjk,
     .         DR, beta_HRS

! Calcul du beta_HRS

      beta2_iii = beta(1,1,1)**2 + beta(2,2,2)**2 + beta(3,3,3)**2

      beta2_iij = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta2_iij = beta2_iij + beta(i,i,j)**2
      endif
      enddo
      enddo

      beta2_jii = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta2_jii = beta2_jii + beta(j,i,i)**2
      endif
      enddo
      enddo

      beta2_ijk = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta2_ijk = beta2_ijk + beta(i,j,k)**2
      endif
      enddo
      enddo
      enddo

      beta_iii_beta_ijj = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta_iii_beta_ijj = beta_iii_beta_ijj + beta(i,i,i)*beta(i,j,j)
      endif
      enddo
      enddo

      beta_jii_beta_iij = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta_jii_beta_iij = beta_jii_beta_iij + beta(j,i,i)*beta(i,i,j)
      endif
      enddo
      enddo

      beta_iii_beta_jji = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta_iii_beta_jji = beta_iii_beta_jji + beta(i,i,i)*beta(j,j,i)
      endif
      enddo
      enddo

      beta_iij_beta_jkk = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_iij_beta_jkk = beta_iij_beta_jkk + beta(i,i,j)*beta(j,k,k)
      endif
      enddo
      enddo
      enddo

      beta_jii_beta_jkk = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_jii_beta_jkk = beta_jii_beta_jkk + beta(j,i,i)*beta(j,k,k)
      endif
      enddo
      enddo
      enddo

      beta_iij_beta_kkj = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_iij_beta_kkj = beta_iij_beta_kkj + beta(i,i,j)*beta(k,k,j)
      endif
      enddo
      enddo
      enddo

      beta_ijk_beta_jik = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_ijk_beta_jik = beta_ijk_beta_jik + beta(i,j,k)*beta(j,i,k)
      endif
      enddo
      enddo
      enddo

      beta2_ijj = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta2_ijj = beta2_ijj + beta(i,j,j)**2
      endif
      enddo
      enddo

      beta_iij_beta_jii = 0.0
      Do i=1,3
      Do j=1,3
      if(i/=j)then
      beta_iij_beta_jii = beta_iij_beta_jii + beta(i,i,j)*beta(j,i,i)
      endif
      enddo
      enddo

      beta_ijj_beta_ikk = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_ijj_beta_ikk = beta_ijj_beta_ikk + beta(i,j,j)*beta(i,k,k)
      endif
      enddo
      enddo
      enddo

      beta_iik_beta_jjk = 0.0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      if(i/=j .AND. i/=k .AND. j/=k)then
      beta_iik_beta_jjk = beta_iik_beta_jjk + beta(i,i,k)*beta(j,j,k)
      endif
      enddo
      enddo
      enddo


      beta2_ZZZ=1./7.*beta2_iii + 4./35.*beta2_iij + 1./35.*beta2_jii
     .          + 2./105.*beta2_ijk + 2./35.*beta_iii_beta_ijj
     .          + 4./35.*beta_jii_beta_iij + 4./35.*beta_iii_beta_jji
     .          + 4./105.*beta_iij_beta_jkk + 1./105.*beta_jii_beta_jkk
     .          + 4./105.*beta_iij_beta_kkj + 4./105.*beta_ijk_beta_jik


      beta2_XZZ=1./35.*beta2_iii + 4./105.*beta_iii_beta_ijj
     .          - 2./35.*beta_iii_beta_jji + 8./105.*beta2_iij
     .          - 2./105.*beta_iij_beta_jkk + 2./35.*beta2_ijk
     .          - 2./105.*beta_ijk_beta_jik + 3./35.*beta2_ijj
     .          - 2./35.*beta_iij_beta_jii + 1./35.*beta_ijj_beta_ikk
     .          - 2./105.*beta_iik_beta_jjk

      DR=beta2_ZZZ/beta2_XZZ
      beta_HRS=sqrt(beta2_ZZZ+beta2_XZZ)

      end subroutine HRS

      subroutine PrintBeta(beta,unit)
      implicit none
c     Arguments
      real*8 beta(3,3,3)
c     Variables and constants
      integer i,j,k,unit
      character*1 FIELDDIR(3)
c     Body of the function
      DATA FIELDDIR /'x','y','z'/
      write(*,9000) (FIELDDIR(i), i=1,3)
      do i = 1,3
         do j = 1,3
            write(*,9001) FIELDDIR(i),FIELDDIR(j),
     &                     (beta(i,j,k), k=1,3)
         write(unit,9002)(beta(i,j,k), k=1,3)
         end do
      end do
 9000 format(6X,          3(2X,7X,A1,7X))
 9001 format(3X,A1,A1,".",3(2X,F15.6))
 9002 format(3(2X,F15.6))
      end



      SUBROUTINE lresp_2PA(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,Yci,nroot)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*4 ::mu_x(nci)
      real*4 ::mu_y(nci)
      real*4 ::mu_z(nci)
      real*4 ::XpY_int(nci,3)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::omega
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)
      real*4 ::apb(nci*(nci+1)/2)
      real*4 ::amb(nci*(nci+1)/2)
      real*4, allocatable ::inv_amb(:)
      real*4, allocatable ::inv_resp(:)
      real*8, allocatable ::XpY(:,:)
      real*8, allocatable ::XmY(:,:)
      real*8, allocatable ::X(:,:)
      real*8, allocatable ::Y(:,:)
      character*1 ::uplo
      integer ::info
      integer, allocatable ::ipiv(:)
      real*4, allocatable ::work (:)

      integer ::ix,iy,iz
      real*8 ::sigma(3,3),A,B,sigma_f,sigma_g,sigma_h

      real*8 ::alpha_xx,alpha_xy,alpha_xz
      real*8 ::alpha_yy,alpha_yz
      real*8 ::alpha_zz


      real*4 ::start_time,end_time

      open(unit=60,file='2PA-abs',status='replace')

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)

      allocate(inv_amb(nci*(nci+1)/2))
      inv_amb=amb
      uplo='U'
      allocate(ipiv(1:nci),work(1:nci))
      call ssptrf(uplo,nci,inv_amb,ipiv,info)
      call ssptri(uplo,nci,inv_amb,ipiv,work,info)
      deallocate(ipiv,work)

      allocate( XpY(nci,3))
      allocate( XmY(nci,3))
      allocate( X(nci,3),
     .          Y(nci,3))

      ! the dipole moment matrix mu_ai
      mu_x=0.0
      mu_y=0.0
      mu_z=0.0
!$omp parallel private(j,io,iv,idum1,idum2,ij)
!$omp do
             Do j=1, nci
                 io=iconf(j,1)
                 iv=iconf(j,2)
                 idum1=max(io,iv)
                 idum2=min(io,iv)
                 ij=idum2+idum1*(idum1-1)/2
           mu_x(j)=-2.0*xl(ij)
           mu_y(j)=-2.0*yl(ij)
           mu_z(j)=-2.0*zl(ij)
             enddo
!$omp end do
!$omp end parallel

      Do ii=1, num_trans
      XpY(:,:)=0.0
      XmY(:,:)=0.0
      X(:,:)=0.0
      Y(:,:)=0.0
      omega=eci(ii)/2.0
      call cpu_time(start_time)
      allocate(inv_resp(nci*(nci+1)/2))
      inv_resp=apb-omega**2.0*inv_amb


      uplo='U'
      XpY_int(:,1)=mu_x(:)
      XpY_int(:,2)=mu_y(:)
      XpY_int(:,3)=mu_z(:)
      allocate(ipiv(1:nci))
      call ssptrf(uplo,nci,inv_resp,ipiv,info)
      call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY_int,nci,info)
      deallocate(ipiv)

      XpY(:,1)=dble(XpY_int(:,1))
      XpY(:,2)=dble(XpY_int(:,2))
      XpY(:,3)=dble(XpY_int(:,3))


!       alpha_xx=0.0
!       alpha_xy=0.0
!       alpha_xz=0.0
!       alpha_yy=0.0
!       alpha_yz=0.0
!       alpha_zz=0.0
! !$omp parallel private(i,io,iv,idum1,idum2,ij)
! !$omp&                 reduction(+:alpha_xx,alpha_xy
! !$omp&                 ,alpha_xz,alpha_yy,alpha_yz
! !$omp&                 ,alpha_zz)
! !$omp do
!       Do j=1, nci
!             io=iconf(j,1)
!             iv=iconf(j,2)
!             idum1=max(io,iv)
!             idum2=min(io,iv)
!             ij=idum2+idum1*(idum1-1)/2
!       alpha_xx=alpha_xx-(xl(ij)*XpY(j,1))*2.0
!       alpha_xy=alpha_xy-(xl(ij)*XpY(j,2))*2.0
!       alpha_xz=alpha_xz-(xl(ij)*XpY(j,3))*2.0
!       alpha_yy=alpha_yy-(yl(ij)*XpY(j,2))*2.0
!       alpha_yz=alpha_yz-(yl(ij)*XpY(j,3))*2.0
!       alpha_zz=alpha_zz-(zl(ij)*XpY(j,3))*2.0
!       enddo
! !$omp end do
! !$omp end parallel
       write(*,*)
       write(*,*)ii
       write(*,*)
!       write(*,3334) 'Polarizability (',-omega*27.21139,
!      .                             ';',omega*27.21139,')'
!       write(*,*)
!       write(*,1111)'x','y','z'
!       write(*,2222)'x',alpha_xx,alpha_xy,alpha_xz
!       write(*,2222)'y',alpha_xy,alpha_yy,alpha_yz
!       write(*,2222)'z',alpha_xz,alpha_yz,alpha_zz
!       write(*,*)
!       write(*,111) 'Mean of alpha',(alpha_xx+alpha_yy+alpha_zz)/3.0
!       write(*,*)


      ! extract X and Y from XpY
      ! (X-Y)=omega*(A-B)^-1 (X+Y)
      ! X=((X+Y)+(X-Y))/2
      ! Y=(X+Y)-X
!$omp parallel private(i,j,ij,ix)
!$omp&                 reduction (+:XmY)
!$omp do
      Do i=1,nci
        Do j=1,nci
        ij=lin8(i,j)
        Do ix=1,3
        XmY(i,ix)=XmY(i,ix)+ dble(omega)
     .               *dble(inv_amb(ij))*XpY(j,ix)
        enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(ix,j)
!$omp&         reduction(+:X,Y)
!$omp do
      Do ix=1,3
      Do j=1,nci
      X(j,ix)=(XpY(j,ix)+XmY(j,ix))/2.0
      Y(j,ix)=XpY(j,ix)-X(j,ix)
      enddo
      enddo
!$omp end do
!$omp end parallel

      deallocate(inv_resp)

      call cpu_time(end_time)
      print '("alpha Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

      if(ii==1)then
      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      call cpu_time(start_time)
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,no,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      endif

c sigma  sigma = -A + B
      call cpu_time(start_time)
      !
      ! Fast version
      !
      sigma(:,:)=0.0
      Do ix=1,3
      Do iy=1,ix
      call TPA_resp_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
      sigma(ix,iy)=2.0**(-1.0/2.0)*(-A+B)
      if(ix/=iy)then
      sigma(iy,ix)=sigma(ix,iy)
      endif
      enddo
      enddo

      sigma_f=0.0
      sigma_g=0.0
      sigma_h=0.0
      Do ix=1,3
      Do iy=1,3
      sigma_f=sigma_f+sigma(ix,ix)*sigma(iy,iy)
      sigma_g=sigma_g+sigma(ix,iy)*sigma(ix,iy)
      sigma_h=sigma_h+sigma(ix,iy)*sigma(iy,ix)
      enddo
      enddo
      sigma_f=sigma_f/30.0
      sigma_g=sigma_g/30.0
      sigma_h=sigma_h/30.0



      call cpu_time(end_time)
      print '("2PA  Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      write(*,3333)'Delta (',eci(ii)*27.21139,')'
      write(*,*)
      write(*,1111)'x','y','z'
      write(*,2222)'x',sigma(1,1),sigma(1,2),sigma(1,3)
      write(*,2222)'y',sigma(2,1),sigma(2,2),sigma(2,3)
      write(*,2222)'z',sigma(3,1),sigma(3,2),sigma(3,3)
      write(*,*)
      write(*,5555)'F =',sigma_f,' G =',sigma_G,' H =',sigma_H
      write(*,4444)'Delta_2PA_//   =',2.0*sigma_f+2.0*sigma_g
     .                               +2.0*sigma_h
      write(*,4444)'Delta_2PA__|_  =',-1.0*sigma_f+4.0*sigma_g
     .                               -1.0*sigma_h
      write(*,4444)'Delta_2PA_circ =',-2.0*sigma_f+3.0*sigma_g
     .                               +3.0*sigma_h
      write(*,4444)'rho = //*(_|_)**-1 =',(2.0*sigma_f+2.0*sigma_g
     .       +2.0*sigma_h)/(-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h)

!       call cpu_time(start_time)
!       !
!       ! regular version
!       !
!       sigma(:,:)=0.0
!       Do ix=1,3
!       Do iy=1,ix
!       call TPA_resp(ix,iy,X,Y,Xci,Yci,nroot,
!      .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
!       sigma(ix,iy)=1.0/sqrt(2.0)*(-A+B)
!       if(ix/=iy)then
!       sigma(iy,ix)=sigma(ix,iy)
!       endif
!       enddo
!       enddo
!
!       sigma_f=0.0
!       sigma_g=0.0
!       sigma_h=0.0
!       Do ix=1,3
!       Do iy=1,3
!       sigma_f=sigma_f+sigma(ix,ix)*sigma(iy,iy)
!       sigma_g=sigma_g+sigma(ix,iy)*sigma(ix,iy)
!       sigma_h=sigma_h+sigma(ix,iy)*sigma(iy,ix)
!       enddo
!       enddo
!       sigma_f=sigma_f/30.0
!       sigma_g=sigma_g/30.0
!       sigma_h=sigma_h/30.0
!
!
!
!       call cpu_time(end_time)
!       print '("2PA  Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0
!       write(*,*)
!       write(*,3333)'Delta (',eci(ii)*27.21139,')'
!       write(*,*)
!       write(*,1111)'x','y','z'
!       write(*,2222)'x',sigma(1,1),sigma(1,2),sigma(1,3)
!       write(*,2222)'y',sigma(2,1),sigma(2,2),sigma(2,3)
!       write(*,2222)'z',sigma(3,1),sigma(3,2),sigma(3,3)
!       write(*,*)
!       write(*,5555)'F =',sigma_f,' G =',sigma_G,' H =',sigma_H
!       !write(*,*)'like in gamess, only F and G, and wrong
!       !.definition for circ. and rho'
!       !write(*,*)'Sigma_2PA_//   =',2.0*sigma_f+4.0*sigma_g
!       !write(*,*)'Sigma_2PA_circ =',-2.0*sigma_f+6.0*sigma_g
!       !write(*,*)'Rho =',(-sigma_f+3.0*sigma_g)/(sigma_f+2.0*sigma_g)
!       !write(*,*)'2/Rho=',2.0/((-sigma_f+3.0*sigma_g)/
!       !.(sigma_f+2.0*sigma_g))
!       !write(*,*)'with F,G, and H'
!       write(*,4444)'Delta_2PA_//   =',2.0*sigma_f+2.0*sigma_g
!      .                               +2.0*sigma_h
!       write(*,4444)'Delta_2PA__|_  =',-1.0*sigma_f+4.0*sigma_g
!      .                               -1.0*sigma_h
!       write(*,4444)'Delta_2PA_circ =',-2.0*sigma_f+3.0*sigma_g
!      .                               +3.0*sigma_h
!       write(*,4444)'rho = //*(_|_)**-1 =',(2.0*sigma_f+2.0*sigma_g
!      .       +2.0*sigma_h)/(-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h)

      write(60,6666)ii,eci(ii)*27.21139,2.0*sigma_f+2.0*sigma_g
     .+2.0*sigma_h,-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h,-2.0*sigma_f
     .+3.0*sigma_g+3.0*sigma_h,(2.0*sigma_f+2.0*sigma_g
     .+2.0*sigma_h)/(-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h)
      enddo
      close(60)
      deallocate(XpY)
      deallocate(XmY)
      deallocate(X,Y)
      deallocate(inv_amb)
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3334  format(A16,F7.3,A1,F7.3,A1)
3333  format(A16,F7.3,A1)
4444  format(A20,F20.3)
5555  format(A3,F20.3,A4,F20.3,A4,F20.3)
6666  format(I3,F7.3,4F20.3)
      end subroutine lresp_2PA

      SUBROUTINE TPA_resp_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .            A_list,B_list,counter_A,counter_B,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
      ! not stable due to the linear system to solve
      use commonresp
      use omp_lib
      implicit none

      integer ::xx,yy,nroot
      integer ::ix,iy,no,nv,nci,ii,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      real*8 ::A,B
      real*8 ::A1,A2,A3,A4
      real*8 ::B1,B2,B3,B4
      integer ::counter_A,counter_B
      integer :: A_list(1:counter_A,1:3)
      integer :: B_list(1:counter_B,1:3)

      A=0.0
      B=0.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      B1=0.0
      B2=0.0
      B3=0.0
      B4=0.0

c A1 ix iy n
      xx=ix
      yy=iy
      call A_2PA_1_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            A_list,counter_A,
     .           mu,nci,moci,ii,A1)
c A2 iy ix n
      xx=iy
      yy=ix
      call A_2PA_1_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            A_list,counter_A,
     .           mu,nci,moci,ii,A2)
c A3 n ix iy
      xx=ix
      yy=iy
      call A_2PA_2_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            A_list,counter_A,
     .           mu,nci,moci,ii,A3)
c A4 n iy ix
      xx=iy
      yy=ix
      call A_2PA_2_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            A_list,counter_A,
     .           mu,nci,moci,ii,A4)
      A=A1+A2+A3+A4

c B1 ix iy n
      xx=ix
      yy=iy
      call B_2PA_1_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            B_list,counter_B,
     .           mu,nci,moci,ii,B1)
c B2 iy ix n
      xx=iy
      yy=ix
      call B_2PA_1_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            B_list,counter_B,
     .           mu,nci,moci,ii,B2)
c B3 n ix iy
      xx=ix
      yy=iy
      call B_2PA_2_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            B_list,counter_B,
     .           mu,nci,moci,ii,B3)
c B4 n iy ix
      xx=iy
      yy=ix
      call B_2PA_2_fast(xx,yy,X,Y,Xci,Yci,nroot,
     .            B_list,counter_B,
     .           mu,nci,moci,ii,B4)
      B=B1+B2+B3+B4

      end subroutine TPA_resp_fast

      Subroutine A_2PA_1_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .           A_list,counter_A,
     .           mu,nci,moci,ifreq,A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,nci,ifreq,moci
      integer ::nroot
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,ii

      integer ::counter_A
      integer ::A_list(1:counter_A,1:3)

      ii=ifreq
      A=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+X(A_list(i,1),ix)*(-mu(A_list(i,2),iy))
     .                 *dble(Yci(A_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_2PA_1_fast

      Subroutine A_2PA_2_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .           A_list,counter_A,
     .           mu,nci,moci,ifreq,A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,nci,ifreq,moci
      integer ::nroot
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,ii

      integer ::counter_A
      integer ::A_list(1:counter_A,1:3)

      ii=ifreq
      A=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+dble(Xci(A_list(i,1),ii))*(-mu(A_list(i,2),ix))
     .                                 *Y(A_list(i,3),iy)
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_2PA_2_fast

      Subroutine B_2PA_1_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .           B_list,counter_B,
     .           mu,nci,moci,ifreq,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,nci,ifreq,moci
      integer ::nroot
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,ii

      integer ::counter_B
      integer ::B_list(1:counter_B,1:3)

      ii=ifreq
      B=0.0
!$omp parallel private(i)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+X(B_list(i,1),ix)*(-mu(B_list(i,2),iy))
     .                 *dble(Yci(B_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_2PA_1_fast

      Subroutine B_2PA_2_fast(ix,iy,X,Y,Xci,Yci,nroot,
     .           B_list,counter_B,
     .           mu,nci,moci,ifreq,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,nci,ifreq,moci
      integer ::nroot
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,ii

      integer ::counter_B
      integer ::B_list(1:counter_B,1:3)

      ii=ifreq
      B=0.0
!$omp parallel private(i)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+dble(Xci(B_list(i,1),ii))*(-mu(B_list(i,2),ix))
     .*Y(B_list(i,3),iy)
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_2PA_2_fast



      SUBROUTINE TPA_resp(ix,iy,X,Y,Xci,Yci,nroot,
     .            mu,maxconf,no,nv,nci,moci,ii,iconf,A,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::xx,yy,nroot
      integer ::ix,iy,no,nv,nci,ii,maxconf,moci
      integer ::iconf(maxconf,2)
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      real*8 ::A,B
      real*8 ::A1,A2,A3,A4
      real*8 ::B1,B2,B3,B4

      A=0.0
      B=0.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      B1=0.0
      B2=0.0
      B3=0.0
      B4=0.0

c A1 ix iy n
      xx=ix
      yy=iy
      call A_2PA_1(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A1)
c A2 iy ix n
      xx=iy
      yy=ix
      call A_2PA_1(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A2)
c A3 n ix iy
      xx=ix
      yy=iy
      call A_2PA_2(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A3)
c A4 n iy ix
      xx=iy
      yy=ix
      call A_2PA_2(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,A4)
      A=A1+A2+A3+A4

c B1 ix iy n
      xx=ix
      yy=iy
      call B_2PA_1(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B1)
c B2 iy ix n
      xx=iy
      yy=ix
      call B_2PA_1(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B2)
c B3 n ix iy
      xx=ix
      yy=iy
      call B_2PA_2(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B3)
c B4 n iy ix
      xx=iy
      yy=ix
      call B_2PA_2(xx,yy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ii,iconf,B4)
      B=B1+B2+B3+B4

      end subroutine TPA_resp




      Subroutine A_2PA_1(ix,iy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2),nroot
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,j,k,ij,ii,kk,io1,io2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      A=0.0
!$omp parallel private(j,i,kk,jwrk,io1,io2,idum1,idum2,ij)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,nci
      Do j=1,no
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      io1=iconf(i,1)
      io2=j
      idum1=max(io1,io2)
      idum2=min(io1,io2)
      ij=idum2+idum1*(idum1-1)/2
      A=A+X(i,ix)*(-mu(ij,iy))*dble(Yci(jwrk,ii))
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_2PA_1

      Subroutine A_2PA_2(ix,iy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,A)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2),nroot
      real*8  ::A
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,j,k,ij,ii,kk,io1,io2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      A=0.0
!$omp parallel private(j,i,kk,jwrk,io1,io2,idum1,idum2,ij)
!$omp&         reduction(+:A)
!$omp do
      Do i=1,nci
      Do j=1,no
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      io1=iconf(i,1)
      io2=j
      idum1=max(io1,io2)
      idum2=min(io1,io2)
      ij=idum2+idum1*(idum1-1)/2
      A=A+dble(Xci(i,ii))*(-mu(ij,ix))*Y(jwrk,iy)
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine A_2PA_2

      Subroutine B_2PA_1(ix,iy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2),nroot
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,j,k,ab,ii,kk,iv1,iv2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      B=0.0
!$omp parallel private(j,i,kk,jwrk,iv1,iv2,idum1,idum2,ab)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      iv1=iconf(i,2)
      iv2=j+no
      idum1=max(iv1,iv2)
      idum2=min(iv1,iv2)
      ab=idum2+idum1*(idum1-1)/2
      B=B+X(i,ix)*(-mu(ab,iy))*dble(Yci(jwrk,ii))
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_2PA_1

      Subroutine B_2PA_2(ix,iy,X,Y,Xci,Yci,nroot,
     .           mu,maxconf,no,nv,nci,moci,ifreq,iconf,B)
      use commonresp
      use omp_lib
      implicit none

      integer ::ix,iy,no,nv,nci,ifreq,maxconf,moci
      integer ::iconf(maxconf,2),nroot
      real*8  ::B
      real*8 ::mu(moci*(moci+1)/2,3)
      real*8 ::X(nci,3)
      real*8 ::Y(nci,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot)
      integer ::i,j,k,ab,ii,kk,iv1,iv2,idum1,idum2
      integer*8 ::lin8
      integer ::jwrk
      logical ::check

      ii=ifreq
      B=0.0
!$omp parallel private(j,i,kk,jwrk,iv1,iv2,idum1,idum2,ab)
!$omp&         reduction(+:B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        check=.false.
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            jwrk=kk
            check=.true.
            endif
        enddo
      if(check .eqv. .true.)then
      iv1=iconf(i,2)
      iv2=j+no
      idum1=max(iv1,iv2)
      idum2=min(iv1,iv2)
      ab=idum2+idum1*(idum1-1)/2
      B=B+dble(Xci(i,ii))*(-mu(ab,ix))*Y(jwrk,iy)
      endif
      enddo
      enddo
!$omp end do
!$omp end parallel
      end subroutine B_2PA_2

      subroutine pol_sos(nroot,nci,eci,Xci,Yci,xl,yl,zl,moci,
     .maxconf,iconf,ak)
      use commonresp
      use omp_lib
      implicit none

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,moci
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::omega
      real*4 ::freq(num_freq+1)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)

      real*8 ::alpha_xx,alpha_xy,alpha_xz
      real*8 ::alpha_yy,alpha_yz
      real*8 ::alpha_zz,ak
      real*8 ::x,y,z

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      open(unit=101,file='wavelength',form='formatted',status='old')
      freq(1)=0.0

      write(*,*)
      write(*,*) 'Wavelengths (nm)'
      write(*,*) '  infinity'
      Do i=1, num_freq
      read(101,*)freq(i+1)
      write(*,*)freq(i+1)
      freq(i+1)=45.56335/freq(i+1)
      enddo
      close(101)

      Do ii=1, num_freq+1
      omega=freq(ii)
      alpha_xx=0.0
      alpha_xy=0.0
      alpha_xz=0.0
      alpha_yy=0.0
      alpha_yz=0.0
      alpha_zz=0.0

      Do i=1,nroot
      x=0.0
      y=0.0
      z=0.0
      Do j=1,nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2

      x=x+mu(ij,1)*(dble(Xci(j,i))+dble(Yci(j,i)))
      y=y+mu(ij,2)*(dble(Xci(j,i))+dble(Yci(j,i)))
      z=z+mu(ij,3)*(dble(Xci(j,i))+dble(Yci(j,i)))
      enddo

      alpha_xx=alpha_xx-(ak*x*x)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      alpha_yy=alpha_yy-(ak*y*y)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      alpha_zz=alpha_zz-(ak*z*z)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      alpha_xy=alpha_xy-(ak*x*y)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      alpha_xz=alpha_xz-(ak*x*z)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      alpha_yz=alpha_yz-(ak*y*z)
     .                 *((1.0/(dble(omega)-dble(eci(i))))
     .                  -(1.0/(dble(omega)+dble(eci(i)))))
      enddo
      write(*,*)
      write(*,3333) 'Polarizability (',-45.56335/omega,
     .                             ';',45.56335/omega,')'
      write(*,*)
      write(*,1111)'x','y','z'
      write(*,2222)'x',alpha_xx,alpha_xy,alpha_xz
      write(*,2222)'y',alpha_xy,alpha_yy,alpha_yz
      write(*,2222)'z',alpha_xz,alpha_yz,alpha_zz
      write(*,*)
      write(*,111) 'Mean of alpha',(alpha_xx+alpha_yy+alpha_zz)/3.0
      write(*,*)



      enddo
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3333  format(A16,F7.1,A1,F7.1,A1)
      end subroutine pol_sos

      subroutine lresp_ESA(nci,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,Yci,nroot,xmolw,thr)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xmolw

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)
      real*4 ::delta_e(nroot)
      real*8 ::thr,umerk(4,nroot)

      real*8 ::mu_s2s(nroot,3),A1,A2,B1,B2,osc_strength(nroot)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)

      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      call cpu_time(start_time)
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,no,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do ii=1, nroot
      if(ii/=state2opt)write(*,*) 'Transition dipole moment from',
     .                                           state2opt,'to',ii
      !if(ii/=state2opt.and.abs(eci(ii)-eci(state2opt))*27.211>0.1)then
      if(ii==state2opt)write(*,12)'<',ii,'| mu - <0|mu|0> |',state2opt
     .                           ,'>'
      Do ix=1, 3
      A1=0.0
      A2=0.0
      B1=0.0
      B2=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A1,A2)
!$omp do
      Do i=1,counter_A
      A1=A1+dble(Xci(A_list(i,1),ii))*(-mu(A_list(i,2),ix))
     .                    *dble(Xci(A_list(i,3),state2opt))
      A2=A2+dble(Yci(A_list(i,1),state2opt))*(-mu(A_list(i,2),ix))
     .                    *dble(Yci(A_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:B1,B2)
!$omp do
      Do i=1,counter_B
      B1=B1+dble(Xci(B_list(i,1),ii))*(-mu(B_list(i,2),ix))
     .                    *dble(Xci(B_list(i,3),state2opt))
      B2=B2+dble(Yci(B_list(i,1),state2opt))*(-mu(B_list(i,2),ix))
     .                    *dble(Yci(B_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
      mu_s2s(ii,ix)=-(A1+A2-B1-B2)/2.0
      write(*,*) ix, mu_s2s(ii,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(state2opt))*27.211
      !mu_s2s(ii,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV      nm       fL'
      Do ii=1, nroot
      delta_e(ii)=eci(ii)-eci(state2opt)
      osc_strength(ii)=2.0/3.0*dble(delta_e(ii))*
     .(mu_s2s(ii,1)**2.0+mu_s2s(ii,2)**2.0+mu_s2s(ii,3)**2.0)

      write(*,11) state2opt, ii, (delta_e(ii))*27.21139,
     .1.d+7/(delta_e(ii)*2.19474625d+5),
     .osc_strength(ii)
      umerk(1,ii)=osc_strength(ii)
      umerk(2:4,ii)=0.0
      enddo
      call print_tdadat(nroot,xmolw,delta_e,umerk(1:4,:),thr,'s2s.dat')
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f9.1,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine lresp_ESA

      subroutine lresp_ESAbis(nci,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,Yci,nroot,mu_s2s)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)

      real*8 ::mu_s2s(nroot,nroot,3),A1,A2,B1,B2,osc_strength

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)

      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      call cpu_time(start_time)
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,no,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do jj=1,nroot
      Do ii=1, nroot
      !if(ii/=jj)write(*,*) 'Transition dipole moment from',
      !.                                           jj,'to',ii
      !if(ii/=jj.and.abs(eci(ii)-eci(jj))*27.211>0.1)then
      !if(ii==jj)write(*,12)'<',ii,'| mu - <0|mu|0> |',jj
      !.                           ,'>'
      Do ix=1, 3
      A1=0.0
      A2=0.0
      B1=0.0
      B2=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A1,A2)
!$omp do
      Do i=1,counter_A
      A1=A1+dble(Xci(A_list(i,1),ii))*(-mu(A_list(i,2),ix))
     .                    *dble(Xci(A_list(i,3),jj))
      A2=A2+dble(Yci(A_list(i,1),jj))*(-mu(A_list(i,2),ix))
     .                    *dble(Yci(A_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:B1,B2)
!$omp do
      Do i=1,counter_B
      B1=B1+dble(Xci(B_list(i,1),ii))*(-mu(B_list(i,2),ix))
     .                    *dble(Xci(B_list(i,3),jj))
      B2=B2+dble(Yci(B_list(i,1),jj))*(-mu(B_list(i,2),ix))
     .                    *dble(Yci(B_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
      mu_s2s(ii,jj,ix)=-(A1+A2-B1-B2)/2.0
      !write(*,*) ix, mu_s2s(ii,jj,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(jj))*27.211
      !mu_s2s(ii,jj,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV       fL'
      Do ii=1, nroot
      osc_strength=2.0/3.0*(eci(ii)-eci(jj))*
     .(mu_s2s(ii,jj,1)**2.0+mu_s2s(ii,jj,2)**2.0+mu_s2s(ii,jj,3)**2.0)
      write(*,11) jj, ii, (eci(ii)-eci(jj))*27.21139,
     .osc_strength
      enddo
      write(*,*)
      write(*,*)
      enddo
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine lresp_ESAbis


      subroutine hyperpol_sos(nroot,nci,eci,Xci,Yci,xl,yl,zl,moci,
     .maxconf,iconf,no,nv)
      use commonresp
      use omp_lib
      implicit none

      ! This subroutine uses <k|mu-mu_0|n> dipole moments that are unrelaxed,
      !           the result should be different that with response functions

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot,no,nv
      integer ::ix,iy,iz,iroot
      integer ::maxconf,moci
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::omega
      real*4 ::freq(num_freq+1)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)

      real*8 ::beta(3,3,3)
      real*8 ::mu_s(nroot,3)
      real*8 ::mu_s2s(nroot,nroot,3)

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      open(unit=101,file='wavelength',form='formatted',status='old')
      freq(1)=0.0
      num_freq=1
      write(*,*)
      write(*,*) 'Wavelengths (nm)'
      write(*,*) '  infinity'
      Do i=1, num_freq
      read(101,*)freq(i+1)
      write(*,*)freq(i+1)
      freq(i+1)=45.56335/freq(i+1)
      enddo
      close(101)

      Do i=1,nroot
      mu_s(i,:)=0.0
      Do j=1,nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2

      mu_s(i,1)=mu_s(i,1)+mu(ij,1)*(dble(Xci(j,i))+dble(Yci(j,i)))
      mu_s(i,2)=mu_s(i,2)+mu(ij,2)*(dble(Xci(j,i))+dble(Yci(j,i)))
      mu_s(i,3)=mu_s(i,3)+mu(ij,3)*(dble(Xci(j,i))+dble(Yci(j,i)))
      enddo
      mu_s(i,:)=mu_s(i,:)*2.0**(1.0/2.0)
      enddo

      call lresp_ESAbis(nci,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,Yci,nroot,mu_s2s)

      Do iroot=1, num_freq+1
      omega=freq(iroot)
      beta=0.0

      Do ix=1,3
      Do iy=1,3
      Do iz=1,3
      Do ii=1,nroot
      Do jj=1,nroot
      beta(ix,iy,iz)=beta(ix,iy,iz)+
     . mu_s(ii,ix)*mu_s2s(ii,jj,iy)*mu_s(jj,iz)
     . /((dble(omega)*2.0-dble(eci(ii)))
     .  *(dble(omega)-dble(eci(jj))))+
     . mu_s(jj,iz)*mu_s2s(jj,ii,iy)*mu_s(ii,ix)
     . /((dble(omega)*2.0+dble(eci(ii)))
     .  *(dble(omega)+dble(eci(jj))))-
     . mu_s(ii,iy)*mu_s2s(ii,jj,ix)*mu_s(jj,iz)
     . /((dble(omega)+dble(eci(ii)))
     .  *(dble(omega)-dble(eci(jj))))+
     . mu_s(ii,ix)*mu_s2s(ii,jj,iz)*mu_s(jj,iy)
     . /((dble(omega)*2.0-dble(eci(ii)))
     .  *(dble(omega)-dble(eci(jj))))+
     . mu_s(jj,iy)*mu_s2s(jj,ii,iz)*mu_s(ii,ix)
     . /((dble(omega)*2.0+dble(eci(ii)))
     .  *(dble(omega)+dble(eci(jj))))-
     . mu_s(ii,iz)*mu_s2s(ii,jj,ix)*mu_s(jj,iy)
     . /((dble(omega)+dble(eci(ii)))
     .  *(dble(omega)-dble(eci(jj))))
      enddo
      enddo
      enddo
      enddo
      enddo
      write(*,*) '_____________________________________________
     .________'
      write(*,*)
      write(*,5555) 'SHG first hyperpolarizability ('
     .,-45.56335/(freq(iroot)*2.0),';',45.56335/freq(iroot),','
     .,45.56335/freq(iroot),')'
      write(*,*)
      call PrintBeta(beta)
      write(*,*) '_____________________________________________
     .________'
      enddo
5555  format(A31,F7.1,A1,F7.1,A1,F7.1,A1)
      end subroutine hyperpol_sos


      subroutine tpa_sos(nroot,nci,eci,Xci,Yci,xl,yl,zl,moci,
     .maxconf,iconf,no,nv)
      use commonresp
      use omp_lib
      implicit none

      ! This subroutine uses <k|mu-mu_0|n> dipole moments that are unrelaxed,
      !           the result should be different that with response functions

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot,no,nv
      integer ::ix,iy,iz,iroot
      integer ::maxconf,moci
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::omega
      real*4 ::freq(num_freq+1)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)

      real*8 ::sigma(3,3)
      real*8 ::mu_s(nroot,3)
      real*8 ::mu_s2s(nroot,nroot,3)

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      Do i=1,nroot
      mu_s(i,:)=0.0
      Do j=1,nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2

      mu_s(i,1)=mu_s(i,1)+mu(ij,1)*(dble(Xci(j,i))+dble(Yci(j,i)))
      mu_s(i,2)=mu_s(i,2)+mu(ij,2)*(dble(Xci(j,i))+dble(Yci(j,i)))
      mu_s(i,3)=mu_s(i,3)+mu(ij,3)*(dble(Xci(j,i))+dble(Yci(j,i)))
      enddo
      mu_s(i,:)=mu_s(i,:)*2.0**(1.0/2.0)
      enddo

      call lresp_ESAbis(nci,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,Yci,nroot,mu_s2s)

      Do iroot=1, nroot
      omega=eci(iroot)/2.0
      sigma=0.0

      Do ix=1,3
      Do iy=1,3
      Do ii=1,nroot
      sigma(ix,iy)=sigma(ix,iy)-
     . (mu_s(ii,ix)*mu_s2s(ii,iroot,iy)/(dble(eci(ii))-dble(omega))+
     . mu_s(ii,iy)*mu_s2s(ii,iroot,ix)/(dble(eci(ii))-dble(omega)))
      enddo
      enddo
      enddo
      write(*,*)
      write(*,3333)'Delta (',eci(iroot)*27.21139,')'
      write(*,*)
      write(*,1111)'x','y','z'
      write(*,2222)'x',sigma(1,1),sigma(1,2),sigma(1,3)
      write(*,2222)'y',sigma(2,1),sigma(2,2),sigma(2,3)
      write(*,2222)'z',sigma(3,1),sigma(3,2),sigma(3,3)
      write(*,*)
      enddo
3333  format(A16,F7.3,A1)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
      end subroutine tpa_sos

      subroutine lresp_ESA_tda(nci,iconf,maxconf,xl,yl,zl,moci,
     .no,nv,eci,Xci,nroot,xmolw,thr)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xmolw

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*8 ::mu(moci*(moci+1)/2,3)
      real*4 ::Xci(nci,nroot),eci(nci)
      real*4 ::delta_e(nroot)
      real*8 ::thr,umerk(4,nroot)

      real*8 ::mu_s2s(nroot,3),A1,B1,osc_strength(nroot)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mu=0.0
      mu(:,1)=xl(:)
      mu(:,2)=yl(:)
      mu(:,3)=zl(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTDA
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)

      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      call cpu_time(start_time)
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,no
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,no,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nv
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+no)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,no,nv,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do ii=1, nroot
      if(ii/=state2opt)write(*,*) 'Transition dipole moment from',
     .                                           state2opt,'to',ii
      !if(ii/=state2opt.and.abs(eci(ii)-eci(state2opt))*27.211>0.1)then
      if(ii==state2opt)write(*,12)'<',ii,'| mu - <0|mu|0> |',state2opt
     .                           ,'>'
      Do ix=1, 3
      A1=0.0
      B1=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A1)
!$omp do
      Do i=1,counter_A
      A1=A1+dble(Xci(A_list(i,1),ii))*(-mu(A_list(i,2),ix))
     .                    *dble(Xci(A_list(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:B1)
!$omp do
      Do i=1,counter_B
      B1=B1+dble(Xci(B_list(i,1),ii))*(-mu(B_list(i,2),ix))
     .                    *dble(Xci(B_list(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel
      mu_s2s(ii,ix)=-(A1-B1)/2.0
      write(*,*) ix, mu_s2s(ii,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(state2opt))*27.211
      !mu_s2s(ii,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV      nm       fL'
      Do ii=1, nroot
      delta_e(ii)=eci(ii)-eci(state2opt)
      osc_strength(ii)=2.0/3.0*dble(delta_e(ii))*
     .(mu_s2s(ii,1)**2.0+mu_s2s(ii,2)**2.0+mu_s2s(ii,3)**2.0)

      write(*,11) state2opt, ii, (delta_e(ii))*27.21139,
     .1.d+7/(delta_e(ii)*2.19474625d+5),
     .osc_strength(ii)
      umerk(1,ii)=osc_strength(ii)
      umerk(2:4,ii)=0.0
      enddo
      call print_tdadat(nroot,xmolw,delta_e,umerk(1:4,:),thr,'s2s.dat')
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTDA
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f9.1,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine lresp_ESA_tda


      subroutine ulresp_ESA(nexa,nexb,nci,iconfa,iconfb,maxconfa,
     .                maxconfb,xla,yla,zla,xlb,ylb,zlb,mocia,mocib,
     .              noa,nob,nva,nvb,eci,Xci,Yci,nroot,xmolw,thr,ak)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::nexa,nexb,maxconfa,maxconfb,mocia,mocib
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::noa,nva,nob,nvb
      integer ::iconfa(maxconfa,2),iconfb(maxconfb,2)
      integer, allocatable :: A_lista(:,:)
      integer ::counter_Aa
      integer, allocatable :: B_lista(:,:)
      integer ::counter_Ba
      integer, allocatable :: A_listb(:,:)
      integer ::counter_Ab
      integer, allocatable :: B_listb(:,:)
      integer ::counter_Bb
      integer*8 ::lin8

      real*8 ::xmolw,ak

      real*8 ::xla(mocia*(mocia+1)/2)
      real*8 ::yla(mocia*(mocia+1)/2)
      real*8 ::zla(mocia*(mocia+1)/2)

      real*8 ::xlb(mocib*(mocib+1)/2)
      real*8 ::ylb(mocib*(mocib+1)/2)
      real*8 ::zlb(mocib*(mocib+1)/2)

      real*8 ::mua(mocia*(mocia+1)/2,3)
      real*8 ::mub(mocib*(mocib+1)/2,3)
      real*4 ::Xci(nci,nroot), Yci(nci,nroot),eci(nci)
      real*4 ::delta_e(nroot)
      real*8 ::thr,umerk(4,nroot)

      real*8 ::mu_s2s(nroot,3),Aa1,Aa2,Ba1,Ba2,osc_strength(nroot)
      real*8 ::Ab1,Ab2,Bb1,Bb2
      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mua=0.0
      mub=0.0
      mua(:,1)=xla(:)
      mua(:,2)=yla(:)
      mua(:,3)=zla(:)
      mub(:,1)=xlb(:)
      mub(:,2)=ylb(:)
      mub(:,3)=zlb(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)



      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !

            ! alpha block


      call cpu_time(start_time)
      counter_Aa=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Aa)
!$omp do
      Do i=1,nexa
      Do j=1,noa
        Do kk=1,nexa
        if(iconfa(kk,1)==j .and. iconfa(kk,2)==iconfa(i,2))then
            counter_Aa=counter_Aa+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_lista(1:counter_Aa,1:3))
      A_lista=-9999
      call List_A(maxconfa,noa,nexa,iconfa,A_lista,counter_Aa)


      counter_Ba=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Ba)
!$omp do
      Do i=1,nexa
      Do j=1,nva
        Do kk=1,nexa
        if(iconfa(kk,1)==iconfa(i,1) .and. iconfa(kk,2)==j+noa)then
            counter_Ba=counter_Ba+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_lista(1:counter_Ba,1:3))
      B_lista=-9999
      call List_B(maxconfa,noa,nva,nexa,iconfa,B_lista,counter_Ba)

            ! beta block


      call cpu_time(start_time)
      counter_Ab=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Ab)
!$omp do
      Do i=1,nexb
      Do j=1,nob
        Do kk=1,nexb
        if(iconfb(kk,1)==j .and. iconfb(kk,2)==iconfb(i,2))then
            counter_Ab=counter_Ab+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_listb(1:counter_Ab,1:3))
      A_listb=-9999
      call List_A(maxconfb,nob,nexb,iconfb,A_listb,counter_Ab)


      counter_Bb=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Bb)
!$omp do
      Do i=1,nexb
      Do j=1,nvb
        Do kk=1,nexb
        if(iconfb(kk,1)==iconfb(i,1) .and. iconfb(kk,2)==j+nob)then
            counter_Bb=counter_Bb+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_listb(1:counter_Bb,1:3))
      B_listb=-9999
      call List_B(maxconfb,nob,nvb,nexb,iconfb,B_listb,counter_Bb)


      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do ii=1, nroot
      if(ii/=state2opt)write(*,*) 'Transition dipole moment from',
     .                                           state2opt,'to',ii
      !if(ii/=state2opt.and.abs(eci(ii)-eci(state2opt))*27.211>0.1)then
      if(ii==state2opt)write(*,12)'<',ii,'| mu - <0|mu|0> |',state2opt
     .                           ,'>'
      Do ix=1, 3
      Aa1=0.0
      Aa2=0.0
      Ba1=0.0
      Ba2=0.0
      Ab1=0.0
      Ab2=0.0
      Bb1=0.0
      Bb2=0.0

      ! alpha


!$omp parallel private(i)
!$omp&         reduction(+:Aa1,Aa2)
!$omp do
      Do i=1,counter_Aa
      Aa1=Aa1+dble(Xci(A_lista(i,1),ii))*(-mua(A_lista(i,2),ix))
     .                    *dble(Xci(A_lista(i,3),state2opt))
      Aa2=Aa2
     . +dble(Yci(A_lista(i,1),state2opt))*(-mua(A_lista(i,2),ix))
     .                    *dble(Yci(A_lista(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:Ba1,Ba2)
!$omp do
      Do i=1,counter_Ba
      Ba1=Ba1+dble(Xci(B_lista(i,1),ii))*(-mua(B_lista(i,2),ix))
     .                    *dble(Xci(B_lista(i,3),state2opt))
      Ba2=Ba2
     . +dble(Yci(B_lista(i,1),state2opt))*(-mua(B_lista(i,2),ix))
     .                    *dble(Yci(B_lista(i,3),ii))
      enddo
!$omp end do
!$omp end parallel


      ! beta


!$omp parallel private(i)
!$omp&         reduction(+:Ab1,Ab2)
!$omp do
      Do i=1,counter_Ab
      Ab1=Ab1+dble(Xci(A_listb(i,1)+nexa,ii))*(-mub(A_listb(i,2),ix))
     .                    *dble(Xci(A_listb(i,3)+nexa,state2opt))
      Ab2=Ab2
     . +dble(Yci(A_listb(i,1)+nexa,state2opt))*(-mub(A_listb(i,2),ix))
     .                    *dble(Yci(A_listb(i,3)+nexa,ii))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:Bb1,Bb2)
!$omp do
      Do i=1,counter_Bb
      Bb1=Bb1+dble(Xci(B_listb(i,1)+nexa,ii))*(-mub(B_listb(i,2),ix))
     .                    *dble(Xci(B_listb(i,3)+nexa,state2opt))
      Bb2=Bb2
     . +dble(Yci(B_listb(i,1)+nexa,state2opt))*(-mub(B_listb(i,2),ix))
     .                    *dble(Yci(B_listb(i,3)+nexa,ii))
      enddo
!$omp end do
!$omp end parallel

      ! divided by two because the beta formula is divided by two in the unrestricted case


      mu_s2s(ii,ix)=-(Aa1+Aa2-Ba1-Ba2+Ab1+Ab2-Bb1-Bb2)/2.0
      write(*,*) ix, mu_s2s(ii,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(state2opt))*27.211
      !mu_s2s(ii,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV      nm       fL'
      Do ii=1, nroot
      delta_e(ii)=eci(ii)-eci(state2opt)
      osc_strength(ii)=2.0/3.0*dble(delta_e(ii))*
     .(mu_s2s(ii,1)**2.0+mu_s2s(ii,2)**2.0+mu_s2s(ii,3)**2.0)

      write(*,11) state2opt, ii, (delta_e(ii))*27.21139,
     .1.d+7/(delta_e(ii)*2.19474625d+5),
     .osc_strength(ii)
      umerk(1,ii)=osc_strength(ii)
      umerk(2:4,ii)=0.0
      enddo
      call print_tdadat(nroot,xmolw,delta_e,umerk(1:4,:),thr,'s2s.dat')
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f9.1,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine ulresp_ESA

      subroutine ulresp_ESA_tda(nexa,nexb,nci,iconfa,iconfb,
     .      maxconfa,maxconfb,xla,yla,zla,xlb,ylb,zlb,mocia,mocib,
     .              noa,nob,nva,nvb,eci,Xci,nroot,xmolw,thr,ak)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::nexa,nexb,maxconfa,maxconfb,mocia,mocib
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::noa,nva,nob,nvb
      integer ::iconfa(maxconfa,2),iconfb(maxconfb,2)
      integer, allocatable :: A_lista(:,:)
      integer ::counter_Aa
      integer, allocatable :: B_lista(:,:)
      integer ::counter_Ba
      integer, allocatable :: A_listb(:,:)
      integer ::counter_Ab
      integer, allocatable :: B_listb(:,:)
      integer ::counter_Bb
      integer*8 ::lin8

      real*8 ::xmolw,ak

      real*8 ::xla(mocia*(mocia+1)/2)
      real*8 ::yla(mocia*(mocia+1)/2)
      real*8 ::zla(mocia*(mocia+1)/2)

      real*8 ::xlb(mocib*(mocib+1)/2)
      real*8 ::ylb(mocib*(mocib+1)/2)
      real*8 ::zlb(mocib*(mocib+1)/2)

      real*8 ::mua(mocia*(mocia+1)/2,3)
      real*8 ::mub(mocib*(mocib+1)/2,3)
      real*4 ::Xci(nci,nroot),eci(nci)
      real*4 ::delta_e(nroot)
      real*8 ::thr,umerk(4,nroot)

      real*8 ::mu_s2s(nroot,3),Aa1,Ba1,osc_strength(nroot)
      real*8 ::Ab1,Bb1
      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      mua=0.0
      mub=0.0
      mua(:,1)=xla(:)
      mua(:,2)=yla(:)
      mua(:,3)=zla(:)
      mub(:,1)=xlb(:)
      mub(:,2)=ylb(:)
      mub(:,3)=zlb(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response sTDA
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)



      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !

            ! alpha block


      call cpu_time(start_time)
      counter_Aa=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Aa)
!$omp do
      Do i=1,nexa
      Do j=1,noa
        Do kk=1,nexa
        if(iconfa(kk,1)==j .and. iconfa(kk,2)==iconfa(i,2))then
            counter_Aa=counter_Aa+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_lista(1:counter_Aa,1:3))
      A_lista=-9999
      call List_A(maxconfa,noa,nexa,iconfa,A_lista,counter_Aa)


      counter_Ba=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Ba)
!$omp do
      Do i=1,nexa
      Do j=1,nva
        Do kk=1,nexa
        if(iconfa(kk,1)==iconfa(i,1) .and. iconfa(kk,2)==j+noa)then
            counter_Ba=counter_Ba+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_lista(1:counter_Ba,1:3))
      B_lista=-9999
      call List_B(maxconfa,noa,nva,nexa,iconfa,B_lista,counter_Ba)

            ! beta block


      call cpu_time(start_time)
      counter_Ab=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Ab)
!$omp do
      Do i=1,nexb
      Do j=1,nob
        Do kk=1,nexb
        if(iconfb(kk,1)==j .and. iconfb(kk,2)==iconfb(i,2))then
            counter_Ab=counter_Ab+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_listb(1:counter_Ab,1:3))
      A_listb=-9999
      call List_A(maxconfb,nob,nexb,iconfb,A_listb,counter_Ab)


      counter_Bb=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_Bb)
!$omp do
      Do i=1,nexb
      Do j=1,nvb
        Do kk=1,nexb
        if(iconfb(kk,1)==iconfb(i,1) .and. iconfb(kk,2)==j+nob)then
            counter_Bb=counter_Bb+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_listb(1:counter_Bb,1:3))
      B_listb=-9999
      call List_B(maxconfb,nob,nvb,nexb,iconfb,B_listb,counter_Bb)


      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do ii=1, nroot
      if(ii/=state2opt)write(*,*) 'Transition dipole moment from',
     .                                           state2opt,'to',ii
      !if(ii/=state2opt.and.abs(eci(ii)-eci(state2opt))*27.211>0.1)then
      if(ii==state2opt)write(*,12)'<',ii,'| mu - <0|mu|0> |',state2opt
     .                           ,'>'
      Do ix=1, 3
      Aa1=0.0
      Ba1=0.0
      Ab1=0.0
      Bb1=0.0

      ! alpha


!$omp parallel private(i)
!$omp&         reduction(+:Aa1)
!$omp do
      Do i=1,counter_Aa
      Aa1=Aa1+dble(Xci(A_lista(i,1),ii))*(-mua(A_lista(i,2),ix))
     .                    *dble(Xci(A_lista(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:Ba1)
!$omp do
      Do i=1,counter_Ba
      Ba1=Ba1+dble(Xci(B_lista(i,1),ii))*(-mua(B_lista(i,2),ix))
     .                    *dble(Xci(B_lista(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel


      ! beta


!$omp parallel private(i)
!$omp&         reduction(+:Ab1)
!$omp do
      Do i=1,counter_Ab
      Ab1=Ab1+dble(Xci(A_listb(i,1)+nexa,ii))*(-mub(A_listb(i,2),ix))
     .                    *dble(Xci(A_listb(i,3)+nexa,state2opt))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:Bb1)
!$omp do
      Do i=1,counter_Bb
      Bb1=Bb1+dble(Xci(B_listb(i,1)+nexa,ii))*(-mub(B_listb(i,2),ix))
     .                    *dble(Xci(B_listb(i,3)+nexa,state2opt))
      enddo
!$omp end do
!$omp end parallel

      ! divided by two because the beta formula is divided by two in the unrestricted case

      mu_s2s(ii,ix)=-(Aa1-Ba1+Ab1-Bb1)/2.0
      write(*,*) ix, mu_s2s(ii,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(state2opt))*27.211
      !mu_s2s(ii,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV      nm       fL'
      Do ii=1, nroot
      delta_e(ii)=eci(ii)-eci(state2opt)
      osc_strength(ii)=2.0/3.0*dble(delta_e(ii))*
     .(mu_s2s(ii,1)**2.0+mu_s2s(ii,2)**2.0+mu_s2s(ii,3)**2.0)

      write(*,11) state2opt, ii, (delta_e(ii))*27.21139,
     .1.d+7/(delta_e(ii)*2.19474625d+5),
     .osc_strength(ii)
      umerk(1,ii)=osc_strength(ii)
      umerk(2:4,ii)=0.0
      enddo
      call print_tdadat(nroot,xmolw,delta_e,umerk(1:4,:),thr,'s2s.dat')
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response sTDA
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f9.1,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine ulresp_ESA_tda

      subroutine sf_lresp_ESA(nci,iconf,maxconf,xla,yla,zla,mocia,
     .                                          xlb,ylb,zlb,mocib,
     .                    noa,nva,nob,nvb,eci,Xci,nroot,xmolw,thr)
      use commonresp
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot
      integer ::maxconf,mocia,noa,nva,mocib,nob,nvb
      integer ::iconf(maxconf,2)
      integer, allocatable :: A_list(:,:)
      integer ::counter_A
      integer, allocatable :: B_list(:,:)
      integer ::counter_B
      integer*8 ::lin8

      real*8 ::xmolw

      real*8 ::xla(mocia*(mocia+1)/2)
      real*8 ::yla(mocia*(mocia+1)/2)
      real*8 ::zla(mocia*(mocia+1)/2)

      real*8 ::xlb(mocib*(mocib+1)/2)
      real*8 ::ylb(mocib*(mocib+1)/2)
      real*8 ::zlb(mocib*(mocib+1)/2)

      real*8 ::mua(mocia*(mocia+1)/2,3)
      real*8 ::mub(mocib*(mocib+1)/2,3)
      real*4 ::Xci(nci,nroot),eci(nci)
      real*4 ::delta_e(nroot)
      real*8 ::thr,umerk(4,nroot)

      real*8 ::mu_s2s(nroot,3),A1,B1,osc_strength(nroot)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      !first spinflip state
      state2opt=1

      mua=0.0
      mua(:,1)=xla(:)
      mua(:,2)=yla(:)
      mua(:,3)=zla(:)
      mub=0.0
      mub(:,1)=xlb(:)
      mub(:,2)=ylb(:)
      mub(:,3)=zlb(:)

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in nonlinear response SF-sTDDFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      write(*,*)'Useful if the first SF state is the ground state'
      write(*,*)
      !
      ! Genarating a list of indexes used in A and B formula to save a great bunch of time
      !
      call cpu_time(start_time)
      counter_A=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_A)
!$omp do
      Do i=1,nci
      Do j=1,noa
        Do kk=1,nci
        if(iconf(kk,1)==j .and. iconf(kk,2)==iconf(i,2))then
            counter_A=counter_A+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(A_list(1:counter_A,1:3))
      A_list=-9999
      call List_A(maxconf,noa,nci,iconf,A_list,counter_A)
      !Do i=1,counter_A
      !write(*,*)i,A_list(i,1:3)
      !enddo

      counter_B=0
!$omp parallel private(j,i,kk)
!$omp&         reduction(+:counter_B)
!$omp do
      Do i=1,nci
      Do j=1,nvb
        Do kk=1,nci
        if(iconf(kk,1)==iconf(i,1) .and. iconf(kk,2)==j+nob)then
            counter_B=counter_B+1
        endif
        enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(B_list(1:counter_B,1:3))
      B_list=-9999
      call List_B(maxconf,nob,nvb,nci,iconf,B_list,counter_B)
      !Do i=1,counter_B
      !write(*,*)i,B_list(i,1:3)
      !enddo
      call cpu_time(end_time)
      print '("A & B indexes list Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*)
      mu_s2s=0.0
      Do ii=1, nroot
!       if(ii/=state2opt)write(*,*) 'Transition dipole moment from',
!      .                                           state2opt,'to',ii
!       !if(ii/=state2opt.and.abs(eci(ii)-eci(state2opt))*27.211>0.1)then
!       if(ii==state2opt)write(*,12)'<',ii,'| mu - <0|mu|0> |',state2opt
!      .                           ,'>'
      Do ix=1, 3
      A1=0.0
      B1=0.0
!$omp parallel private(i)
!$omp&         reduction(+:A1)
!$omp do
      Do i=1,counter_A
      A1=A1+dble(Xci(A_list(i,1),ii))*(-mua(A_list(i,2),ix))
     .                    *dble(Xci(A_list(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i)
!$omp&         reduction(+:B1)
!$omp do
      Do i=1,counter_B
      B1=B1+dble(Xci(B_list(i,1),ii))*(-mub(B_list(i,2),ix))
     .                    *dble(Xci(B_list(i,3),state2opt))
      enddo
!$omp end do
!$omp end parallel
      mu_s2s(ii,ix)=-(A1-B1)*sqrt(2.0)/2.0 ! Unrestricted formula so 1/2, with alpha -> beta, missing beta -> alpha so sqrt(2)
!       write(*,*) ix, mu_s2s(ii,ix)
      enddo
      !else
      !write(*,*) 'transition rejected'
      !write(*,*) 'delta_E', abs(eci(ii)-eci(state2opt))*27.211
      !mu_s2s(ii,ix)=0.0
      !endif
      enddo
      write(*,*)
      write(*,*) 'State to state      eV      nm       fL'
      Do ii=1, nroot
      delta_e(ii)=eci(ii)-eci(state2opt)
      osc_strength(ii)=2.0/3.0*dble(delta_e(ii))*
     .(mu_s2s(ii,1)**2.0+mu_s2s(ii,2)**2.0+mu_s2s(ii,3)**2.0)

      write(*,11) state2opt, ii, (delta_e(ii))*27.21139,
     .1.d+7/(delta_e(ii)*2.19474625d+5),
     .osc_strength(ii)
      umerk(1,ii)=osc_strength(ii)
      umerk(2:4,ii)=0.0
      enddo
      call print_tdadat(nroot,xmolw,delta_e,umerk(1:4,:),thr,'s2s.dat')
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               end of     nonlinear response SF-sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
 11   format(i6,i9,f9.3,f9.1,f11.4)
 12   format(a1,i3,a17,i3,a1)

      end subroutine sf_lresp_ESA

      SUBROUTINE optrot(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,no,nv,
     .                  xm,ym,zm,xmass)
      use commonresp
      use commonlogicals
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xl(moci*(moci+1)/2)
      real*8 ::yl(moci*(moci+1)/2)
      real*8 ::zl(moci*(moci+1)/2)

      real*4 ::mu_x(nci)
      real*4 ::mu_y(nci)
      real*4 ::mu_z(nci)

      real*8 ::xm(moci*(moci+1)/2)
      real*8 ::ym(moci*(moci+1)/2)
      real*8 ::zm(moci*(moci+1)/2)
      real*8 ::xmass,refindex,vorfaktor,reffaktor,refval
      logical::da

      real*4 ::apb(nci*(nci+1)/2)
      real*4 ::amb(nci*(nci+1)/2)
      real*4, allocatable ::inv_amb(:)
      real*4, allocatable ::inv_resp(:)
      real*4, allocatable ::XpY(:,:),XmY(:,:)
      real*4 ::omega
      real*4, allocatable ::freq(:)
      real*4 ::alpha_xx,alpha_xy,alpha_xz
      real*4 ::alpha_yy,alpha_yz
      real*4 ::alpha_zz
      character*1 ::uplo
      integer ::info
      integer, allocatable ::ipiv(:)
      real*4, allocatable ::work (:)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      character(len=14):: dummy

      logical :: file_exists

      integer :: nlines

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in linear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='

      INQUIRE(FILE="wavelength", EXIST=file_exists)
      if(file_exists==.false.)then
      num_freq=1
      allocate(freq(1))
      freq(1)=45.56335/589.30
      else
      nlines = 0
      open(unit=101,file='wavelength',form='formatted')
      DO
      READ (101,*, END=10)
      nlines = nlines + 1
      END DO
   10 rewind (101)
      num_freq=nlines
      allocate(freq(num_freq))
      Do i=1, num_freq
      read(101,*)freq(i)
      write(*,*)freq(i)
      freq(i)=45.56335/freq(i)
      enddo
      close(101)
      endif

      write(*,*) 'Optical rotation alpha[grad*cm^3*g^-1*dm^-1]'
      write(*,*) 'including Lorentz factor for common solvent (n=1.4)'
      write(*,*) 'lambda [eV]       alpha       alpha(n=1.0) Phi
     .         Phi(n=1.0)'
      ! shift


c refractive index of solvent
      refindex=1.4d0
      reffaktor=(refindex**2+2.0d0)/3.0d0

      allocate(inv_amb(nci*(nci+1)/2))
      inv_amb=amb
      uplo='U'
      allocate(ipiv(1:nci),work(1:nci))
      call ssptrf(uplo,nci,inv_amb,ipiv,info)
      call ssptri(uplo,nci,inv_amb,ipiv,work,info)
      deallocate(ipiv,work)

      allocate( XpY(nci,3))
      allocate( XmY(nci,3))
      ! the dipole moment matrix mu_ai
      mu_x=0.0
      mu_y=0.0
      mu_z=0.0
!$omp parallel private(j,io,iv,idum1,idum2,ij)
!$omp do
              Do j=1, nci
                  io=iconf(j,1)
                  iv=iconf(j,2)
                  idum1=max(io,iv)
                  idum2=min(io,iv)
                  ij=idum2+idum1*(idum1-1)/2
            mu_x(j)=-2.0*xl(ij)
            mu_y(j)=-2.0*yl(ij)
            mu_z(j)=-2.0*zl(ij)
              enddo
!$omp end do
!$omp end parallel

      Do ii=1, num_freq
      omega=freq(ii)
      XpY(:,:)=0.0
      call cpu_time(start_time)
      allocate(inv_resp(nci*(nci+1)/2))
      inv_resp=apb-omega**2.0*inv_amb

      uplo='U'
      XpY(:,1)=mu_x(:)
      XpY(:,2)=mu_y(:)
      XpY(:,3)=mu_z(:)
      allocate(ipiv(1:nci))
      call ssptrf(uplo,nci,inv_resp,ipiv,info)
      call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY,nci,info)
      deallocate(ipiv)

      !(X-Y)=omega*(A-B)^-1 (X+Y)
      XmY=0.0
!$omp parallel private(i,j,ij,ix)
!$omp&                 reduction (+:XmY)
!$omp do
      Do i=1,nci
        Do j=1,nci
        ij=lin8(i,j)
        Do ix=1,3
        XmY(i,ix)=XmY(i,ix)+
     .     inv_amb(ij)*XpY(j,ix)!*dble(omega) !because beta=-Tr(G)/omega/3
        enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
      if(nto)then
      write(dummy,'(a,i0)')'1-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,1)*omega
      close(14)
      write(dummy,'(a,i0)')'2-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,2)*omega
      close(14)
      write(dummy,'(a,i0)')'3-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,3)*omega
      close(14)
      endif
      alpha_xx=0.0
      alpha_yy=0.0
      alpha_zz=0.0
!$omp parallel private(i,io,iv,idum1,idum2,ij)
!$omp&                 reduction(+:alpha_xx,alpha_yy,alpha_zz)
!$omp do
      Do j=1, nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      alpha_xx=alpha_xx+(xm(ij)*XmY(j,1))*2.0 ! minus included
      alpha_yy=alpha_yy+(ym(ij)*XmY(j,2))*2.0
      alpha_zz=alpha_zz+(zm(ij)*XmY(j,3))*2.0
      enddo
!$omp end do
!$omp end parallel
      alpha_xx=alpha_xx/3.0
      alpha_yy=alpha_yy/3.0
      alpha_zz=alpha_zz/3.0
      !G_xx= omega*-alpha_xx
      !beta=-Tr(G)/omega/3=Tr(alpha)/3

      vorfaktor=(38652./xmass)*(589.3/(45.56335/omega))**2
      if(ii==8)then
      refval=0
      inquire(file='.ref',exist=da)
      if(da)then
      open(unit=33,file='.ref')
      read(33,*)refval
      close(33)
      endif
      write(*,142)45.56335/omega,omega*27.21139,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor*xmass/100.0,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor*xmass/100.0,
     .refval,' Na D-line '
      else
      write(*,143)45.56335/omega,omega*27.21139,
     .237.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor,
     .237.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor*xmass/100.0,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor*xmass/100.0
      endif

      deallocate(inv_resp)
      enddo
      deallocate(XpY,XmY)
      deallocate(inv_amb)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3333  format(A16,F7.1,A1,F7.1,A1)
 142  format(f6.1,f6.2,5f12.2,a11)
 143  format(f6.1,f6.2,4f12.2)
      end subroutine optrot

      SUBROUTINE optrot_velo(nci,apb,amb,iconf,maxconf,xv,yv,zv,moci,
     .                  no,nv,
     .                  xm,ym,zm,xmass)
      use commonresp
      use commonlogicals
      use omp_lib
      IMPLICIT NONE

      integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
      integer ::io1,io2,iv1,iv2,iwrk,jwrk
      integer ::maxconf,moci,no,nv
      integer ::iconf(maxconf,2)
      integer*8 ::lin8

      real*8 ::xv(moci*(moci+1)/2)
      real*8 ::yv(moci*(moci+1)/2)
      real*8 ::zv(moci*(moci+1)/2)

      real*8 ::xvelo(moci*(moci+1)/2)
      real*8 ::yvelo(moci*(moci+1)/2)
      real*8 ::zvelo(moci*(moci+1)/2)

      real*4 ::mu_x(nci)
      real*4 ::mu_y(nci)
      real*4 ::mu_z(nci)

      real*8 ::xm(moci*(moci+1)/2)
      real*8 ::ym(moci*(moci+1)/2)
      real*8 ::zm(moci*(moci+1)/2)
      real*8 ::xmass,refindex,vorfaktor,reffaktor,refval
      logical::da

      real*4 ::apb(nci*(nci+1)/2)
      real*4 ::amb(nci*(nci+1)/2)
      real*4, allocatable ::inv_amb(:)
      real*4, allocatable ::inv_resp(:)
      real*4, allocatable ::XpY(:,:),XmY(:,:)
      real*4 ::omega
      real*4, allocatable ::freq(:)
      real*4 ::alpha_xx,alpha_xy,alpha_xz
      real*4 ::alpha_yy,alpha_yz
      real*4 ::alpha_zz
      character*1 ::uplo
      integer ::info
      integer, allocatable ::ipiv(:)
      real*4, allocatable ::work (:)

      integer ::ix,iy,iz

      real*4 ::start_time,end_time

      character(len=14):: dummy

      logical :: file_exists

      integer :: nlines

      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in linear response sTD-DFT
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)'VELOCITY REPRESENTATION (ORIGIN INDEPENDENT)'
      INQUIRE(FILE="wavelength", EXIST=file_exists)
      if(file_exists==.false.)then
      num_freq=1
      allocate(freq(1))
      freq(1)=45.56335/589.30
      else
      nlines = 0
      open(unit=101,file='wavelength',form='formatted')
      DO
      READ (101,*, END=10)
      nlines = nlines + 1
      END DO
   10 rewind (101)
      num_freq=nlines
      allocate(freq(num_freq))
      Do i=1, num_freq
      read(101,*)freq(i)
      write(*,*)freq(i)
      freq(i)=45.56335/freq(i)
      enddo
      close(101)
      endif
      write(*,*) 'Optical rotation alpha[grad*cm^3*g^-1*dm^-1]'
      write(*,*) 'including Lorentz factor for common solvent (n=1.4)'
      write(*,*) 'lambda [eV]       alpha       alpha(n=1.0) Phi
     .         Phi(n=1.0)'

c refractive index of solvent
      refindex=1.4d0
      reffaktor=(refindex**2+2.0d0)/3.0d0

      allocate(inv_amb(nci*(nci+1)/2))
      inv_amb=amb
      uplo='U'
      allocate(ipiv(1:nci),work(1:nci))
      call ssptrf(uplo,nci,inv_amb,ipiv,info)
      call ssptri(uplo,nci,inv_amb,ipiv,work,info)
      deallocate(ipiv,work)

      ! mu_v= -(A-B)**(-1)*nabla

      xvelo=0.0
      yvelo=0.0
      zvelo=0.0
!$omp parallel private(i,j,jk,io,iv,idum1,idum2,ij,ab)
!$omp&                 reduction (+:xvelo,yvelo,zvelo)
!$omp do
      Do i=1, nci
            io=iconf(i,1)
            iv=iconf(i,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ab=idum2+idum1*(idum1-1)/2
        Do j=1, nci
            jk=lin8(i,j)
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      xvelo(ab)=xvelo(ab)+xv(ij)*dble(inv_amb(jk))
      yvelo(ab)=yvelo(ab)+yv(ij)*dble(inv_amb(jk))
      zvelo(ab)=zvelo(ab)+zv(ij)*dble(inv_amb(jk))
        enddo
      enddo
!$omp end do
!$omp end parallel

! the dipole moment matrix mu_ai
      mu_x=0.0
      mu_y=0.0
      mu_z=0.0
!$omp parallel private(j,io,iv,idum1,idum2,ij)
!$omp do
        Do j=1, nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      mu_x(j)=-2.0*xvelo(ij)
      mu_y(j)=-2.0*yvelo(ij)
      mu_z(j)=-2.0*zvelo(ij)
        enddo
!$omp end do
!$omp end parallel

      allocate( XpY(nci,3))
      allocate( XmY(nci,3))

      Do ii=1, num_freq
      omega=freq(ii)
      XpY(:,:)=0.0
      call cpu_time(start_time)
      allocate(inv_resp(nci*(nci+1)/2))
      inv_resp=apb-omega**2.0*inv_amb


      uplo='U'
      XpY(:,1)=mu_x(:)
      XpY(:,2)=mu_y(:)
      XpY(:,3)=mu_z(:)
      allocate(ipiv(1:nci))
      call ssptrf(uplo,nci,inv_resp,ipiv,info)
      call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY,nci,info)
      deallocate(ipiv)

      !(X-Y)=omega*(A-B)^-1 (X+Y)
      XmY=0.0
!$omp parallel private(i,j,ij,ix)
!$omp&                 reduction (+:XmY)
!$omp do
      Do i=1,nci
        Do j=1,nci
        ij=lin8(i,j)
        Do ix=1,3
        XmY(i,ix)=XmY(i,ix)+
     .     inv_amb(ij)*XpY(j,ix)!*dble(omega) because beta=-Tr(G)/omega/3
        enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
      if(nto)then
      write(dummy,'(a,i0)')'1-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,1)*omega
      close(14)
      write(dummy,'(a,i0)')'2-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,2)*omega
      close(14)
      write(dummy,'(a,i0)')'3-',ii
      open(unit=14,file=dummy)
      write(14,*)XmY(:,3)*omega
      close(14)
      endif
      alpha_xx=0.0
      alpha_yy=0.0
      alpha_zz=0.0
!$omp parallel private(i,io,iv,idum1,idum2,ij)
!$omp&                 reduction(+:alpha_xx,alpha_yy,alpha_zz)
!$omp do
      Do j=1, nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
      alpha_xx=alpha_xx+(xm(ij)*XmY(j,1))*2.0 ! minus included
      alpha_yy=alpha_yy+(ym(ij)*XmY(j,2))*2.0
      alpha_zz=alpha_zz+(zm(ij)*XmY(j,3))*2.0
      enddo
!$omp end do
!$omp end parallel
      alpha_xx=alpha_xx/3.0
      alpha_yy=alpha_yy/3.0
      alpha_zz=alpha_zz/3.0
      !G_xx= omega*-alpha_xx
      !beta=-Tr(G)/omega/3=Tr(alpha)/3

      vorfaktor=(38652./xmass)*(589.3/(45.56335/omega))**2
      if(ii==8)then
      refval=0
      inquire(file='.ref',exist=da)
      if(da)then
      open(unit=33,file='.ref')
      read(33,*)refval
      close(33)
      endif
      write(*,142)45.56335/omega,omega*27.21139,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor*xmass/100.0,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor*xmass/100.0,
     .refval,' Na D-line '
      else
      write(*,143)45.56335/omega,omega*27.21139,
     .237.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor,
     .237.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*reffaktor*vorfaktor*xmass/100.0,
     .235.722*(alpha_xx+alpha_yy+alpha_zz)
     ./64604.8*(2.0*137.036/3.0)*vorfaktor*xmass/100.0
      endif

      deallocate(inv_resp)
      enddo
      deallocate(XpY,XmY)
      deallocate(inv_amb)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3333  format(A16,F7.1,A1,F7.1,A1)
 142  format(f6.1,f6.2,5f12.2,a11)
 143  format(f6.1,f6.2,4f12.2)
      end subroutine optrot_velo
