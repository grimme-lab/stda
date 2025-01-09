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
!!!   (|) integrals are computed exactly, for testing purpose only, not optimized !!!!!!!!
      subroutine rrpamat_full(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,
     .                  apb,ambsqr,ca,nao,moci,nprims,epsi)
      use commonlogicals
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,a,b,k,l,ii
      real*4 ek,ej,sdot,ak,ax,de
      real*8 :: ca(nao*moci)
      real*4, allocatable :: integral2(:,:,:,:)
      integer :: nao, moci, nprims
      integer*8 :: lin8
      real*4 :: start_time,end_time,start,finished
      integer :: m,n,o,p
      real*4, allocatable :: iatemp1(:,:),itemp(:,:,:)
      real*4, allocatable :: iatemp2(:,:)
      real*4, allocatable :: iajtemp1(:),iajtemp2(:)
      real*4, allocatable :: K_iajb(:,:,:,:)
      real*4, allocatable :: J_ijab(:,:,:,:)
      real*8 :: ddot
      integer :: canon
      real*8 :: wtime
      integer :: ino,nno,inv,nnv
      integer*8 :: counter
      real*8 :: epsi(moci)

      start_time=wtime()
!       call cpu_time(start)
! Compute (  |  ) integrals with lbcint
      ! reduce the mo range to match the configuration space
      ino=minval(iconf(1:nci,1))
      nno=maxval(iconf(1:nci,1))
      inv=minval(iconf(1:nci,2))-no
      nnv=maxval(iconf(1:nci,2))-no
!       write(*,*)'MO integrals computed for occ.',ino,'to',nno
!       write(*,*)'and for unocc.',inv,'to',nnv

      allocate(integral2(nno+1-ino,nao,nao,nao),stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for 2-e integrals'

      call two_elec_int_i(ncent,nprims,nao,integral2,
     .ino,nno,ca,moci)

!       call cpu_time(finished)
      end_time=wtime()
!       print '("cpu time = ",f12.2," minutes.")'
!      .      ,(finished-start)/60.0
      print '("time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      start_time=wtime()
      call cpu_time(start)


      allocate(itemp(nao,nao,nao))
      allocate(iatemp1(nao,nao))
      allocate(iatemp2(nao,nao))
      allocate(iajtemp1(nao))
      allocate(iajtemp2(nao))



      if(abs(dax).lt.1.0d-6) then
      allocate(K_iajb(ino:nno,inv:nnv,ino:nno,inv:nnv))
!$omp parallel private(ii,p,j,a,b,itemp,
!$omp& iatemp1)
!$omp& shared(ca,integral2,nao,no)
!$omp do
      Do ii=ino, nno
      !(i.|..)
      itemp(:,:,:)=integral2(ii,:,:,:)

      Do a=inv, nnv
      iatemp1=0.0
      Do p=1,nao
      !(ia|..)
      call sgemv('T',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp1(:,p),1)

      Do j=ino, nno
      iajtemp1=0.0
      !(ia|j.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp1,1)
      enddo

      Do b=inv, nnv
      !(ia|jb)
      K_iajb(ii,a,j,b)=sdot(nao,iajtemp1,1,
     .           real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)

      enddo
      enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      else

      allocate(K_iajb(ino:nno,inv:nnv,ino:nno,inv:nnv))
      allocate(J_ijab(ino:nno,ino:nno,inv:nnv,inv:nnv))

!$omp parallel private(ii,p,j,a,b,itemp,
!$omp& iatemp1,iatemp2,iajtemp1,iajtemp2)
!$omp& shared(ca,integral2,nao,no)
!$omp do
      Do ii=ino, nno
      !(i.|..)
      itemp(:,:,:)=integral2(ii,:,:,:) ! this is increasing the speed drastically :)

      Do a=inv, nnv
      iatemp1=0.0
      iatemp2=0.0
      Do p=1,nao
      !(ia|..)
      call sgemv('T',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp1(:,p),1)
      !(i.|a.)
      call sgemv('N',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp2(:,p),1)
      enddo

      Do j=ino, nno
      iajtemp1=0.0
      iajtemp2=0.0
      !(ia|j.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp1,1)
      !(ij|a.)
      call sgemv('T',nao,nao,1.0,iatemp2,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp2,1)

      Do b=inv, nnv
      !(ia|jb)     Since (ia|jb)=(jb|ia), I could decrease memory by collapsing the array... but with respect to integral
      K_iajb(ii,a,j,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      !(ij|ab)
      J_ijab(ii,j,a,b)=sdot(nao,iajtemp2,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      enddo
      enddo
      enddo
      enddo
!$omp end do
!$omp end parallel
      endif

      deallocate(itemp)
      deallocate(iatemp1)
      deallocate(iatemp2)
      deallocate(iajtemp1)
      deallocate(iajtemp2)
      deallocate(integral2)

      write(*,*)'MO (  |  ) computed'
!       call cpu_time(finished)
      end_time=wtime()
!       print '("cpu time = ",f12.2," minutes.")'
!      .      ,(finished-start)/60.0
      print '("time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0


      ak=real(dak)
      ax=real(dax)
! calculate A+B and A-B
      apb=0.0e0
      ambsqr=0.0e0

! if ax=0, A-B is diagonal and its off-diagonal elements do not need to be calculated
      if(abs(dax).lt.1.0d-6) then
        ij=0
!$omp parallel private(ij,i,j,ek,de)
!$omp do
        do i=1,nci
           do j=1,i-1
              ij=lin(i,j)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(j,1),iconf(j,2)-no) ! ek = (ia|jb)
              apb(ij)=2.0*ak*ek
              ambsqr(ij)=0.0
           enddo ! j
           de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           ij=lin(i,i)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(i,1),iconf(i,2)-no)
           if(aresp.or.resp.or.optrota) then
           ambsqr(ij)=de ! diagonal element of (A-B)
           else
           ambsqr(ij)=sqrt(de) ! diagonal element of (A-B)^0.5
           endif
           apb(ij)=de+ak*ek*2.0      ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel

        deallocate(K_iajb)
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) ambsqr
        close(53)

      else

        ij=0
        ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,ek,ej)
!$omp do
        do i=1,nci
           do j=1,i-1
              ij=lin(i,j)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(j,1),iconf(j,2)-no) ! ek = (ia|jb)
        ej=J_ijab(iconf(i,1),iconf(j,1),iconf(i,2)-no,iconf(j,2)-no)*ax !  ej = (ij|ab)
              ambsqr(ij)=2.0*ak*ek
        ek=K_iajb(iconf(i,1),iconf(j,2)-no,iconf(j,1),iconf(i,2)-no) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ambsqr(ij)=ambsqr(ij)-ax*ek-ej
              apb(ij)=ax*ek-ej
           enddo ! j
           ij=lin(i,i)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(i,1),iconf(i,2)-no)
        ej=J_ijab(iconf(i,1),iconf(i,1),iconf(i,2)-no,iconf(i,2)-no)*ax
        de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           apb(ij)=de+ax*ek  -ej  ! diagonal element of A-B
           ambsqr(ij)=de-ax*ek+2.0*ak*ek -ej! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel

      deallocate(K_iajb)
      deallocate(J_ijab)


!        call prmat4(6,apb,nci,0,'A-B')
!        call prmat4(6,ambsqr,nci,0,'A+B')


      if(aresp.or.resp.or.optrota) then
        write(*,*) ' calculating (A-B)^0.5 not necessary...'
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) apb
        close(53)
        apb=ambsqr
      else
        open(unit=52,file='apbmat',form='unformatted',status='replace')
        write(52) ambsqr
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) apb
        write(*,*) ' calculating (A-B)^0.5 ...'
        write(*,'('' estimated time (min) '',f8.2)')
     .            float(nci)**2*float(nci)/4.d+8/60.
        call smatpow(nci,apb) ! calculate sqrt(a-b), termed ambsqr

        ambsqr=apb

        rewind(52)
        read(52) apb
        close(52,status='delete')
        close(53)
      endif
      endif ! GGA/hybrid case


      return

      end subroutine rrpamat_full

        function wtime ( )
        implicit none

        integer ( kind = 4 ) clock_max
        integer ( kind = 4 ) clock_rate
        integer ( kind = 4 ) clock_reading
        real ( kind = 8 ) wtime

        call system_clock ( clock_reading, clock_rate, clock_max )

        wtime = real ( clock_reading, kind = 8 )
     .        / real ( clock_rate, kind = 8 )

        return
        end


      subroutine two_elec_int_i(ncent,nprims,nbf,integral,
     .ino,nno,ca,moci)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l
      integer,target :: n,m,o,p,mm,nn
      integer,pointer :: ddj,ddk,ddl,dddl
      logical :: extra_step
      integer :: ncent,nprims,nbf,canon
      integer :: ino,nno,moci
      real*8 :: ca(nbf*moci)
      integer :: shls(4)
      integer, target :: di,dj,dk,dl
      integer :: orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:,:,:)
      real*4 :: integral(nno+1-ino,nbf,nbf,nbf)
      double precision :: norm_cart(1:10,1:10),thresh

      integer,external :: CINTcgto_cart

      integer*8 :: opt
      ! enforce not RSH integrals, putting mu to zero
      env(9)=0.0d0

      ! see Theoret. Chim. Acta 33 1-6 (1974) Diercksen
      ! and Methods in computational chemistry Vol 1, 1987, page 257
      ! It uses integral symmetries
      ! number of integrals nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2
      ! including vanishing ones
      ! integral canonical index for (ij|kl):
      ! i>=j k>=l
      ! ij=i(i-1)/2+j   kl=k(k-1)/2+l
      ! ij>=kl
      ! canonical_index=ij(ij-1)/2+kl

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
      norm_cart=1.0
      norm_cart(6,1)=1.0
      norm_cart(6,2)=dsqrt(3.0d0)
      norm_cart(6,3)=dsqrt(3.0d0)
      norm_cart(6,4)=1.0
      norm_cart(6,5)=dsqrt(3.0d0)
      norm_cart(6,6)=1.0
      norm_cart(10,1)=1.0
      norm_cart(10,2)=dsqrt(5.0d0)
      norm_cart(10,3)=dsqrt(5.0d0)
      norm_cart(10,4)=dsqrt(5.0d0)
      norm_cart(10,5)=dsqrt(15.0d0)
      norm_cart(10,6)=dsqrt(5.0d0)
      norm_cart(10,7)=1.0
      norm_cart(10,8)=dsqrt(5.0d0)
      norm_cart(10,9)=dsqrt(5.0d0)
      norm_cart(10,10)=1.0
      endif
      integral=0.0d0
      thresh=1.0d-7

      write(*,*)'number of unique AO integrals',
     .int(nbf,8)*(int(nbf,8)+1)/2*
     .       (int(nbf,8)*(int(nbf,8)+1)/2+1)/2
      write(*,*)'number of (i.|..) that will be stored temporarily',
     .int((nno+1-ino),8)*int(nbf,8)*int(nbf,8)*int(nbf,8)
      !compute (ij|kl)
      call cint2e_cart_optimizer(opt,atm,ncent,bas,nbas,env)
!$omp parallel private(i,j,k,l,m,n,o,p,shls,di,dj,dk,dl,buf,
!$omp&                 ddj,ddk,ddl,mm,extra_step)
!$omp& reduction(+:integral)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)
      Do j=1,i
      shls(2)=j-1 ; dj=di_all(j)

      Do k=1,i-1
      shls(3)=k-1 ; dk=di_all(k)

      Do l=1,k
      shls(4)=l-1 ; dl=di_all(l)
      allocate(buf(di,dj,dk,dl))
      call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
      ! select case for writing into integral
      nullify(ddj,ddk,ddl)

      if(i==j)then
      ddj=>m
      else
      ddj=>dj
      endif

      ddk=>dk

      if(k==l)then
      ddl=>o
      else
      ddl=>dl
      endif

      Do m=1,di
      Do n=1,ddj
      Do o=1,ddk
      Do p=1,ddl
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ino,nno)
      endif
      enddo
      enddo
      enddo
      enddo

      deallocate(buf)

      enddo!l

      enddo!k

      k=i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      shls(3)=k-1 ; dk=di_all(k)
      Do l=1,j
      shls(4)=l-1 ; dl=di_all(l)

      allocate(buf(di,dj,dk,dl))
      call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
      ! select case for writing into integral
      extra_step=.false.
      nullify(ddj,ddk,ddl)

      if(i==j)then
      ddj=>m
      else
      ddj=>dj
      endif

      if(k==l)then
      ddk=>mm
      ddl=>o
      extra_step=.true.
      else

      if(j==l)then
      ddk=>mm
      ddl=>dl
      extra_step=.true.
      else
      ddk=>dk
      ddl=>dl
      endif
      endif


      Do m=1,di
      mm=m-1
      Do n=1,ddj
      Do o=1,ddk
      Do p=1,ddl
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ino,nno)
      endif
      enddo
      enddo
      if(extra_step)then
      o=m
      Do p=1,n
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ino,nno)
      endif
      enddo
      endif
      enddo
      enddo

      deallocate(buf)
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      enddo
!$omp end do
!$omp end parallel

      call CINTdel_optimizer(opt)

      write(*,*)'(i.|..) integrals done'
      end subroutine two_elec_int_i

      subroutine select_ijkl(i,j,k,l,ca,nao,moci,itemp,value,ino,nno)
       implicit none
       integer :: i,j,k,l,p
       integer :: ino,nno,moci,nao
       real*4 :: itemp((nno+1-ino),nao,nao,nao)
       real*4 :: value
       real*8 :: ca(nao*moci)


      if(i==j.and.i==k.and.i==l)then!iiii
      Do p=ino,nno
      itemp(p,i,i,i)=itemp(p,i,i,i)+value*ca(i+(p-1)*nao)
      enddo
      elseif(i==j.and.i==k.and.k/=l)then!iiil
      Do p=ino,nno
      itemp(p,i,i,l)=itemp(p,i,i,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,l,i)=itemp(p,i,l,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,l,i,i)=itemp(p,l,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,i,i)=itemp(p,i,i,i)+value*real(ca(l+(p-1)*nao),4)
      enddo
      elseif(i==j.and.i/=k.and.i==l)then!iiki
      Do p=ino,nno
      itemp(p,i,i,k)=itemp(p,i,i,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,k,i)=itemp(p,i,k,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,k,i,i)=itemp(p,k,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,i,i)=itemp(p,i,i,i)+value*real(ca(k+(p-1)*nao),4)
      enddo
      elseif(i/=j.and.i==k.and.i==l)then!ijii
      Do p=ino,nno
      itemp(p,i,i,j)=itemp(p,i,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,j,i)=itemp(p,i,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,j,i,i)=itemp(p,j,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,i,i)=itemp(p,i,i,i)+value*real(ca(j+(p-1)*nao),4)
      enddo
      elseif(i/=j.and.j==k.and.j==l)then!ijjj
      Do p=ino,nno
      itemp(p,j,j,i)=itemp(p,j,j,i)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,j,i,j)=itemp(p,j,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,i,j,j)=itemp(p,i,j,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,j,j,j)=itemp(p,j,j,j)+value*real(ca(i+(p-1)*nao),4)
      enddo
      elseif(i==j.and.i/=k.and.k==l)then!iikk
      Do p=ino,nno
      itemp(p,i,k,k)=itemp(p,i,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,k,i,i)=itemp(p,k,i,i)+value*real(ca(k+(p-1)*nao),4)
      enddo
      elseif(i==j.and.i/=k.and.k/=l)then!iikl
      Do p=ino,nno
      itemp(p,i,k,l)=itemp(p,i,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,l,k)=itemp(p,i,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,l,i,i)=itemp(p,l,i,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(p,k,i,i)=itemp(p,k,i,i)+value*real(ca(l+(p-1)*nao),4)
      enddo
      elseif(i/=j.and.i/=k.and.k==l)then!ijkk
      Do p=ino,nno
      itemp(p,j,k,k)=itemp(p,j,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,k,k)=itemp(p,i,k,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,k,i,j)=itemp(p,k,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(p,k,j,i)=itemp(p,k,j,i)+value*real(ca(k+(p-1)*nao),4)
      enddo
      elseif(i/=j.and.i==k.and.j==l)then!ijij or jiji
      Do p=ino,nno
      itemp(p,j,i,j)=itemp(p,j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,i,j)=itemp(p,i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,j,j,i)=itemp(p,j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,j,i)=itemp(p,i,j,i)+value*real(ca(j+(p-1)*nao),4)
      enddo
      elseif(i/=j.and.i==l.and.j==k)then!ijji or jiij
      Do p=ino,nno
      itemp(p,j,i,j)=itemp(p,j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,i,j)=itemp(p,i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,j,j,i)=itemp(p,j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,j,i)=itemp(p,i,j,i)+value*real(ca(j+(p-1)*nao),4)
      enddo
      else!ijkl ijil
      Do p=ino,nno
      itemp(p,j,k,l)=itemp(p,j,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,j,l,k)=itemp(p,j,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(p,i,k,l)=itemp(p,i,k,l)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,i,l,k)=itemp(p,i,l,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(p,l,i,j)=itemp(p,l,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(p,l,j,i)=itemp(p,l,j,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(p,k,i,j)=itemp(p,k,i,j)+value*real(ca(l+(p-1)*nao),4)
      itemp(p,k,j,i)=itemp(p,k,j,i)+value*real(ca(l+(p-1)*nao),4)
      enddo
      endif


      end subroutine select_ijkl

      subroutine rrpamat_full_direct(nci,ncent,no,nv,mxcnf,iconf,
     .                  dak,dax,ed,apb,ambsqr,ca,nao,moci,nprims,epsi)
      use commonlogicals
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,a,b,k,l,ii
      real*4 ek,ej,sdot,ak,ax,de
      real*8 :: ca(nao*moci)
      integer :: nao, moci, nprims
      integer*8 :: lin8
      real*4 :: start_time,end_time,start,finished
      integer :: m,n,o,p
      real*4, allocatable :: iatemp1(:,:),itemp(:,:,:)
      real*4, allocatable :: iajtemp1(:)
      real*4, allocatable :: K_iajb(:,:,:,:)
      real*4, allocatable :: J_ijab(:,:,:,:)
      real*8 :: ddot
      integer :: canon
      real*8 :: wtime
      integer :: ino,nno,inv,nnv,step
      integer*8 :: counter
      real*8 :: epsi(moci)

      start_time=wtime()
!       call cpu_time(start)
! Compute (  |  ) integrals with lbcint
      ! reduce the mo range to match the configuration space
      ino=minval(iconf(1:nci,1))
      nno=maxval(iconf(1:nci,1))
      inv=minval(iconf(1:nci,2))-no
      nnv=maxval(iconf(1:nci,2))-no
      write(*,*)'MO integrals computed from occ.',ino,'to',nno
      write(*,*)'and from unocc.',inv,'to',nnv


      allocate(itemp(nao,nao,nao), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for itemp'
      allocate(iatemp1(nao,nao), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for iatemp1'
      allocate(iajtemp1(nao), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for iajtemp1'


      if(abs(dax).lt.1.0d-6) then
      allocate(K_iajb(ino:nno,inv:nnv,ino:nno,inv:nnv), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for K_iajb'
!$omp parallel private(i,p,j,a,b,itemp,
!$omp& iatemp1)
!$omp& shared(ca,nao,no)
!$omp do
      Do i=ino, nno
      !(i.|..)
      call two_elec_int_i_direct(ncent,nprims,nao,itemp,
     .ca,moci,i)

      Do a=inv, nnv
      iatemp1=0.0
      Do p=1,nao
      !(ia|..)
      call sgemv('T',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp1(:,p),1)
      enddo

      Do j=ino, i-1
      iajtemp1=0.0
      !(ia|j.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp1,1)

      Do b=inv, nnv
      !(ia|jb)
      K_iajb(i,a,j,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      K_iajb(j,b,i,a)=K_iajb(i,a,j,b)
      enddo
      enddo

      iajtemp1=0.0
      !(ia|i.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(i-1)*nao:nao+(i-1)*nao),4),1,0.0,iajtemp1,1)
      Do b=inv, a-1
      !(ia|ib)
      K_iajb(i,a,i,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      K_iajb(i,b,i,a)=K_iajb(i,a,i,b)
      enddo
      !(ia|ia)
      K_iajb(i,a,i,a)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1)
      enddo
      enddo
!$omp end do
!$omp end parallel
      else

      allocate(K_iajb(ino:nno,inv:nnv,ino:nno,inv:nnv), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for K_iajb'
      allocate(J_ijab(ino:nno,ino:nno,inv:nnv,inv:nnv), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for J_ijab'

!$omp parallel private(i,p,j,a,b,itemp,
!$omp& iatemp1,iajtemp1)
!$omp& shared(ca,nao,no)
!$omp do
      Do i=ino, nno
      !(i.|..)
      call two_elec_int_i_direct(ncent,nprims,nao,itemp,
     .ca,moci,i)

      Do a=inv, nnv

      iatemp1=0.0
      Do p=1,nao
      !(ia|..)
      call sgemv('T',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp1(:,p),1)
      enddo

      Do j=ino, i-1
      iajtemp1=0.0
      !(ia|j.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp1,1)

      Do b=inv, nnv
      !(ia|jb)
      K_iajb(i,a,j,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      K_iajb(j,b,i,a)=K_iajb(i,a,j,b)
      enddo
      enddo

      iajtemp1=0.0
      !(ia|i.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(i-1)*nao:nao+(i-1)*nao),4),1,0.0,iajtemp1,1)
      Do b=inv, a-1
      !(ia|ib)
      K_iajb(i,a,i,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      K_iajb(i,b,i,a)=K_iajb(i,a,i,b)
      enddo
      !(ia|ia)
      K_iajb(i,a,i,a)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1)



      iatemp1=0.0
      Do p=1,nao
      !(i.|a.)
      call sgemv('N',nao,nao,1.0,itemp(:,:,p),nao,
     .real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1,0.0,iatemp1(:,p),1)
      enddo

      Do j=ino, i-1
      iajtemp1=0.0
      !(ij|a.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(j-1)*nao:nao+(j-1)*nao),4),1,0.0,iajtemp1,1)

      Do b=inv, nnv
      !(ij|ab)
      J_ijab(i,j,a,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      J_ijab(j,i,b,a)=J_ijab(i,j,a,b)
      enddo
      enddo

      iajtemp1=0.0
      !(ii|a.)
      call sgemv('T',nao,nao,1.0,iatemp1,nao,
     .real(ca(1+(i-1)*nao:nao+(i-1)*nao),4),1,0.0,iajtemp1,1)

      Do b=inv, a-1
      !(ii|ab)
      J_ijab(i,i,a,b)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(b+no-1)*nao:nao+(b+no-1)*nao),4),1)
      J_ijab(i,i,b,a)=J_ijab(i,i,a,b)
      enddo
      !(ii|aa)
      J_ijab(i,i,a,a)=sdot(nao,iajtemp1,1,
     .                real(ca(1+(a+no-1)*nao:nao+(a+no-1)*nao),4),1)

      enddo
      enddo
!$omp end do
!$omp end parallel
      endif

      deallocate(itemp)
      deallocate(iatemp1)
      deallocate(iajtemp1)

      write(*,*)'MO (  |  ) computed'
!       call cpu_time(finished)
      end_time=wtime()
!       print '("cpu time = ",f12.2," minutes.")'
!      .      ,(finished-start)/60.0
      print '("time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0


      ak=real(dak)
      ax=real(dax)
! calculate A+B and A-B
      apb=0.0e0
      ambsqr=0.0e0

! if ax=0, A-B is diagonal and its off-diagonal elements do not need to be calculated
      if(abs(dax).lt.1.0d-6) then
        ij=0
!$omp parallel private(ij,i,j,ek,de)
!$omp do
        do i=1,nci
           do j=1,i-1
              ij=lin(i,j)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(j,1),iconf(j,2)-no) ! ek = (ia|jb)
              apb(ij)=2.0*ak*ek
              ambsqr(ij)=0.0
           enddo ! j
           de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           ij=lin(i,i)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(i,1),iconf(i,2)-no)
           if(aresp.or.resp.or.optrota) then
           ambsqr(ij)=de ! diagonal element of (A-B)
           else
           ambsqr(ij)=sqrt(de) ! diagonal element of (A-B)^0.5
           endif
           apb(ij)=de+ak*ek*2.0       ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel

        deallocate(K_iajb)
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) ambsqr
        close(53)

      else

        ij=0
        ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,ek,ej)
!$omp do
        do i=1,nci
           do j=1,i-1
              ij=lin(i,j)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(j,1),iconf(j,2)-no) ! ek = (ia|jb)
        ej=J_ijab(iconf(i,1),iconf(j,1),iconf(i,2)-no,iconf(j,2)-no)*ax !  ej = (ij|ab)
              ambsqr(ij)=2.0*ak*ek
        ek=K_iajb(iconf(i,1),iconf(j,2)-no,iconf(j,1),iconf(i,2)-no) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ambsqr(ij)=ambsqr(ij)-ax*ek-ej
              apb(ij)=ax*ek-ej
           enddo ! j
           ij=lin(i,i)
        ek=k_iajb(iconf(i,1),iconf(i,2)-no,iconf(i,1),iconf(i,2)-no)
        ej=J_ijab(iconf(i,1),iconf(i,1),iconf(i,2)-no,iconf(i,2)-no)*ax
        de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           apb(ij)=de+ax*ek-ej    ! diagonal element of A-B
           ambsqr(ij)=de-ax*ek+2.0*ak*ek-ej ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel

      deallocate(K_iajb)
      deallocate(J_ijab)


!        call prmat4(6,apb,nci,0,'A-B')
!        call prmat4(6,ambsqr,nci,0,'A+B')


      if(aresp.or.resp.or.optrota) then
        write(*,*) ' calculating (A-B)^0.5 not necessary...'
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) apb
        close(53)
        apb=ambsqr
      else
        open(unit=52,file='apbmat',form='unformatted',status='replace')
        write(52) ambsqr
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) apb
        write(*,*) ' calculating (A-B)^0.5 ...'
        write(*,'('' estimated time (min) '',f8.2)')
     .            float(nci)**2*float(nci)/4.d+8/60.
        call smatpow(nci,apb) ! calculate sqrt(a-b), termed ambsqr

        ambsqr=apb

        rewind(52)
        read(52) apb
        close(52,status='delete')
        close(53)
      endif
      endif ! GGA/hybrid case


      return

      end subroutine rrpamat_full_direct

      subroutine two_elec_int_i_direct(ncent,nprims,nbf,integral,
     .ca,moci,ii)
      use stdacommon
      use commonlibcint
      use omp_lib
      implicit none
      integer :: i,j,k,l,ii
      integer,target :: n,m,o,p,mm,nn
      integer,pointer :: ddj,ddk,ddl,dddl
      logical :: extra_step
      integer :: ncent,nprims,nbf,canon
      integer :: moci
      real*8 :: ca(nbf*moci)
      integer :: shls(4)
      integer, target :: di,dj,dk,dl
      integer :: orb_cart(1:10,1:10)
      double precision, allocatable :: buf(:,:,:,:)
      real*4 :: integral(nbf,nbf,nbf)
      double precision :: norm_cart(1:10,1:10),thresh

      integer,external :: CINTcgto_cart

      integer*8 :: opt



      ! see Theoret. Chim. Acta 33 1-6 (1974) Diercksen
      ! and Methods in computational chemistry Vol 1, 1987, page 257
      ! It uses integral symmetries
      ! number of integrals nbf*(nbf+1)/2*(nbf*(nbf+1)/2+1)/2
      ! including vanishing ones
      ! integral canonical index for (ij|kl):
      ! i>=j k>=l
      ! ij=i(i-1)/2+j   kl=k(k-1)/2+l
      ! ij>=kl
      ! canonical_index=ij(ij-1)/2+kl

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
      norm_cart=1.0
      norm_cart(6,1)=1.0
      norm_cart(6,2)=dsqrt(3.0d0)
      norm_cart(6,3)=dsqrt(3.0d0)
      norm_cart(6,4)=1.0
      norm_cart(6,5)=dsqrt(3.0d0)
      norm_cart(6,6)=1.0
      norm_cart(10,1)=1.0
      norm_cart(10,2)=dsqrt(5.0d0)
      norm_cart(10,3)=dsqrt(5.0d0)
      norm_cart(10,4)=dsqrt(5.0d0)
      norm_cart(10,5)=dsqrt(15.0d0)
      norm_cart(10,6)=dsqrt(5.0d0)
      norm_cart(10,7)=1.0
      norm_cart(10,8)=dsqrt(5.0d0)
      norm_cart(10,9)=dsqrt(5.0d0)
      norm_cart(10,10)=1.0
      endif
      integral=0.0d0
      thresh=1.0d-7

      !compute (ij|kl)
      call cint2e_cart_optimizer(opt,atm,ncent,bas,nbas,env)
!$omp parallel private(i,j,k,l,m,n,o,p,shls,di,dj,dk,dl,buf,
!$omp&                 ddj,ddk,ddl,mm,extra_step)
!$omp& reduction(+:integral)
!$omp do
      Do i=1,nbas
      shls(1)=i-1 ; di=di_all(i)
      Do j=1,i
      shls(2)=j-1 ; dj=di_all(j)

      Do k=1,i-1
      shls(3)=k-1 ; dk=di_all(k)

      Do l=1,k
      shls(4)=l-1 ; dl=di_all(l)
      allocate(buf(di,dj,dk,dl))
      !!!! this is the step that cost too much
      call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
      ! select case for writing into integral
      nullify(ddj,ddk,ddl)

      if(i==j)then
      ddj=>m
      else
      ddj=>dj
      endif

      ddk=>dk

      if(k==l)then
      ddl=>o
      else
      ddl=>dl
      endif

      Do m=1,di
      Do n=1,ddj
      Do o=1,ddk
      Do p=1,ddl
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl_directA(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ii)
      endif
      enddo
      enddo
      enddo
      enddo

      deallocate(buf)

      enddo!l

      enddo!k

      k=i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      shls(3)=k-1 ; dk=di_all(k)
      Do l=1,j
      shls(4)=l-1 ; dl=di_all(l)

      allocate(buf(di,dj,dk,dl))
      call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
      ! select case for writing into integral
      extra_step=.false.
      nullify(ddj,ddk,ddl)

      if(i==j)then
      ddj=>m
      else
      ddj=>dj
      endif

      if(k==l)then
      ddk=>mm
      ddl=>o
      extra_step=.true.
      else

      if(j==l)then
      ddk=>mm
      ddl=>dl
      extra_step=.true.
      else
      ddk=>dk
      ddl=>dl
      endif
      endif


      Do m=1,di
      mm=m-1
      Do n=1,ddj
      Do o=1,ddk
      Do p=1,ddl
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl_direct(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ii)
      endif
      enddo
      enddo
      if(extra_step)then
      o=m
      Do p=1,n
      if(dabs(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p))>thresh)then
      call select_ijkl_directB(sum_di(i)+orb_cart(di,m)
     .,sum_di(j)+orb_cart(dj,n)
     .,sum_di(k)+orb_cart(dk,o)
     .,sum_di(l)+orb_cart(dl,p),ca,nbf,moci,integral,real(buf(m,n,o,p)
     .*norm_cart(di,m)
     .*norm_cart(dj,n)
     .*norm_cart(dk,o)
     .*norm_cart(dl,p),4),ii)
      endif
      enddo
      endif
      enddo
      enddo

      deallocate(buf)
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      enddo
!$omp end do
!$omp end parallel

      call CINTdel_optimizer(opt)


      end subroutine two_elec_int_i_direct

      subroutine select_ijkl_direct(i,j,k,l,ca,nao,moci,itemp,value,p)
       implicit none
       integer :: i,j,k,l,p
       integer :: moci,nao
       real*4 :: itemp(nao,nao,nao)
       real*4 :: value
       real*8 :: ca(nao*moci)


      if(i==j.and.i==k.and.i==l)then!iiii

      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(i+(p-1)*nao),4)

      elseif(i==j.and.i==k.and.k/=l)then!iiil

      itemp(i,i,l)=itemp(i,i,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,l,i)=itemp(i,l,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(l,i,i)=itemp(l,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(l+(p-1)*nao),4)

      elseif(i==j.and.i/=k.and.i==l)then!iiki

      itemp(i,i,k)=itemp(i,i,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,i)=itemp(i,k,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i/=j.and.i==k.and.i==l)then!ijii

      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(j,i,i)=itemp(j,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(j+(p-1)*nao),4)

      elseif(i/=j.and.j==k.and.j==l)then!ijjj

      itemp(j,j,i)=itemp(j,j,i)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,i,j)=itemp(j,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(i,j,j)=itemp(i,j,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,j,j)=itemp(j,j,j)+value*real(ca(i+(p-1)*nao),4)

      elseif(i==j.and.i/=k.and.k==l)then!iikk

      itemp(i,k,k)=itemp(i,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i==j.and.i/=k.and.k/=l)then!iikl

      itemp(i,k,l)=itemp(i,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,l,k)=itemp(i,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(l,i,i)=itemp(l,i,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(l+(p-1)*nao),4)

      elseif(i/=j.and.i/=k.and.k==l)then!ijkk

      itemp(j,k,k)=itemp(j,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,k)=itemp(i,k,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(k,i,j)=itemp(k,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,j,i)=itemp(k,j,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i/=j.and.i==k.and.j==l)then!ijij or jiji

      itemp(j,i,j)=itemp(j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,j,i)=itemp(j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(j+(p-1)*nao),4)

      elseif(i/=j.and.i==l.and.j==k)then!ijji or jiij

      itemp(j,i,j)=itemp(j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,j,i)=itemp(j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(j+(p-1)*nao),4)

      else!ijkl ijil

      itemp(j,k,l)=itemp(j,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(j,l,k)=itemp(j,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,l)=itemp(i,k,l)+value*real(ca(j+(p-1)*nao),4)
      itemp(i,l,k)=itemp(i,l,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(l,i,j)=itemp(l,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(l,j,i)=itemp(l,j,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,i,j)=itemp(k,i,j)+value*real(ca(l+(p-1)*nao),4)
      itemp(k,j,i)=itemp(k,j,i)+value*real(ca(l+(p-1)*nao),4)

      endif


      end subroutine select_ijkl_direct
      subroutine select_ijkl_directA(i,j,k,l,ca,nao,moci,itemp,value,p)
       implicit none
       integer :: i,j,k,l,p
       integer :: moci,nao
       real*4 :: itemp(nao,nao,nao)
       real*4 :: value
       real*8 :: ca(nao*moci)


      if(i==j.and.i/=k.and.i==l)then!iiki

      itemp(i,i,k)=itemp(i,i,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,i)=itemp(i,k,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i==j.and.i/=k.and.k==l)then!iikk

      itemp(i,k,k)=itemp(i,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i==j.and.i/=k.and.k/=l)then!iikl

      itemp(i,k,l)=itemp(i,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,l,k)=itemp(i,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(l,i,i)=itemp(l,i,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,i,i)=itemp(k,i,i)+value*real(ca(l+(p-1)*nao),4)

      elseif(i/=j.and.i/=k.and.k==l)then!ijkk

      itemp(j,k,k)=itemp(j,k,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,k)=itemp(i,k,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(k,i,j)=itemp(k,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,j,i)=itemp(k,j,i)+value*real(ca(k+(p-1)*nao),4)

      elseif(i/=j.and.i==l.and.j==k)then!ijji or jiij

      itemp(j,i,j)=itemp(j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,j,i)=itemp(j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(j+(p-1)*nao),4)

      else!ijkl ijil

      itemp(j,k,l)=itemp(j,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(j,l,k)=itemp(j,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,l)=itemp(i,k,l)+value*real(ca(j+(p-1)*nao),4)
      itemp(i,l,k)=itemp(i,l,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(l,i,j)=itemp(l,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(l,j,i)=itemp(l,j,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,i,j)=itemp(k,i,j)+value*real(ca(l+(p-1)*nao),4)
      itemp(k,j,i)=itemp(k,j,i)+value*real(ca(l+(p-1)*nao),4)

      endif


      end subroutine select_ijkl_directA

       subroutine select_ijkl_directB(i,j,k,l,ca,nao,moci,itemp,value,p)
       implicit none
       integer :: i,j,k,l,p
       integer :: moci,nao
       real*4 :: itemp(nao,nao,nao)
       real*4 :: value
       real*8 :: ca(nao*moci)


      if(i==j.and.i==k.and.i==l)then!iiii

      itemp(i,i,i)=itemp(i,i,i)+value*ca(i+(p-1)*nao)

      elseif(i==j.and.i==k.and.k/=l)then!iiil

      itemp(i,i,l)=itemp(i,i,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,l,i)=itemp(i,l,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(l,i,i)=itemp(l,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(l+(p-1)*nao),4)

      elseif(i/=j.and.i==k.and.i==l)then!ijii

      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(j,i,i)=itemp(j,i,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,i)=itemp(i,i,i)+value*real(ca(j+(p-1)*nao),4)

      elseif(i/=j.and.i==k.and.j==l)then!ijij or jiji

      itemp(j,i,j)=itemp(j,i,j)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,i,j)=itemp(i,i,j)+value*real(ca(j+(p-1)*nao),4)
      itemp(j,j,i)=itemp(j,j,i)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,j,i)=itemp(i,j,i)+value*real(ca(j+(p-1)*nao),4)

      else!ijkl ijil

      itemp(j,k,l)=itemp(j,k,l)+value*real(ca(i+(p-1)*nao),4)
      itemp(j,l,k)=itemp(j,l,k)+value*real(ca(i+(p-1)*nao),4)
      itemp(i,k,l)=itemp(i,k,l)+value*real(ca(j+(p-1)*nao),4)
      itemp(i,l,k)=itemp(i,l,k)+value*real(ca(j+(p-1)*nao),4)
      itemp(l,i,j)=itemp(l,i,j)+value*real(ca(k+(p-1)*nao),4)
      itemp(l,j,i)=itemp(l,j,i)+value*real(ca(k+(p-1)*nao),4)
      itemp(k,i,j)=itemp(k,i,j)+value*real(ca(l+(p-1)*nao),4)
      itemp(k,j,i)=itemp(k,j,i)+value*real(ca(l+(p-1)*nao),4)

      endif


      end subroutine select_ijkl_directB
