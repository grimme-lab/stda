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
subroutine xstd_rpamat(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),P_ij(:,:),AABB_integral(:,:)

real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv




if(abs(dax).lt.1.0d-6) then

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

else

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao),P_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ij,AABB_integral)

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

endif

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(dax)
! calculate A+B and A-B
apb=0.0e0
ambsqr=0.0e0

! if ax=0, A-B is diagonal and its off-diagonal elements do not need to be calculated
if(abs(dax).lt.1.0d-6) then
  ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,ek,de)
!$omp do
do i=1,nci
   io=iconf(i,1)
   iv=iconf(i,2)
   iiv=iv-no
   iwrk=(io-1)*nv + iiv
   do j=1,i-1
      ij=lin(i,j)
      jo=iconf(j,1)
      jv=iconf(j,2)
      jjv=jv-no
      jwrk=(jo-1)*nv+jjv
      ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
      apb(ij)=2.0*ak*ek
      ambsqr(ij)=0.0
   enddo ! j
   de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
   ij=lin(i,i)
   ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
   if(aresp.or.resp.or.optrota) then
   ambsqr(ij)=de ! diagonal element of (A-B)
   else
   ambsqr(ij)=sqrt(de) ! diagonal element of (A-B)^0.5
   endif
   apb(ij)=de+ak*ek*2.0       ! diagonal element of A+B
enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia)
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) ambsqr
        close(53)
else
  ij=0
  ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,iwrk2,ek,ej)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     iiv=iv-no
     iwrk=(io-1)*nv + iiv
     do j=1,i-1
        ij=lin(i,j)
        jo=iconf(j,1)
        jv=iconf(j,2)
        jjv=jv-no
        jwrk=(jo-1)*nv+jjv
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        ambsqr(ij)=2.0*ak*ek
        iwrk2=(io-1)*nv+jjv
        jwrk=(jo-1)*nv+iiv
        ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
ek=sdot(nao,P_ia(io,jv,:),1,Q_ia(jo,iv,:),1)
        ambsqr(ij)=ambsqr(ij)-ax*ek-ej
        apb(ij)=ax*ek-ej
     enddo ! j
     ij=lin(i,i)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           apb(ij)=de+ax*ek  -ej  ! diagonal element of A-B
           ambsqr(ij)=de-ax*ek+2.0*ak*ek -ej! diagonal element of A+B
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,P_ij,Q_ab)


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
  write(*,'('' estimated time (min) '',f8.2)')float(nci)**2*float(nci)/4.d+8/60.
  call smatpow(nci,apb) ! calculate sqrt(a-b), termed ambsqr

  ambsqr=apb

  rewind(52)
  read(52) apb
  close(52,status='delete')
  close(53)
endif


endif ! GGA/hybrid case

return

end subroutine xstd_rpamat

subroutine Xstda_mat(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: hci(nci,nci)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),P_ij(:,:),AABB_integral(:,:)

real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao),P_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ij,AABB_integral)

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(dax)
! calculate A+B and A-B
hci=0.0e0
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     do j=1,i-1
        jo=iconf(j,1)
        jv=iconf(j,2)
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        hci(i,j)=ak*ek-ej
        hci(j,i)=hci(i,j)
     enddo ! j
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           hci(i,i)=de+ak*ek-ej
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,P_ij,Q_ab)

return

end subroutine Xstda_mat

subroutine two_elec_int(ncent,nbfA,nbfB,integral,start_basA,end_basA,start_basB,end_basB)
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer :: i,j,k,l,start_basA,end_basA,start_basB,end_basB
integer,target :: n,m,o,p,mm,nn
integer,pointer :: ddj,ddk,ddl,dddl
logical :: extra_step
integer :: ncent,nbfA,nbfB,canon


integer :: shls(4),counter
integer, target :: di,dj,dk,dl
integer :: orb_cart(1:10,1:10)
double precision, allocatable :: buf(:,:,:,:)
real*4 :: integral(nbfA,nbfB)
double precision :: norm_cart(1:10,1:10),thresh

integer,external :: CINTcgto_cart

integer*8 :: opt

integer :: sm_di(start_basA:end_basA),sm_dj(start_basB:end_basB)

! enforce not RSH integrals, putting mu to zero
env(9)=0.0d0

!only (aa|aa) to compute

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

integral=0.0
thresh=1.0d-7
write(*,*)"number of AO's (aa|bb)",nbfA*nbfB
!compute (ij|kl)
call cint2e_cart_optimizer(opt,atm,ncent,bas,nbas,env)

sm_di(start_basA)=0
Do i=start_basA,end_basA-1
sm_di(i+1)=sm_di(i)+di_all(i)
enddo
sm_dj(start_basB)=0
Do i=start_basB,end_basB-1
sm_dj(i+1)=sm_dj(i)+di_all(i)
enddo

!$omp parallel private(i,j,m,o,shls,di,dj,dk,dl,buf)
!$omp do
Do i=start_basA,end_basA
shls(1)=i-1 ; di=di_all(i)
shls(2)=i-1 ; dj=di_all(i)
Do j=start_basB,end_basB
shls(3)=j-1 ; dk=di_all(j)
shls(4)=j-1 ; dl=di_all(j)
allocate(buf(di,dj,dk,dl))
call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
Do m=1,di
Do o=1,dk
integral(sm_di(i)+orb_cart(di,m),sm_dj(j)+orb_cart(dk,o))&
&=real(buf(m,m,o,o)*norm_cart(di,m)*norm_cart(di,m)*norm_cart(dk,o)*norm_cart(dk,o),4)
enddo
enddo

deallocate(buf)

enddo
enddo
!$omp end do
!$omp end parallel

call CINTdel_optimizer(opt)
write(*,*)"AO's (aa|bb) done"
end subroutine two_elec_int

subroutine two_elec_int_RSH(ncent,nbfA,nbfB,integral,start_basA,end_basA,start_basB,end_basB)
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer :: i,j,k,l,start_basA,end_basA,start_basB,end_basB
integer,target :: n,m,o,p,mm,nn
integer,pointer :: ddj,ddk,ddl,dddl
logical :: extra_step
integer :: ncent,nbfA,nbfB,canon


integer :: shls(4),counter
integer, target :: di,dj,dk,dl
integer :: orb_cart(1:10,1:10)
double precision, allocatable :: buf(:,:,:,:)
real*4 :: integral(nbfA,nbfB)
double precision :: norm_cart(1:10,1:10),thresh

integer,external :: CINTcgto_cart

integer*8 :: opt

integer :: sm_di(start_basA:end_basA),sm_dj(start_basB:end_basB)

! Range separated hybrid functional

env(9)=RSH

!only (aa|aa) to compute

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

integral=0.0
thresh=1.0d-7
write(*,*)"number of AO's (aa|bb)",nbfA*nbfB
!compute (ij|kl)
call cint2e_cart_optimizer(opt,atm,ncent,bas,nbas,env)

sm_di(start_basA)=0
Do i=start_basA,end_basA-1
sm_di(i+1)=sm_di(i)+di_all(i)
enddo
sm_dj(start_basB)=0
Do i=start_basB,end_basB-1
sm_dj(i+1)=sm_dj(i)+di_all(i)
enddo

!$omp parallel private(i,j,m,o,shls,di,dj,dk,dl,buf)
!$omp do
Do i=start_basA,end_basA
shls(1)=i-1 ; di=di_all(i)
shls(2)=i-1 ; dj=di_all(i)
Do j=start_basB,end_basB
shls(3)=j-1 ; dk=di_all(j)
shls(4)=j-1 ; dl=di_all(j)
allocate(buf(di,dj,dk,dl))
call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
Do m=1,di
Do o=1,dk
integral(sm_di(i)+orb_cart(di,m),sm_dj(j)+orb_cart(dk,o))&
&=real(buf(m,m,o,o)*norm_cart(di,m)*norm_cart(di,m)*norm_cart(dk,o)*norm_cart(dk,o),4)
enddo
enddo

deallocate(buf)

enddo
enddo
!$omp end do
!$omp end parallel

call CINTdel_optimizer(opt)
write(*,*)"AO's (aa|bb) done"
end subroutine two_elec_int_RSH

subroutine two_elec_int_RSH_SR(ncent,nbfA,nbfB,integral,start_basA,end_basA,start_basB,end_basB)
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer :: i,j,k,l,start_basA,end_basA,start_basB,end_basB
integer,target :: n,m,o,p,mm,nn
integer,pointer :: ddj,ddk,ddl,dddl
logical :: extra_step
integer :: ncent,nbfA,nbfB,canon


integer :: shls(4),counter
integer, target :: di,dj,dk,dl
integer :: orb_cart(1:10,1:10)
double precision, allocatable :: buf(:,:,:,:)
real*4 :: integral(nbfA,nbfB)
double precision :: norm_cart(1:10,1:10),thresh

integer,external :: CINTcgto_cart

integer*8 :: opt

integer :: sm_di(start_basA:end_basA),sm_dj(start_basB:end_basB)

! Range separated hybrid functional

env(9)=-RSH2

!only (aa|aa) to compute

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

integral=0.0
thresh=1.0d-7
write(*,*)"number of AO's (aa|bb)",nbfA*nbfB
!compute (ij|kl)
call cint2e_cart_optimizer(opt,atm,ncent,bas,nbas,env)

sm_di(start_basA)=0
Do i=start_basA,end_basA-1
sm_di(i+1)=sm_di(i)+di_all(i)
enddo
sm_dj(start_basB)=0
Do i=start_basB,end_basB-1
sm_dj(i+1)=sm_dj(i)+di_all(i)
enddo

!$omp parallel private(i,j,m,o,shls,di,dj,dk,dl,buf)
!$omp do
Do i=start_basA,end_basA
shls(1)=i-1 ; di=di_all(i)
shls(2)=i-1 ; dj=di_all(i)
Do j=start_basB,end_basB
shls(3)=j-1 ; dk=di_all(j)
shls(4)=j-1 ; dl=di_all(j)
allocate(buf(di,dj,dk,dl))
call cint2e_cart(buf,shls,atm,ncent,bas,nbas,env,opt)
Do m=1,di
Do o=1,dk
integral(sm_di(i)+orb_cart(di,m),sm_dj(j)+orb_cart(dk,o))&
&=real(buf(m,m,o,o)*norm_cart(di,m)*norm_cart(di,m)*norm_cart(dk,o)*norm_cart(dk,o),4)
enddo
enddo

deallocate(buf)

enddo
enddo
!$omp end do
!$omp end parallel

call CINTdel_optimizer(opt)
write(*,*)"AO's (aa|bb) done"
end subroutine two_elec_int_RSH_SR

subroutine xstd_rpamat_RSH(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj,ek_RSH

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),P_ij(:,:),AABB_integral(:,:),P_ia_RSH(:,:,:),P_ij_RSH(:,:)
real*4,allocatable :: P_ia_RSH_SR(:,:,:),P_ij_RSH_SR(:,:)
real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)
allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!RSH integrals short range

allocate(AABB_integral(nao,nao))
call two_elec_int_RSH_SR(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ia_RSH_SR(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia_RSH_SR(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ij_RSH_SR(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH_SR(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

!RSH integrals long range

allocate(AABB_integral(nao,nao))
call two_elec_int_RSH(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ia_RSH(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia_RSH(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ij_RSH(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel


deallocate(Q_ij,AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(RSH_ax,4)
! calculate A+B and A-B
apb=0.0e0
ambsqr=0.0e0


  ij=0
  ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,iwrk2,ek,ej)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     iiv=iv-no
     iwrk=(io-1)*nv + iiv
     do j=1,i-1
        ij=lin(i,j)
        jo=iconf(j,1)
        jv=iconf(j,2)
        jjv=jv-no
        jwrk=(jo-1)*nv+jjv
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij_RSH_SR(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        ambsqr(ij)=2.0*ak*ek
        iwrk2=(io-1)*nv+jjv
        jwrk=(jo-1)*nv+iiv
        ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
ek=ax*sdot(nao,P_ia_RSH_SR(io,jv,:),1,Q_ia(jo,iv,:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,jv,:),1,Q_ia(jo,iv,:),1)
        ambsqr(ij)=ambsqr(ij)-ek-ej
        apb(ij)=ek-ej
     enddo ! j
     ij=lin(i,i)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ek_RSH=ax*sdot(nao,P_ia_RSH_SR(io,iv,:),1,Q_ia(io,iv,:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij_RSH_SR(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           apb(ij)=de+ek_RSH  -ej  ! diagonal element of A-B
           ambsqr(ij)=de-ek_RSH+2.0*ak*ek -ej! diagonal element of A+B
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,Q_ab,P_ia_RSH,P_ij_RSH,P_ia_RSH_SR,P_ij_RSH_SR)


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
  write(*,'('' estimated time (min) '',f8.2)')float(nci)**2*float(nci)/4.d+8/60.
  call smatpow(nci,apb) ! calculate sqrt(a-b), termed ambsqr

  ambsqr=apb

  rewind(52)
  read(52) apb
  close(52,status='delete')
  close(53)
endif




return

end subroutine xstd_rpamat_RSH

subroutine xstd_rpamat_RSH2(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj,ek_RSH

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),P_ij(:,:),AABB_integral(:,:),P_ia_RSH(:,:,:),P_ij_RSH(:,:)
real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

!RSH integrals long range

allocate(AABB_integral(nao,nao))
call two_elec_int_RSH(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ia_RSH(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia_RSH(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ij_RSH(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel


deallocate(Q_ij,AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(RSH_ax,4)
! calculate A+B and A-B
apb=0.0e0
ambsqr=0.0e0


  ij=0
  ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,iwrk2,ek,ej)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     iiv=iv-no
     iwrk=(io-1)*nv + iiv
     do j=1,i-1
        ij=lin(i,j)
        jo=iconf(j,1)
        jv=iconf(j,2)
        jjv=jv-no
        jwrk=(jo-1)*nv+jjv
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        ambsqr(ij)=2.0*ak*ek
        iwrk2=(io-1)*nv+jjv
        jwrk=(jo-1)*nv+iiv
        ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
ek=ax*sdot(nao,P_ia(io,jv,:),1,Q_ia(jo,iv,:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,jv,:),1,Q_ia(jo,iv,:),1)
        ambsqr(ij)=ambsqr(ij)-ek-ej
        apb(ij)=ek-ej
     enddo ! j
     ij=lin(i,i)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ek_RSH=ax*sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           apb(ij)=de+ek_RSH  -ej  ! diagonal element of A-B
           ambsqr(ij)=de-ek_RSH+2.0*ak*ek -ej! diagonal element of A+B
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,Q_ab,P_ia_RSH,P_ij_RSH,P_ij)


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
  write(*,'('' estimated time (min) '',f8.2)')float(nci)**2*float(nci)/4.d+8/60.
  call smatpow(nci,apb) ! calculate sqrt(a-b), termed ambsqr

  ambsqr=apb

  rewind(52)
  read(52) apb
  close(52,status='delete')
  close(53)
endif




return

end subroutine xstd_rpamat_RSH2

subroutine Xstda_mat_RSH(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: hci(nci,nci)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj,ek_RSH

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),AABB_integral(:,:),P_ij_RSH(:,:)
real*4,allocatable :: P_ij_RSH_SR(:,:)

real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel


!RSH integrals short range

allocate(AABB_integral(nao,nao))
call two_elec_int_RSH_SR(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ij_RSH_SR(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH_SR(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

!RSH integrals long range

allocate(AABB_integral(nao,nao))
call two_elec_int_RSH(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ij_RSH(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ij,AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(RSH_ax,4)
! calculate A+B and A-B
hci=0.0e0
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     do j=1,i-1
        jo=iconf(j,1)
        jv=iconf(j,2)
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij_RSH_SR(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        hci(i,j)=ak*ek-ej
        hci(j,i)=hci(i,j)
     enddo ! j
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij_RSH_SR(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           hci(i,i)=de+ak*ek-ej
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,P_ij_RSH,P_ij_RSH_SR,Q_ab)

return

end subroutine Xstda_mat_RSH

subroutine Xstda_mat_RSH2(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alpha,beta,epsi) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dak,dax,ed(mxcnf),alpha,beta
real*4, intent(out) :: hci(nci,nci)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ak,ax,de,alphak,betaj,ek_RSH

real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),P_ij(:,:),AABB_integral(:,:),P_ij_RSH(:,:)

real*4 :: value1,value2,value3,rabx
real*8 :: epsi(moci)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
Do k=ino, i
Do j=1, nao
    Q_ij(lin(i,k),j)=clow(j+(i-1)*nao)*clow(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

allocate(Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inv,nnv
Do l=inv, j
Do i=1,nao
    Q_ab(lin(j,l),i)=clow(i+(j-1)*nao)*clow(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel



!RSH integrals long range

allocate(AABB_integral(nao,nao))

call two_elec_int_RSH(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)

allocate(P_ij_RSH(ino*(ino+1)/2:nno*(nno+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=ino, nno
    Do k=ino, nno

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ij(lin(i,k),:),1,0.0,P_ij_RSH(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ij,AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ak=real(dak)
ax=real(RSH_ax,4)
! calculate A+B and A-B
hci=0.0e0
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
  do i=1,nci
     io=iconf(i,1)
     iv=iconf(i,2)
     do j=1,i-1
        jo=iconf(j,1)
        jv=iconf(j,2)
        ! ek = (ia|jb)
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ij(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,jo),:),1,Q_ab(lin(iv,jv),:),1)
        hci(i,j)=ak*ek-ej
        hci(j,i)=hci(i,j)
     enddo ! j
ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
ej=ax*sdot(nao,P_ij(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)&
&+real(RSH_beta,4)*sdot(nao,P_ij_RSH(lin(io,io),:),1,Q_ab(lin(iv,iv),:),1)
de=real(epsi(iconf(i,2))-epsi(iconf(i,1)),4)
           hci(i,i)=de+ak*ek-ej
  enddo ! i
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,P_ij_RSH,P_ij,Q_ab)

return

end subroutine Xstda_mat_RSH2


subroutine Xsrtdacorr(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,clow,nao,moci)
use omp_lib
use commonlogicals
use stdacommon
use commonlibcint
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: dak,dax,ed(mxcnf)
real*4, allocatable :: bmat(:)
integer i,j,k,io,iv,jo,jv,ierr,iiv,jjv,iwrk,jwrk,lin,ij
real*4 ek,ej,sdot,ak,ax,de,fact


real*4 :: start_time,end_time,start
integer :: nao,moci
integer :: ino,nno,inv,nnv
real*8  :: clow(nao*moci)
real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),AABB_integral(:,:)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0

ij=nci
ij=ij*(ij+1)/2
allocate(bmat(ij), stat=ierr)
if(ierr.ne.0)stop 'allocation for bmat crashed'
ak=real(dak)
ax=real(dax)
! calculate 0.5*B
bmat=0.0e0
fact=0.50d0 ! this is the scaling of the B-contribution
open(unit=52,file='bmat',form='unformatted',status='replace')
ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,ek,ej)
!$omp do
      do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)

              ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
              bmat(ij)=(fact)*ak*ek

              ek=sdot(nao,P_ia(io,jv,:),1,Q_ia(jo,iv,:),1)
              bmat(ij)=bmat(ij)-fact*ax*ek ! scaled by ax
           enddo
           ij=lin(i,i)
           ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
           bmat(ij)=fact*(ak*ek-ax*ek) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
write(52)bmat
close(52)
deallocate(bmat,P_ia,Q_ia)
return

end subroutine XSrtdacorr

subroutine uXstda_mat(nci,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,mxcnfb,iconfa,iconfb,dax,eda,edb,hci,alpha,beta,&
&xyz,nao,mocia,mocib,clowa,clowb,epsia,epsib) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,nexa,nexb,ncent,noa,nva,mxcnfa,iconfa(mxcnfa,2),nob,nvb,mxcnfb,iconfb(mxcnfb,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dax,eda(mxcnfa),edb(mxcnfb),alpha,beta
real*4, intent(out) :: hci(nci,nci)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ax,de,alphak,betaj

real*4 :: start_time,end_time,start
integer :: inoa,nnoa,inva,nnva
integer :: inob,nnob,invb,nnvb

integer :: nao,mocia,mocib
real*8  :: clowa(nao*mocia),clowb(nao*mocib)

real*4,allocatable :: Q_iaA(:,:,:),P_iaA(:,:,:),Q_ijA(:,:),Q_abA(:,:),P_ijA(:,:),AABB_integral(:,:)
real*4,allocatable :: Q_iaB(:,:,:),P_iaB(:,:,:),Q_ijB(:,:),Q_abB(:,:),P_ijB(:,:)
real*4 :: value1,value2,value3,rabx
real*8 :: epsia(mocia),epsib(mocib)

alphaK=real(alpha,4)
betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
inoa=minval(iconfa(1:nexa,1))
nnoa=maxval(iconfa(1:nexa,1))
inva=minval(iconfa(1:nexa,2)) ! start with no+1
nnva=maxval(iconfa(1:nexa,2))
write(*,*)'alpha MOs'
write(*,*)'Q transition charges computed for occ.',inoa,'to',nnoa
write(*,*)'and for unocc.',inva,'to',nnva

allocate(Q_iaA(inoa:nnoa,inva:nnva,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inoa, nnoa
Do j=inva,nnva
Do k=1,nao
    Q_iaA(i,j,k)=clowa(k+(i-1)*nao)*clowa(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_iaA(inoa:nnoa,inva:nnva,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inoa, nnoa
    Do j=inva,nnva

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_iaA(i,j,:),1,0.0,P_iaA(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao),P_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
Do k=inoa, i
Do j=1, nao
    Q_ijA(lin(i,k),j)=clowa(j+(i-1)*nao)*clowa(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
    Do k=inoa, nnoa

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ijA(lin(i,k),:),1,0.0,P_ijA(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ijA)

allocate(Q_abA(inva*(inva+1)/2:nnva*(nnva+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inva,nnva
Do l=inva, j
Do i=1,nao
    Q_abA(lin(j,l),i)=clowa(i+(j-1)*nao)*clowa(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0
call cpu_time(start_time)

! reduce the mo range to match the configuration space
inob=minval(iconfb(1:nci-nexa,1))
nnob=maxval(iconfb(1:nci-nexa,1))
invb=minval(iconfb(1:nci-nexa,2)) ! start with no+1
nnvb=maxval(iconfb(1:nci-nexa,2))
write(*,*)'beta MOs'
write(*,*)'Q transition charges computed for occ.',inob,'to',nnob
write(*,*)'and for unocc.',invb,'to',nnvb

allocate(Q_iaB(inob:nnob,invb:nnvb,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inob, nnob
Do j=invb,nnvb
Do k=1,nao
    Q_iaB(i,j,k)=clowb(k+(i-1)*nao)*clowb(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_iaB(inob:nnob,invb:nnvb,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inob, nnob
    Do j=invb,nnvb

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_iaB(i,j,:),1,0.0,P_iaB(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ijB(inob*(inob+1)/2:nnob*(nnob+1)/2,nao),P_ijB(inob*(inob+1)/2:nnob*(nnob+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=inob, nnob
Do k=inob, i
Do j=1, nao
    Q_ijB(lin(i,k),j)=clowb(j+(i-1)*nao)*clowb(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=inob, nnob
    Do k=inob, nnob

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ijB(lin(i,k),:),1,0.0,P_ijB(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ijB,AABB_integral)

allocate(Q_abB(invb*(invb+1)/2:nnvb*(nnvb+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=invb,nnvb
Do l=invb, j
Do i=1,nao
    Q_abB(lin(j,l),i)=clowb(i+(j-1)*nao)*clowb(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ax=real(dax)
! calculate A+B and A-B
hci=0.0e0

!alpha-alpha block
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
  do i=1,nexa
     io=iconfa(i,1)
     iv=iconfa(i,2)
     do j=1,i-1
        jo=iconfa(j,1)
        jv=iconfa(j,2)
        ! ek = (ia|jb)
ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(jo,jv,:),1)
        !  ej = (ij|ab)
ej=ax*sdot(nao,P_ijA(lin(io,jo),:),1,Q_abA(lin(iv,jv),:),1)
        hci(i,j)=ek-ej
        hci(j,i)=hci(i,j)
     enddo ! j
ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(io,iv,:),1)
ej=ax*sdot(nao,P_ijA(lin(io,io),:),1,Q_abA(lin(iv,iv),:),1)
de=real(epsia(iconfa(i,2))-epsia(iconfa(i,1)),4)
           hci(i,i)=de+ek-ej
  enddo ! i
!$omp end do
!$omp end parallel

! beta blocks
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
      do i = nexa+1,nci
         io=iconfb(i-nexa,1)
         iv=iconfb(i-nexa,2)
! alpha-beta
         do j = 1,nexa
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaA(jo,jv,:),1)
            hci(j,i)=ek
            hci(i,j)=ek
         enddo
! beta-beta
         do j = nexa+1,i-1
            jo=iconfb(j-nexa,1)
            jv=iconfb(j-nexa,2)
            ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(jo,jv,:),1)
            ej=ax*sdot(nao,P_ijB(lin(io,jo),:),1,Q_abB(lin(iv,jv),:),1)
            hci(j,i)=ek-ej
            hci(i,j)=hci(j,i)
         enddo
         ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(io,iv,:),1)
         ej=ax*sdot(nao,P_ijB(lin(io,io),:),1,Q_abB(lin(iv,iv),:),1)
         de=real(epsib(iconfb(i-nexa,2))-epsib(iconfb(i-nexa,1)),4)
         hci(i,i)=de+ek-ej
      enddo
!$omp end do
!$omp end parallel

deallocate(P_iaA,Q_iaA,P_ijA,Q_abA,P_iaB,Q_iaB,P_ijB,Q_abB)

return

end subroutine uXstda_mat

subroutine SF_Xstda_mat(nci,nexb,ncent,noa,nva,nob,nvb,mxcnfb,iconfb,dax,edb,hci,beta,&
&xyz,nao,mocia,mocib,clowa,clowb,epsia,epsib) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,nexb,ncent,noa,nva,nob,nvb,mxcnfb,iconfb(mxcnfb,2)
real*8, intent(in)  :: xyz(4,ncent)
real*8, intent(in)  :: dax,edb(mxcnfb),beta
real*4, intent(out) :: hci(nci,nci)
integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ax,de,betaj

real*4 :: start_time,end_time,start
integer :: inoa,nnoa,inva,nnva
integer :: inob,nnob,invb,nnvb

integer :: nao,mocia,mocib
real*8  :: clowa(nao*mocia),clowb(nao*mocib)

real*4,allocatable :: Q_ijA(:,:),P_ijA(:,:),AABB_integral(:,:)
real*4,allocatable :: Q_abB(:,:)
real*4 :: value1,value2,value3,rabx
real*8 :: epsia(mocia),epsib(mocib)

betaJ=real(beta,4)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
inoa=minval(iconfb(1:nci,1))
nnoa=maxval(iconfb(1:nci,1))
write(*,*)'alpha MOs'
write(*,*)'Q transition charges computed for occ.',inoa,'to',nnoa

allocate(Q_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao),P_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
Do k=inoa, i
Do j=1, nao
    Q_ijA(lin(i,k),j)=clowa(j+(i-1)*nao)*clowa(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
    Do k=inoa, nnoa

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ijA(lin(i,k),:),1,0.0,P_ijA(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ijA)


write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0
call cpu_time(start_time)

! reduce the mo range to match the configuration space
invb=minval(iconfb(1:nci,2)) ! start with no+1
nnvb=maxval(iconfb(1:nci,2))
write(*,*)'beta MOs'
write(*,*)'and for unocc.',invb,'to',nnvb

allocate(Q_abB(invb*(invb+1)/2:nnvb*(nnvb+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=invb,nnvb
Do l=invb, j
Do i=1,nao
    Q_abB(lin(j,l),i)=clowb(i+(j-1)*nao)*clowb(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ax=real(dax)
! calculate A+B and A-B
hci=0.0e0

! beta blocks
!$omp parallel private(i,j,io,iv,jo,jv,ek,ej,de)
!$omp do
      do i = 1,nci
         io=iconfb(i,1)
         iv=iconfb(i,2)
! beta-beta
         do j = 1,i-1
            jo=iconfb(j,1)
            jv=iconfb(j,2)
            ej=ax*sdot(nao,P_ijA(lin(io,jo),:),1,Q_abB(lin(iv,jv),:),1)
            hci(j,i)=-ej
            hci(i,j)=hci(j,i)
         enddo
         ej=ax*sdot(nao,P_ijA(lin(io,io),:),1,Q_abB(lin(iv,iv),:),1)
         de=real(epsib(iconfb(i,2))-epsia(iconfb(i,1)),4)
         hci(i,i)=de-ej
      enddo
!$omp end do
!$omp end parallel

deallocate(P_ijA,Q_abB)

return

end subroutine SF_Xstda_mat

subroutine uXsrpamat(nex,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,&
                     &mxcnfb,iconfa,iconfb,dax,eda,edb,ambsqr,&
                     &alpha,beta,xyz,nao,mocia,mocib,clowa,clowb,epsia,epsib)
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nex,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,nao,mocia,mocib
real*8, intent(in)  :: xyz(4,ncent)
integer, intent(in) :: mxcnfb,iconfa(mxcnfa,2),iconfb(mxcnfb,2)
real*4, intent(out) :: ambsqr(nex*(nex+1)/2)
real*8, intent(in)  :: dax,eda(mxcnfa),edb(mxcnfb),alpha,beta
real*4, allocatable :: apb(:)
integer i,j,ij,io,iv,jo,jv,ierr,lin,k,iiv,jjv,iwrk,jwrk,iwrk2,l,a,b,m,n,o,p,aa,bb,ia,jb,nci
real*4 ek,ej,sdot,ax,de,alphak,betaj
real*8  :: clowa(nao*mocia),clowb(nao*mocib)
real*8 :: epsia(mocia),epsib(mocib)
real*4 :: start_time,end_time,start
integer :: inoa,nnoa,inva,nnva
integer :: inob,nnob,invb,nnvb
real*4,allocatable :: Q_iaA(:,:,:),P_iaA(:,:,:),Q_ijA(:,:),Q_abA(:,:),P_ijA(:,:),AABB_integral(:,:)
real*4,allocatable :: Q_iaB(:,:,:),P_iaB(:,:,:),Q_ijB(:,:),Q_abB(:,:),P_ijB(:,:)
real*4 :: value1,value2,value3,rab
nci=nex
alphaK=real(alpha,4)
betaJ=real(beta,4)
write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
inoa=minval(iconfa(1:nexa,1))
nnoa=maxval(iconfa(1:nexa,1))
inva=minval(iconfa(1:nexa,2)) ! start with no+1
nnva=maxval(iconfa(1:nexa,2))
write(*,*)'alpha MOs'
write(*,*)'Q transition charges computed for occ.',inoa,'to',nnoa
write(*,*)'and for unocc.',inva,'to',nnva

allocate(Q_iaA(inoa:nnoa,inva:nnva,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inoa, nnoa
Do j=inva,nnva
Do k=1,nao
    Q_iaA(i,j,k)=clowa(k+(i-1)*nao)*clowa(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_iaA(inoa:nnoa,inva:nnva,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inoa, nnoa
    Do j=inva,nnva

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_iaA(i,j,:),1,0.0,P_iaA(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao),P_ijA(inoa*(inoa+1)/2:nnoa*(nnoa+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
Do k=inoa, i
Do j=1, nao
    Q_ijA(lin(i,k),j)=clowa(j+(i-1)*nao)*clowa(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=inoa, nnoa
    Do k=inoa, nnoa

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ijA(lin(i,k),:),1,0.0,P_ijA(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ijA)

allocate(Q_abA(inva*(inva+1)/2:nnva*(nnva+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=inva,nnva
Do l=inva, j
Do i=1,nao
    Q_abA(lin(j,l),i)=clowa(i+(j-1)*nao)*clowa(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0
call cpu_time(start_time)

! reduce the mo range to match the configuration space
inob=minval(iconfb(1:nci-nexa,1))
nnob=maxval(iconfb(1:nci-nexa,1))
invb=minval(iconfb(1:nci-nexa,2)) ! start with no+1
nnvb=maxval(iconfb(1:nci-nexa,2))
write(*,*)'beta MOs'
write(*,*)'Q transition charges computed for occ.',inob,'to',nnob
write(*,*)'and for unocc.',invb,'to',nnvb

allocate(Q_iaB(inob:nnob,invb:nnvb,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inob, nnob
Do j=invb,nnvb
Do k=1,nao
    Q_iaB(i,j,k)=clowb(k+(i-1)*nao)*clowb(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_iaB(inob:nnob,invb:nnvb,nao))

!$omp parallel private(i,j)
!$omp do
Do i=inob, nnob
    Do j=invb,nnvb

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_iaB(i,j,:),1,0.0,P_iaB(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

allocate(Q_ijB(inob*(inob+1)/2:nnob*(nnob+1)/2,nao),P_ijB(inob*(inob+1)/2:nnob*(nnob+1)/2,nao))

!$omp parallel private(i,k)
!$omp do
Do i=inob, nnob
Do k=inob, i
Do j=1, nao
    Q_ijB(lin(i,k),j)=clowb(j+(i-1)*nao)*clowb(j+(k-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,k)
!$omp do
Do i=inob, nnob
    Do k=inob, nnob

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ijB(lin(i,k),:),1,0.0,P_ijB(lin(i,k),:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(Q_ijB,AABB_integral)

allocate(Q_abB(invb*(invb+1)/2:nnvb*(nnvb+1)/2,nao))

!$omp parallel private(j,l)
!$omp do
Do j=invb,nnvb
Do l=invb, j
Do i=1,nao
    Q_abB(lin(j,l),i)=clowb(i+(j-1)*nao)*clowb(i+(l-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ax=real(dax)

open(unit=52,file='apbmat',form='unformatted',status='replace')
ambsqr=0.0e0

if(abs(dax).lt.1.0d-6) then
! calculate A+B and (A-B)^0.5
  allocate(apb(nex*(nex+1)/2),stat=ierr)
  if(ierr.ne.0)stop 'allocation for A+B crashed'
  apb=0.0e0
  ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,ek,de)
!$omp do
! alpha-alpha block
  do i = 1,nexa
     io=iconfa(i,1)
     iv=iconfa(i,2)
     do j=1,i-1
        ij=lin(i,j)
        jo=iconfa(j,1)
        jv=iconfa(j,2)
        ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(jo,jv,:),1)! ek = (ia|jb)
        apb(ij)=2.0*ek
        ambsqr(ij)=0.0
     enddo
     ij=lin(i,i)
     ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(io,iv,:),1)
     de=real(epsia(iconfa(i,2))-epsia(iconfa(i,1)),4)
     ambsqr(ij)=sqrt(de) ! diagonal element of (A-B)^0.5
     apb(ij)=de+ek*2.0 ! diagonal element of A+B
  enddo
!$omp end do
!$omp end parallel
!$omp parallel private(ij,i,j,io,iv,jo,jv,de,ek)
!$omp do
! beta block
  do i = nexa+1,nex
     io=iconfb(i-nexa,1)
     iv=iconfb(i-nexa,2)
!... with alpha
     do j = 1,nexa
        ij=lin(i,j)
        jo=iconfa(j,1)
        jv=iconfa(j,2)
        ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaA(jo,jv,:),1)! ek = (ia|jb)
        apb(ij)=2.0*ek
     enddo
!... with beta
     do j=nexa+1,i-1
        ij=lin(i,j)
        jo=iconfb(j-nexa,1)
        jv=iconfb(j-nexa,2)
        ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(jo,jv,:),1) ! ek = (ia|jb)
        apb(ij)=2.0*ek
     enddo
     ij=lin(i,i)
     ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(io,iv,:),1)
     de=real(epsib(iconfb(i-nexa,2))-epsib(iconfb(i-nexa,1)),4)
     ambsqr(ij)=sqrt(de) ! diagonal element of (A-B)^0.5
     apb(ij)=de+ek*2.0 ! diagonal element of A+B
  enddo
!$omp end do
!$omp end parallel
  write(52)apb
else

!********************
! now the hybrid case
! this is tedious, since we have to diagonalize A-B
! for alpha and beta seperately (to save time)
! since we diagonalize alpha and beta blocks of A-B separately, we temporarily use apb for this and allocate for each space separately
!********************
  allocate(apb(nexa*(nexa+1)/2),stat=ierr)
  apb=0.0e0
  if(ierr.ne.0)stop 'allocation failed for A-B vector'
! calculate A+B and A-B
! alpha-alpha block
  ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,de,ek,ej)
!$omp do
  do i = 1,nexa
     io=iconfa(i,1)
     iv=iconfa(i,2)
     do j=1,i-1
        ij=lin(i,j)
        jo=iconfa(j,1)
        jv=iconfa(j,2)
        ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(jo,jv,:),1) ! ek = (ia|jb)
        ej=ax*sdot(nao,P_ijA(lin(io,jo),:),1,Q_abA(lin(iv,jv),:),1) !  ej = (ij|ab)
        ambsqr(ij)=2.0*ek-ej
        ek=sdot(nao,P_iaA(io,jv,:),1,Q_iaA(jo,iv,:),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
        ! A+B part
        ambsqr(ij)=ambsqr(ij)-ax*ek
        ! we first use apb as A-B
        apb(ij)=ax*ek-ej
     enddo
     ij=lin(i,j)
     ek=sdot(nao,P_iaA(io,iv,:),1,Q_iaA(io,iv,:),1)
     ej=ax*sdot(nao,P_ijA(lin(io,io),:),1,Q_abA(lin(iv,iv),:),1)
     de=real(epsia(iconfa(i,2))-epsia(iconfa(i,1)),4)
     apb(ij)=de+ax*ek-ej! diagonal element of A-B
     ambsqr(ij)=de-ax*ek+2.0*ek-ej ! diagonal element of A+B
  enddo
!$omp end do
!$omp end parallel

! beta...
!$omp parallel private(ij,i,j,io,iv,jo,jv,de,ek,ej)
!$omp do
  do i = nexa+1,nex
     io=iconfb(i-nexa,1)
     iv=iconfb(i-nexa,2)
! ...alpha block
     do j = 1,nexa
        ij=lin(i,j)
        jo=iconfa(j,1)
        jv=iconfa(j,2)
        ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaA(jo,jv,:),1)! ek = (ia|jb)
        ambsqr(ij)=2.0*ek
     enddo
! ...beta block
     do j = nexa+1,i-1
        ij=ij+1
        jo=iconfb(j-nexa,1)
        jv=iconfb(j-nexa,2)
        ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(jo,jv,:),1) ! ek = (ia|jb)
        ej=ax*sdot(nao,P_ijB(lin(io,jo),:),1,Q_abB(lin(iv,jv),:),1) !  ej = (ij|ab)
        ambsqr(ij)=2.0*ek-ej
        ek=sdot(nao,P_iaB(io,jv,:),1,Q_iaB(jo,iv,:),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
        ambsqr(ij)=ambsqr(ij)-ax*ek ! scaled by ax
     enddo
     ij=lin(i,j)
     ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(io,iv,:),1)
     ej=ax*sdot(nao,P_ijB(lin(io,io),:),1,Q_abB(lin(iv,iv),:),1)
     de=real(epsib(iconfb(i-nexa,2))-epsib(iconfb(i-nexa,1)),4)
     ambsqr(ij)=de+2.0*ek-ax*ek-ej ! diagonal element of A+0.5*B  (beta part)
  enddo
!$omp end do
!$omp end parallel

  write(52)ambsqr

! this the time determining step, since A-B needs to be diagonalized in full nci space
  write(*,*) ' calculating (A-B)^0.5 ...'
  write(*,'('' estimated time (min) '',f8.2)') (float(nexb)**2*float(nexb)+float(nexa)**2*float(nexa))/4.d+8/60.

! take power of alpha-alpha block of A-B
  call smatpow(nexa,apb)
! blow up alpha-alpha block of (A-B)^0.5
  ij=0
  do i=1,nexa
    do j=1,i
      ij=ij+1
      ambsqr(ij)=apb(ij)
    enddo
  enddo

! free memory anf reallocate for beta-beta part
  deallocate(apb)
  allocate(apb(nexb*(nexb+1)/2),stat=ierr)
  if(ierr.ne.0)stop 'allocation failed for A-B vector'

! beta-beta block
  ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,de,ek,ej)
!$omp do
  do i = 1,nexb
     io=iconfb(i,1)
     iv=iconfb(i,2)
     do j = 1,i-1
        ij=lin(i,j)
        jo=iconfb(j,1)
        jv=iconfb(j,2)
        ej=ax*sdot(nao,P_ijB(lin(io,jo),:),1,Q_abB(lin(iv,jv),:),1) !  ej = (ij|ab)
        ek=sdot(nao,P_iaB(io,jv,:),1,Q_iaB(jo,iv,:),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
        ! we first use apb as A-B
        apb(ij)=ax*ek-ej
     enddo
     ij=lin(i,i)
     ek=sdot(nao,P_iaB(io,iv,:),1,Q_iaB(io,iv,:),1)
     ej=ax*sdot(nao,P_ijB(lin(io,io),:),1,Q_abB(lin(iv,iv),:),1)
     de=real(epsib(iconfb(i,2))-epsib(iconfb(i,1)),4)
     apb(ij)=de+ax*ek-ej ! diagonal element of A-B
  enddo
!$omp end do
!$omp end parallel

! take power of beta-beta block of A-B
  call smatpow(nexb,apb)
! blow up beta-beta block of (A-B)^0.5

  k=nexa*(nexa+1)
  k=k/2
  ij=0
  do i=nexa+1,nex
    do j=1,nexa
      k=k+1
      ambsqr(k)=0.0e0
    enddo
    do j=nexa+1,i
      k=k+1
      ij=ij+1
      ambsqr(k)=apb(ij)
    enddo
  enddo
! free memory and reallocate for beta-beta part
  deallocate(apb)
endif ! GGA/hybrid case
close(52)

deallocate(P_iaA,Q_iaA,P_ijA,Q_abA,P_iaB,Q_iaB,P_ijB,Q_abB)

return

end subroutine uXsrpamat


subroutine RSH_Xsrtdacorr(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,clow,nao,moci)
use omp_lib
use commonlogicals
use stdacommon
use commonlibcint
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: dak,dax,ed(mxcnf)
real*4, allocatable :: bmat(:)
integer i,j,k,io,iv,jo,jv,ierr,iiv,jjv,iwrk,jwrk,lin,ij
real*4 ek,ej,sdot,ak,ax,de,fact,ek_RSH


real*4 :: start_time,end_time,start
integer :: nao,moci
integer :: ino,nno,inv,nnv
real*8  :: clow(nao*moci)
real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),AABB_integral(:,:),P_ia_RSH(:,:,:)

write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals read'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
write(*,*)'Q transition charges computed for occ.',ino,'to',nno
write(*,*)'and for unocc.',inv,'to',nnv

allocate(Q_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
Do j=inv,nnv
Do k=1,nao
    Q_ia(i,j,k)=clow(k+(i-1)*nao)*clow(k+(j-1)*nao)
enddo
enddo
enddo
!$omp end do
!$omp end parallel

allocate(P_ia(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

!RSH integrals long range

allocate(AABB_integral(nao,nao))

call two_elec_int_RSH(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)
allocate(P_ia_RSH(ino:nno,inv:nnv,nao))

!$omp parallel private(i,j)
!$omp do
Do i=ino, nno
    Do j=inv,nnv

    call sgemv('T',nao,nao,1.0,AABB_integral,nao,Q_ia(i,j,:),1,0.0,P_ia_RSH(i,j,:),1)

    enddo
enddo
!$omp end do
!$omp end parallel

deallocate(AABB_integral)

write(*,*)"Intermediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0

ij=nci
ij=ij*(ij+1)/2
allocate(bmat(ij), stat=ierr)
if(ierr.ne.0)stop 'allocation for bmat crashed'
ak=real(dak)
ax=real(dax)
! calculate 0.5*B
bmat=0.0e0
fact=0.50d0 ! this is the scaling of the B-contribution
open(unit=52,file='bmat',form='unformatted',status='replace')
ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,ek,ej)
!$omp do
      do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)

              ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(jo,jv,:),1)
              bmat(ij)=(fact)*ak*ek

              ek=ax*sdot(nao,P_ia(io,jv,:),1,Q_ia(jo,iv,:),1)&
              &+real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,jv,:),1,Q_ia(jo,iv,:),1)
              bmat(ij)=bmat(ij)-fact*ek ! scaled by ax
           enddo
           ij=lin(i,i)
           ek=sdot(nao,P_ia(io,iv,:),1,Q_ia(io,iv,:),1)
           ek_RSH=real(RSH_beta,4)*sdot(nao,P_ia_RSH(io,iv,:),1,Q_ia(io,iv,:),1)
           bmat(ij)=fact*(ak*ek-ax*ek-ek_RSH) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
write(52)bmat
close(52)
deallocate(bmat,P_ia,Q_ia)
return

end subroutine RSH_XSrtdacorr
