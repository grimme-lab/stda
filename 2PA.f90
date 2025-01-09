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
SUBROUTINE lresp_2PA_full(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,no,nv,eci,Xci,Yci,nroot,&
&ncent,dax,nao,clow)
use commonresp
use omp_lib
IMPLICIT NONE

integer ::i,j,k,ii,jj,kk,ij,jk,ab,io,iv,idum1,idum2,nci
integer ::io1,io2,iv1,iv2,iwrk,jwrk,nroot,ncent,nao
integer ::maxconf,moci,no,nv,ino,nno,inv,nnv
integer ::iconf(maxconf,2)
integer, allocatable :: A_list(:,:)
integer ::counter_A
integer, allocatable :: B_list(:,:)
integer ::counter_B
integer*8 ::lin8
real*8 :: dax,clow(nao*moci)

real*8 ::xl(moci*(moci+1)/2)
real*8 ::yl(moci*(moci+1)/2)
real*8 ::zl(moci*(moci+1)/2)

real*4 ::mu_x(nci)
real*4 ::mu_y(nci)
real*4 ::mu_z(nci)
real*4 ::XpY_int(nci,3)
real*4 ::XpYci(nci)

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


real*4 ::start_time,end_time,sdot

real*4,allocatable :: f_ijka(:,:),f_abic(:,:),F_ij(:,:),F_ab(:,:)

open(unit=60,file='2PA-abs',status='replace')

mu=0.0
mu(:,1)=xl(:)
mu(:,2)=yl(:)
mu(:,3)=zl(:)

write(*,*)
write(*,*)'======================================================================'
write(*,*)'               Welcome in nonlinear response sTD-DFT program'
write(*,*)'======================================================================'
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
allocate( X(nci,3),Y(nci,3))

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
           mu_x(j)=-xl(ij)
           mu_y(j)=-yl(ij)
           mu_z(j)=-zl(ij)
             enddo
!$omp end do
!$omp end parallel

! Generating Hartree XC kernel

ino=minval(iconf(1:nci,1))
nno=maxval(iconf(1:nci,1))
inv=minval(iconf(1:nci,2)) ! start with no+1
nnv=maxval(iconf(1:nci,2))
allocate(f_ijka(ino*(ino+1)/2:nno*(nno+1)/2,nci),f_abic(inv*(inv+1)/2:nnv*(nnv+1)/2,nci))
call HXC(nci,ncent,no,nv,maxconf,iconf,dax,nao,moci,clow,f_ijka,f_abic,ino,nno,inv,nnv)

allocate(F_ij(ino*(ino+1)/2:nno*(nno+1)/2,4),F_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,4))

Do ii=1, num_trans

XpY(:,:)=0.0
XmY(:,:)=0.0
X(:,:)=0.0
Y(:,:)=0.0
omega=-eci(ii)/2.0
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

 write(*,*)
 write(*,*)ii
 write(*,*)

! extract X and Y from XpY
! (X-Y)=omega*(A-B)^-1 (X+Y)
! X=((X+Y)+(X-Y))/2
! Y=(X+Y)-X
!$omp parallel private(i,j,ij,ix) reduction (+:XmY)
!$omp do
      Do i=1,nci
        Do j=1,nci
        ij=lin8(i,j)
        Do ix=1,3
        XmY(i,ix)=XmY(i,ix)+ dble(omega)*dble(inv_amb(ij))*XpY(j,ix)
        enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(ix,j) reduction(+:X,Y)
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
print '("alpha Time = ",f12.2," minutes.")',(end_time-start_time)/60.0

XpYci(:)=Xci(:,ii)+Yci(:,ii)
!$omp parallel private(i,j,ix)
!$omp do
Do i=ino,nno
Do j=ino,i
Do ix=1,3
F_ij(lin8(i,j),ix)=sdot(nci,f_ijka(lin8(i,j),:),1,XpY_int(:,ix),1)
enddo
F_ij(lin8(i,j),4)=sdot(nci,f_ijka(lin8(i,j),:),1,XpYci(:),1)
enddo
enddo
!$omp end do
!$omp end parallel
!$omp parallel private(i,j,ix)
!$omp do
Do i=inv,nnv
Do j=inv,i
Do ix=1,3
F_ab(lin8(i,j),ix)=sdot(nci,f_abic(lin8(i,j),:),1,XpY_int(:,ix),1)
enddo
F_ab(lin8(i,j),4)=sdot(nci,f_abic(lin8(i,j),:),1,XpYci(:),1)
enddo
enddo
!$omp end do
!$omp end parallel
write(*,*)'F_ij and F_ab computed'

if(ii==1)then
!
! Genarating a list of indexes used in A and B formula to save a great bunch of time
!
call cpu_time(start_time)
counter_A=0
!$omp parallel private(j,i,kk) reduction(+:counter_A)
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
!$omp parallel private(j,i,kk) reduction(+:counter_B)
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
print '("A & B indexes list Time = ",f12.2," minutes.")',(end_time-start_time)/60.0
endif

! sigma  sigma = -A + B
call cpu_time(start_time)
!
! Fast version
!
sigma(:,:)=0.0
Do ix=1,3
Do iy=1,ix
call TPA_resp_fast_full(ix,iy,X,Y,Xci,Yci,nroot,A_list,B_list,counter_A,counter_B,mu,&
&maxconf,no,nv,nci,moci,ii,iconf,A,B,ino,nno,inv,nnv,F_ij,F_ab)
sigma(ix,iy)=(-A+B)/2.0d0
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
print '("2PA  Time = ",f12.2," minutes.")',(end_time-start_time)/60.0
write(*,*)
write(*,3333)'Delta (',eci(ii)*27.21139,')'
write(*,*)
write(*,1111)'x','y','z'
write(*,2222)'x',sigma(1,1),sigma(1,2),sigma(1,3)
write(*,2222)'y',sigma(2,1),sigma(2,2),sigma(2,3)
write(*,2222)'z',sigma(3,1),sigma(3,2),sigma(3,3)
write(*,*)
write(*,5555)'F =',sigma_f,' G =',sigma_G,' H =',sigma_H
write(*,4444)'Delta_2PA_//   =',2.0*sigma_f+2.0*sigma_g+2.0*sigma_h
write(*,4444)'Delta_2PA__|_  =',-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h
write(*,4444)'Delta_2PA_circ =',-2.0*sigma_f+3.0*sigma_g+3.0*sigma_h
write(*,4444)'rho = //*(_|_)**-1 =',(2.0*sigma_f+2.0*sigma_g+2.0*sigma_h)/&
&(-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h)
write(60,6666)ii,eci(ii)*27.21139,2.0*sigma_f+2.0*sigma_g+2.0*sigma_h,&
&-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h,-2.0*sigma_f&
&+3.0*sigma_g+3.0*sigma_h,(2.0*sigma_f+2.0*sigma_g&
&+2.0*sigma_h)/(-1.0*sigma_f+4.0*sigma_g-1.0*sigma_h)
enddo
close(60)
deallocate(XpY)
deallocate(XmY)
deallocate(X,Y)
deallocate(inv_amb)
write(*,*)
write(*,*)'======================================================================'
write(*,*)'               end of     nonlinear response sTD-DFT program'
write(*,*)'======================================================================'
write(*,*)
111   format(A15,F20.6)
1111  format(A22,2A20)
2222  format(A2,3F20.6)
3334  format(A16,F7.3,A1,F7.3,A1)
3333  format(A16,F7.3,A1)
4444  format(A20,F20.3)
5555  format(A3,F20.3,A4,F20.3,A4,F20.3)
6666  format(I3,F7.3,4F20.3)
end subroutine lresp_2PA_full





subroutine HXC(nci,ncent,no,nv,mxcnf,iconf,dax,nao,moci,clow,f_ijka,f_abic,ino,nno,inv,nnv) ! (AA|BB) only
use commonlogicals
use stdacommon
use commonlibcint
use omp_lib
implicit none
integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
real*8, intent(in)  :: dax
real*4, intent(out) :: f_ijka(ino*(ino+1)/2:nno*(nno+1)/2,nci),f_abic(inv*(inv+1)/2:nnv*(nnv+1)/2,nci)
integer i,j,ij,io,iv,ko,kv,ierr,lin,iiv,jjv,iwrk,jwrk,iwrk2,k,l,a,b,m,n,o,p,aa,bb,ia,jb
real*4 ek,ej,sdot,ax

real*8 :: wtime
real*4 :: start_time,end_time,start
integer :: ino,nno,inv,nnv

integer :: nao,moci
real*8  :: clow(nao*moci)

real*4,allocatable :: Q_ia(:,:,:),P_ia(:,:,:),Q_ij(:,:),Q_ab(:,:),AABB_integral(:,:)

real*4 :: value1,value2,value3,rabx
Write(*,*)'Compute f_HXC once'
write(*,*)'using on_site AO integrals'

! Read AO integrals data

allocate(AABB_integral(nao,nao))
call two_elec_int(ncent,nao,nao,AABB_integral,1,nbas,1,nbas)


write(*,*)'AO integrals computed'

call cpu_time(start_time)

! reduce the mo range to match the configuration space
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

allocate(Q_ij(ino*(ino+1)/2:nno*(nno+1)/2,nao),Q_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,nao))

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


deallocate(AABB_integral)


write(*,*)"Intemediates computed"
call cpu_time(end_time)
print '("time = ",f12.2," minutes.")',(end_time-start_time)/60.0


ax=real(dax)
! calculate f_ijka and f_abic
f_ijka=0.0e0
f_abic=0.0e0
!$omp parallel private(i,j,k,ko,kv,ek,ej)
!$omp do
Do i=ino,nno
Do j=ino,i
Do k=1,nci
ko=iconf(k,1)
kv=iconf(k,2)
ej=sdot(nao,Q_ij(lin(i,j),:),1,P_ia(ko,kv,:),1)
ek=-ax*sdot(nao,Q_ij(lin(i,ko),:),1,P_ia(j,kv,:),1)
f_ijka(lin(i,j),k)=ej+ek
enddo
enddo
enddo
!$omp end do
!$omp end parallel

!$omp parallel private(i,j,k,ko,kv,ek,ej)
!$omp do
Do i=inv,nnv
Do j=inv,i
Do k=1,nci
ko=iconf(k,1)
kv=iconf(k,2)
ej=sdot(nao,Q_ab(lin(i,j),:),1,P_ia(ko,kv,:),1)
ek=-ax*sdot(nao,Q_ab(lin(j,kv),:),1,P_ia(ko,i,:),1)
f_abic(lin(i,j),k)=ej+ek
enddo
enddo
enddo
!$omp end do
!$omp end parallel

deallocate(P_ia,Q_ia,Q_ij,Q_ab)
write(*,*)'f_HXC computed'

return

end subroutine HXC

SUBROUTINE TPA_resp_fast_full(ix,iy,X,Y,Xci,Yci,nroot,A_list,B_list,counter_A,counter_B,mu,maxconf,no,nv,nci,moci,ii,&
&iconf,A,B,ino,nno,inv,nnv,F_ij,F_ab)
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
real*8 ::A1,A2,A3,A4,A5,A6
real*8 ::B1,B2,B3,B4,B5,B6
integer ::counter_A,counter_B
integer :: A_list(1:counter_A,1:3)
integer :: B_list(1:counter_B,1:3)

integer :: ino,nno,inv,nnv
real*4 :: F_ij(ino*(ino+1)/2:nno*(nno+1)/2,4),F_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,4)

A=0.0d0
B=0.0d0
A1=0.0d0
A2=0.0d0
A3=0.0d0
A4=0.0d0
A5=0.0d0
A6=0.0d0
B1=0.0d0
B2=0.0d0
B3=0.0d0
B4=0.0d0
B5=0.0d0
B6=0.0d0

! A1 ix iy n
      xx=ix
      yy=iy
call A_2PA_1_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A1,F_ij,ino,nno,iconf,maxconf)
! A2 iy ix n
      xx=iy
      yy=ix
call A_2PA_1_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A2,F_ij,ino,nno,iconf,maxconf)
! A3 n ix iy
      xx=ix
      yy=iy
call A_2PA_2_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A3,F_ij,ino,nno,iconf,maxconf)
! A4 n iy ix
      xx=iy
      yy=ix
call A_2PA_2_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A4,F_ij,ino,nno,iconf,maxconf)
! A5 ix n iy
      xx=ix
      yy=iy
call A_2PA_3_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A5,F_ij,ino,nno,iconf,maxconf)
! A6 iy n ix
      xx=iy
      yy=ix
call A_2PA_3_fast_full(xx,yy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ii,A6,F_ij,ino,nno,iconf,maxconf)
      A=A1+A2+A3+A4+A5+A6

! B1 ix iy n
      xx=ix
      yy=iy
call B_2PA_1_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B1,F_ab,inv,nnv,iconf,maxconf)
! B2 iy ix n
      xx=iy
      yy=ix
call B_2PA_1_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B2,F_ab,inv,nnv,iconf,maxconf)
! B3 n ix iy
      xx=ix
      yy=iy
call B_2PA_2_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B3,F_ab,inv,nnv,iconf,maxconf)
! B4 n iy ix
      xx=iy
      yy=ix
call B_2PA_2_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B4,F_ab,inv,nnv,iconf,maxconf)
! B5 ix n iy
      xx=ix
      yy=iy
call B_2PA_3_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B5,F_ab,inv,nnv,iconf,maxconf)
! B6 iy n ix
      xx=iy
      yy=ix
call B_2PA_3_fast_full(xx,yy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ii,B6,F_ab,inv,nnv,iconf,maxconf)
      B=B1+B2+B3+B4+B5+B6

end subroutine TPA_resp_fast_full

Subroutine A_2PA_1_fast_full(ix,iy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ifreq,A,F_ij,ino,nno,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,ino,nno,maxconf
integer ::nroot
real*8  ::A
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ij(ino*(ino+1)/2:nno*(nno+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_A
integer ::A_list(1:counter_A,1:3)

ii=ifreq
A=0.0
!$omp parallel private(i) reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+X(A_list(i,1),ix)*(-mu(A_list(i,2),iy)&
     &+dble(F_ij(lin8(iconf(A_list(i,1),1),iconf(A_list(i,3),1)),iy))&
     &)*dble(Yci(A_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
end subroutine A_2PA_1_fast_full

Subroutine A_2PA_2_fast_full(ix,iy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ifreq,A,F_ij,ino,nno,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,ino,nno,maxconf
integer ::nroot
real*8  ::A
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ij(ino*(ino+1)/2:nno*(nno+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_A
integer ::A_list(1:counter_A,1:3)

ii=ifreq
A=0.0
!$omp parallel private(i) reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+dble(Xci(A_list(i,1),ii))*(-mu(A_list(i,2),ix)&
      &+dble(F_ij(lin8(iconf(A_list(i,1),1),iconf(A_list(i,3),1)),ix))&
      &)*Y(A_list(i,3),iy)
      enddo
!$omp end do
!$omp end parallel
end subroutine A_2PA_2_fast_full

Subroutine A_2PA_3_fast_full(ix,iy,X,Y,Xci,Yci,nroot,A_list,counter_A,mu,nci,moci,ifreq,A,F_ij,ino,nno,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,ino,nno,maxconf
integer ::nroot
real*8  ::A
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ij(ino*(ino+1)/2:nno*(nno+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_A
integer ::A_list(1:counter_A,1:3)

ii=ifreq
A=0.0
!$omp parallel private(i) reduction(+:A)
!$omp do
      Do i=1,counter_A
      A=A+X(A_list(i,1),ix)*(&
     &dble(F_ij(lin8(iconf(A_list(i,1),1),iconf(A_list(i,3),1)),4))&
     &)*dble(Y(A_list(i,3),iy))
      enddo
!$omp end do
!$omp end parallel
end subroutine A_2PA_3_fast_full

Subroutine B_2PA_1_fast_full(ix,iy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ifreq,B,F_ab,inv,nnv,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,inv,nnv,maxconf
integer ::nroot
real*8  ::B
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_B
integer ::B_list(1:counter_B,1:3)

ii=ifreq
B=0.0
!$omp parallel private(i) reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+X(B_list(i,1),ix)*(-mu(B_list(i,2),iy)&
      &+dble(F_ab(lin8(iconf(B_list(i,1),2),iconf(B_list(i,3),2)),iy))&
      &)*dble(Yci(B_list(i,3),ii))
      enddo
!$omp end do
!$omp end parallel
end subroutine B_2PA_1_fast_full

Subroutine B_2PA_2_fast_full(ix,iy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ifreq,B,F_ab,inv,nnv,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,inv,nnv,maxconf
integer ::nroot
real*8  ::B
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_B
integer ::B_list(1:counter_B,1:3)

ii=ifreq
B=0.0
!$omp parallel private(i) reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+dble(Xci(B_list(i,1),ii))*(-mu(B_list(i,2),ix)&
      &+dble(F_ab(lin8(iconf(B_list(i,1),2),iconf(B_list(i,3),2)),ix))&
      &)*Y(B_list(i,3),iy)
      enddo
!$omp end do
!$omp end parallel
end subroutine B_2PA_2_fast_full

Subroutine B_2PA_3_fast_full(ix,iy,X,Y,Xci,Yci,nroot,B_list,counter_B,mu,nci,moci,ifreq,B,F_ab,inv,nnv,iconf,maxconf)
use commonresp
use omp_lib
implicit none

integer ::ix,iy,nci,ifreq,moci,inv,nnv,maxconf
integer ::nroot
real*8  ::B
real*8 ::mu(moci*(moci+1)/2,3)
real*8 ::X(nci,3)
real*8 ::Y(nci,3)
real*4 ::Xci(nci,nroot), Yci(nci,nroot)
real*4 :: F_ab(inv*(inv+1)/2:nnv*(nnv+1)/2,4)
integer ::i,ii
integer*8 :: lin8
integer ::iconf(maxconf,2)
integer ::counter_B
integer ::B_list(1:counter_B,1:3)

ii=ifreq
B=0.0
!$omp parallel private(i) reduction(+:B)
!$omp do
      Do i=1,counter_B
      B=B+X(B_list(i,1),ix)*(&
      &+dble(F_ab(lin8(iconf(B_list(i,1),2),iconf(B_list(i,3),2)),4))&
      &)*dble(Y(B_list(i,3),iy))
      enddo
!$omp end do
!$omp end parallel
end subroutine B_2PA_3_fast_full
