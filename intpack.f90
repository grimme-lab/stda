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
!
! Adapted from xtb4stda by MdeW
!! ------------------------------------------------------------------------
module intpack
   use iso_fortran_env, only : wp => real64

contains
! propa split in propa0 and propa1 to avoid computing several time the same things and neglect if efact below thr
subroutine propa0(aname,c,va,nt,ij,kl,efact)
   use stdacommon
   implicit real(wp)(a-h,o-z)
   external aname
   ! aufpunkte,ref point,intarray
   real(wp) c(3),va(nt)
   ! local
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib,gama, &
   & d(3),dd(84),e(3),aa(20),bb(20),a(3),b(3)
   dimension v(nt),val(nt)

   a(1:3)=eta(ij,1:3)
   b(1:3)=eta(kl,1:3)
   etaij4=eta(ij,4)
   etakl4=eta(kl,4)
   iff1=dint(eta(ij,5))
   iff2=dint(eta(kl,5))
   v=0.0
   ! --- a,b are centres of gaussians
   ra = a
   rb = b
   ! --- ga,gb are their exponents
   ga=etaij4
   gb=etakl4
   aa = 0
   bb = 0
   ia=iff1
   ib=iff2
   aa(ia)=1.0d0
   bb(ib)=1.0d0
   ! --- ia,ib are the canonical indices for monom
   ! --- apply product theorem
   cij=0.0d0
   ckl=0.0d0

   call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)
   !     if(mprp.eq.16) goto 200

   ! --- calculate cartesian prefactor for first gaussian
   call rhftce(aa,a,e,iff1)
   ! --- calculate cartesian prefactor for second gaussian
   call rhftce(bb,b,e,iff2)
   !     --- form their product
   call prod(aa,bb,dd,iff1,iff2)
   ! ----- e is center of product gaussian with exponent gama
   ! ----- c is reference point
   d = e - c
   end subroutine propa0

subroutine propa1(aname,c,va,nt,ij,kl,efact)
   use stdacommon
   implicit real(wp)(a-h,o-z)
   external aname
   ! aufpunkte,ref point,intarray
   real(wp) c(3),va(nt)
   ! local
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib,gama, &
   & d(3),dd(84),e(3),aa(20),bb(20),a(3),b(3)
   common/ prptyp / mprp
   dimension v(nt),val(nt)
   integer,parameter :: lin(84) = &
      & (/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0,2,2,0, &
      &   2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,3,0,5,5, &
      &   1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/)
   integer,parameter :: min(84) = &
      & (/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2, &
      &   1,2,1,0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0, &
      &   0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2/)
   integer,parameter :: nin(84) = &
      & (/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2, &
      &   1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1, &
      &    5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2/)

   iff1=ia
   iff2=ib
   val = 0

   !     d = e
   !   aname represents an external function
   if(mprp.eq.16) goto 200 !for magnetic moment
   if(iff1.gt.10.or.iff2.gt.10) goto 110
   if(iff1.gt.4.or.iff2.gt.4) goto 120
   if(iff1.gt.1.or.iff2.gt.1) goto 130
   ! s-s - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   l=lin(1)
   m=min(1)
   n=nin(1)
   call aname(l,m,n,gama,v,d)
   do j=1,nt
      val(j)=dd(1)*v(j)+val(j)
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   130 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 131
   ! s-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,4
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   131 continue
   ! p-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,10
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   120 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 121
   ! s-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,10
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   121 continue
   if(iff1.gt.4.and.iff2.gt.4) goto 122
   ! p-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,20
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   122 continue
   ! d-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,35
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   110 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 111
   ! s-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,20
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   111 continue
   if(iff1.gt.4.and.iff2.gt.4) goto 112
   ! p-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,35
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   112 continue
   if(iff1.gt.10.and.iff2.gt.10) goto 113
   ! d-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,56
      if(dabs(dd(i)).le.1.d-8) cycle
      call aname(lin(i),min(i),nin(i),gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   113 continue
   ! f-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,84
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return

   ! ang mom
   !cha
   200 do i=1,3
          va(i)=0.d0
          d(i)=e(i)-c(i)
       enddo
       if (efact.eq.0.d0) return
       ! lll, mmm, nnn are dummy
       call aname(lll,mmm,nnn,gama,v,d)
       efact=aa(ia)*bb(ib)*efact
       do i=1,nt
          va(i)=v(i)*efact
       enddo
       return
   !cha
end subroutine propa1



!=======================================================================
! cartesian gaussian functions (6d,10f...)
! iff :
! s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
! 1 2   3   4   5     6     7     8   9  10
! fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
!   11   12    13    14    15    16    17    18   19    20
! a, are aufpunkte
! nt : # of returns
! va : integral
!=======================================================================

!! --------------------------------------------------------------[SAW1710]-
!     changed do loops, replaced spaceships
subroutine propa(aname,c,va,nt,ij,kl)
   use stdacommon
   implicit real(wp)(a-h,o-z)
   external aname
   ! aufpunkte,ref point,intarray
   real(wp) a(3),b(3),c(3),va(nt)
   ! local
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
   common/ prptyp / mprp
   dimension d(3),dd(84),v(nt),val(nt) !,v(3),val(3)
   dimension e(3),aa(20),bb(20)
   integer,parameter :: lin(84) = &
      & (/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0,2,2,0, &
      &   2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,3,0,5,5, &
      &   1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/)
   integer,parameter :: min(84) = &
      & (/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2, &
      &   1,2,1,0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0, &
      &   0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2/)
   integer,parameter :: nin(84) = &
      & (/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2, &
      &   1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1, &
      &    5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2/)

   a(1:3)=eta(ij,1:3)
   b(1:3)=eta(kl,1:3)
   etaij4=eta(ij,4)
   etakl4=eta(kl,4)
   iff1=dint(eta(ij,5))
   iff2=dint(eta(kl,5))
   v=0.0
   ! --- a,b are centres of gaussians
   ra = a
   rb = b
   ! --- ga,gb are their exponents
   ga=etaij4
   gb=etakl4
   aa = 0
   bb = 0
   ia=iff1
   ib=iff2
   aa(ia)=1.0d0
   bb(ib)=1.0d0
   ! --- ia,ib are the canonical indices for monom
   ! --- apply product theorem
   cij=0.0d0
   ckl=0.0d0

   call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)
   !     if(mprp.eq.16) goto 200

   ! --- calculate cartesian prefactor for first gaussian
   call rhftce(aa,a,e,iff1)
   ! --- calculate cartesian prefactor for second gaussian
   call rhftce(bb,b,e,iff2)
   !     --- form their product
   call prod(aa,bb,dd,iff1,iff2)
   val = 0
   ! ----- e is center of product gaussian with exponent gama
   ! ----- c is reference point
   d = e - c
   !     d = e
   !   aname represents an external function
   if(mprp.eq.16) goto 200 !for magnetic moment
   if(iff1.gt.10.or.iff2.gt.10) goto 110
   if(iff1.gt.4.or.iff2.gt.4) goto 120
   if(iff1.gt.1.or.iff2.gt.1) goto 130
   ! s-s - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   l=lin(1)
   m=min(1)
   n=nin(1)
   call aname(l,m,n,gama,v,d)
   do j=1,nt
      val(j)=dd(1)*v(j)+val(j)
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   130 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 131
   ! s-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,4
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   131 continue
   ! p-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,10
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   120 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 121
   ! s-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,10
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   121 continue
   if(iff1.gt.4.and.iff2.gt.4) goto 122
   ! p-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,20
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   122 continue
   ! d-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,35
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   110 continue
   if(iff1.ne.1.and.iff2.ne.1) goto 111
   ! s-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,20
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   111 continue
   if(iff1.gt.4.and.iff2.gt.4) goto 112
   ! p-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,35
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   112 continue
   if(iff1.gt.10.and.iff2.gt.10) goto 113
   ! d-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,56
      if(dabs(dd(i)).le.1.d-8) cycle
      call aname(lin(i),min(i),nin(i),gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return
   113 continue
   ! f-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   do i=1,84
      if(dabs(dd(i)).le.1.d-8) cycle
   l=lin(i)
   m=min(i)
   n=nin(i)
      call aname(l,m,n,gama,v,d)
      do j=1,nt
         val(j)=dd(i)*v(j)+val(j)
      enddo
   enddo
   do i=1,nt
      va(i)=efact*val(i)
   enddo
   return

   ! ang mom
   !cha
   200 do i=1,3
          va(i)=0.d0
          d(i)=e(i)-c(i)
       enddo
       if (efact.eq.0.d0) return
       ! lll, mmm, nnn are dummy
       call aname(lll,mmm,nnn,gama,v,d)
       efact=aa(ia)*bb(ib)*efact
       do i=1,nt
          va(i)=v(i)*efact
       enddo
       return
   !cha
end subroutine propa


!! --------------------------------------------------------------[SAW1710]-
!     changed do loops
pure subroutine divpt(a,alpha,b,beta,ca,cb,e,fgama,ffact)
   implicit none
   !     this computes the center, exponent, and multiplying factor of
   !     a single gaussian which can replace the product of two gaussia
   !     centers a and b, and exponents alpha and beta.
   real(wp), intent(in)  :: a(3),alpha,b(3),beta,ca,cb
   real(wp), intent(out) :: fgama,ffact,e(3)
   integer :: i
   real(wp)  :: absqd,abprod,tol
   real(wp), parameter :: toluol = 20.7232660d0
   !      abexp=alpha+beta
   do i=1,3
      e(i)=(alpha*a(i)+beta*b(i))/(alpha+beta)
   enddo
   absqd=0.0d0
   do i=1,3
      absqd=absqd+(b(i)-a(i))*(b(i)-a(i))
   enddo
   abprod=absqd*alpha*beta/(alpha+beta)
   fgama=(alpha+beta)
   ffact=0.d0
   tol=toluol+ca+cb
   if(abprod.gt.tol) return
   ffact=dexp(-abprod)
   return
end subroutine divpt

!! --------------------------------------------------------------[SAW1710]-
!     made pure, made explicit
pure subroutine rhftce(cfs,a,e,iff)
   implicit none
   integer,intent(in)  :: iff
   real(wp), intent(in)  :: a(*),e(*)
   real(wp), intent(inout) :: cfs(*)
   real(wp), parameter   :: c2 = 2.0d0
   real(wp), parameter   :: c3 = 3.0d0
   real(wp)  :: aex,aey,aez
   ! ---- e = center of product function, a = center of single gaussian
   aex = e(1)-a(1)
   aey = e(2)-a(2)
   aez = e(3)-a(3)

   select case(iff)
   case(1)
      continue
   case(2)
      cfs(1)=aex*cfs(2)
   case(3)
      cfs(1)=aey*cfs(3)
   case(4)
      cfs(1)=aez*cfs(4)
   case(5)
      cfs(1)=aex*aex*cfs(5)
      cfs(2)=c2*aex*cfs(5)
   case(6)
      cfs(1)=aey*aey*cfs(6)
      cfs(3)=c2*aey*cfs(6)
   case(7)
      cfs(1)=aez*aez*cfs(7)
      cfs(4)=c2*aez*cfs(7)
   case(8)
      cfs(1)=aex*aey*cfs(8)
      cfs(2)=aey*cfs(8)
      cfs(3)=aex*cfs(8)
   case(9)
      cfs(1)=aex*aez*cfs(9)
      cfs(2)=aez*cfs(9)
      cfs(4)=aex*cfs(9)
   case(10)
      cfs(1)=aey*aez*cfs(10)
      cfs(3)=aez*cfs(10)
      cfs(4)=aey*cfs(10)
   case(11)
      cfs(1)=aex*aex*aex*cfs(11)
      cfs(2)=c3*aex*aex*cfs(11)
      cfs(5)=c3*aex*cfs(11)
   case(12)
      cfs(1)=aey*aey*aey*cfs(12)
      cfs(3)=c3*aey*aey*cfs(12)
      cfs(6)=c3*aey*cfs(12)
   case(13)
      cfs(1)=aez*aez*aez*cfs(13)
      cfs(4)=c3*aez*aez*cfs(13)
      cfs(7)=c3*aez*cfs(13)
   case(14)
      cfs(1)=aex*aex*aey*cfs(14)
      cfs(2)=c2*aex*aey*cfs(14)
      cfs(3)=aex*aex*cfs(14)
      cfs(5)=aey*cfs(14)
      cfs(8)=c2*aex*cfs(14)
   case(15)
      cfs(1)=aex*aex*aez*cfs(15)
      cfs(2)=c2*aex*aez*cfs(15)
      cfs(4)=aex*aex*cfs(15)
      cfs(5)=aez*cfs(15)
      cfs(9)=c2*aex*cfs(15)
   case(16)
      cfs(1)=aey*aey*aex*cfs(16)
      cfs(2)=aey*aey*cfs(16)
      cfs(3)=c2*aey*aex*cfs(16)
      cfs(6)=aex*cfs(16)
      cfs(8)=c2*aey*cfs(16)
   case(17)
      cfs(1)=aey*aey*aez*cfs(17)
      cfs(3)=c2*aey*aez*cfs(17)
      cfs(4)=aey*aey*cfs(17)
      cfs(6)=aez*cfs(17)
      cfs(10)=c2*aey*cfs(17)
   case(18)
      cfs(1)=aez*aez*aex*cfs(18)
      cfs(2)=aez*aez*cfs(18)
      cfs(4)=c2*aez*aex*cfs(18)
      cfs(7)=aex*cfs(18)
      cfs(9)=c2*aez*cfs(18)
   case(19)
      cfs(1)=aez*aez*aey*cfs(19)
      cfs(3)=aez*aez*cfs(19)
      cfs(4)=c2*aez*aey*cfs(19)
      cfs(7)=aey*cfs(19)
      cfs(10)=c2*aez*cfs(19)
   case(20)
      cfs(1)=aex*aey*aez*cfs(20)
      cfs(2)=aez*aey*cfs(20)
      cfs(3)=aex*aez*cfs(20)
      cfs(4)=aex*aey*cfs(20)
      cfs(8)=aez*cfs(20)
      cfs(9)=aey*cfs(20)
      cfs(10)=aex*cfs(20)
   case default
      continue
   end select

   return
end subroutine rhftce


!! --------------------------------------------------------------[SAW1710]-
!     made explicit, made pure
pure subroutine prod(c,d,s,iff1,iff2)
   implicit none
   integer,intent(in)  :: iff1,iff2
   real(wp), intent(in)  :: c(*),d(*)
   real(wp), intent(out) :: s(*)
   if(iff1.gt.10.or.iff2.gt.10) goto 30
   if(iff1.gt.4.or.iff2.gt.4) goto 20
   ! s-s - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s( 1)=c( 1)*d( 1)
   !              end of s - s
   if(iff1.eq.1.and.iff2.eq.1) return
   ! s-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s( 2)=c( 1)*d( 2)+c( 2)*d( 1)
   s( 3)=c( 1)*d( 3)+c( 3)*d( 1)
   s( 4)=c( 1)*d( 4)+c( 4)*d( 1)
   !              end of s - p
   if(iff1.eq.1.or.iff2.eq.1) return
   ! p-p - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s( 5)=c( 2)*d( 2)
   s( 6)=c( 3)*d( 3)
   s( 7)=c( 4)*d( 4)
   s( 8)=c( 2)*d( 3)+c( 3)*d( 2)
   s( 9)=c( 2)*d( 4)+c( 4)*d( 2)
   s(10)=c( 3)*d( 4)+c( 4)*d( 3)
   !              end of p - p
   return
   20  continue
   ! s-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s( 1)=c( 1)*d( 1)
   s( 2)=c( 1)*d( 2)+c( 2)*d( 1)
   s( 3)=c( 1)*d( 3)+c( 3)*d( 1)
   s( 4)=c( 1)*d( 4)+c( 4)*d( 1)
   s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)
   s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)
   s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)
   s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)
   s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)
   s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)
   !              end of s - d
   if(iff1.eq.1.or.iff2.eq.1) return
   ! p-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s(11)=c( 2)*d( 5)+c( 5)*d( 2)
   s(12)=c( 3)*d( 6)+c( 6)*d( 3)
   s(13)=c( 4)*d( 7)+c( 7)*d( 4)
   s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)
   s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)
   s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)
   s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)
   s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)
   s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)
   s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9) &
      & +c( 9)*d( 3)+c( 4)*d( 8)+c( 8)*d( 4)
   !              end of p - d
   if(iff1.lt.5.or.iff2.lt.5) return
   ! d-d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s(21)=c( 5)*d( 5)
   s(22)=c( 6)*d( 6)
   s(23)=c( 7)*d( 7)
   s(24)=c( 5)*d( 8)+c( 8)*d( 5)
   s(25)=c( 5)*d( 9)+c( 9)*d( 5)
   s(26)=c( 6)*d( 8)+c( 8)*d( 6)
   s(27)=c( 6)*d(10)+c(10)*d( 6)
   s(28)=c( 7)*d( 9)+c( 9)*d( 7)
   s(29)=c( 7)*d(10)+c(10)*d( 7)
   s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)
   s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)
   s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)
   s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)
   s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)
   s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)
   !              end of d - d
   return
   30  continue
   ! s-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s( 1)=c( 1)*d( 1)
   s( 2)=c( 1)*d( 2)+c( 2)*d( 1)
   s( 3)=c( 1)*d( 3)+c( 3)*d( 1)
   s( 4)=c( 1)*d( 4)+c( 4)*d( 1)
   s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)
   s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)
   s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)
   s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)
   s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)
   s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)
   s(11)=c( 2)*d( 5)+c( 5)*d( 2)+c(1)*d(11)+c(11)*d(1)
   s(12)=c( 3)*d( 6)+c( 6)*d( 3)+c(1)*d(12)+c(12)*d(1)
   s(13)=c( 4)*d( 7)+c( 7)*d( 4)+c(1)*d(13)+c(13)*d(1)
   s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)+c(1)*d(14)+ &
      & c(14)*d(1)
   s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)+c(1)*d(15)+ &
      & c(15)*d(1)
   s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)+c(1)*d(16)+ &
      & c(16)*d(1)
   s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)+c(1)*d(17)+ &
      & c(17)*d(1)
   s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)+c(1)*d(18)+ &
      & c(18)*d(1)
   s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)+c(1)*d(19)+ &
      & c(19)*d(1)
   s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9)+c(9)*d(3)+c(4)*d(8)+ &
      & c(8)*d(4)+c(1)*d(20)+c(20)*d(1)
   !              end of s - f
   if(iff1.eq.1.or.iff2.eq.1) return
   ! p-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s(21)=c( 5)*d( 5)+c(2)*d(11)+c(11)*d(2)
   s(22)=c( 6)*d( 6)+c(3)*d(12)+c(12)*d(3)
   s(23)=c( 7)*d( 7)+c(4)*d(13)+c(13)*d(4)
   s(24)=c( 5)*d( 8)+c( 8)*d( 5)+c(3)*d(11)+c(11)*d(3)+c(2)*d(14)+ &
      & c(14)*d(2)
   s(25)=c( 5)*d( 9)+c( 9)*d( 5)+c(2)*d(15)+c(15)*d(2)+c(4)*d(11)+ &
      & c(11)*d(4)
   s(26)=c( 6)*d( 8)+c( 8)*d( 6)+c(2)*d(12)+c(12)*d(2)+c(3)*d(16)+ &
      & c(16)*d(3)
   s(27)=c( 6)*d(10)+c(10)*d( 6)+c(3)*d(17)+c(17)*d(3)+c(4)*d(12)+ &
      & c(12)*d(4)
   s(28)=c( 7)*d( 9)+c( 9)*d( 7)+c(2)*d(13)+c(13)*d(2)+c(4)*d(18)+ &
      & c(18)*d(4)
   s(29)=c( 7)*d(10)+c(10)*d( 7)+c(3)*d(13)+c(13)*d(3)+c(4)*d(19)+ &
      & c(19)*d(4)
   s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)+c(2)*d(16)+c(16)*d(2)+ &
      & c(3)*d(14)+c(14)*d(3)
   s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)+c(2)*d(18)+c(18)*d(2)+ &
      & c(4)*d(15)+c(15)*d(4)
   s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)+c(3)*d(19)+c(19)*d(3)+ &
      & c(4)*d(17)+c(17)*d(4)
   s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)+c(3)*d(15)+ &
      & c(15)*d(3)+c(4)*d(14)+c(14)*d(4)+c(2)*d(20)+c(20)*d(2)
   s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)+c(2)*d(17)+ &
      & d(2)*c(17)+c(3)*d(20)+c(20)*d(3)+c(4)*d(16)+c(16)*d(4)
   s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)+c(2)*d(19)+ &
      & c(19)*d(2)+c(3)*d(18)+c(18)*d(3)+c(4)*d(20)+c(20)*d(4)
   !              end of p - f
   if(iff1.eq.2.or.iff2.eq.2) return
   ! d-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s(36)=c(5)*d(11)+c(11)*d(5)
   s(37)=c(6)*d(12)+c(12)*d(6)
   s(38)=c(7)*d(13)+c(13)*d(7)
   s(39)=c(6)*d(11)+c(11)*d(6)+c(5)*d(16)+c(16)*d(5)+c(8)*d(14)+ &
      & c(14)*d(8)
   s(40)=c(7)*d(11)+c(11)*d(7)+c(5)*d(18)+c(18)*d(5)+c(9)*d(15)+ &
      & c(15)*d(9)
   s(41)=c(5)*d(12)+c(12)*d(5)+c(6)*d(14)+c(14)*d(6)+c(8)*d(16)+ &
      & c(16)*d(8)
   s(42)=c(5)*d(13)+c(13)*d(5)+c(7)*d(15)+c(15)*d(7)+c(9)*d(18)+ &
      & c(18)*d(9)
   s(43)=c(7)*d(12)+c(12)*d(7)+c(6)*d(19)+c(19)*d(6)+c(10)*d(17)+ &
      & c(17)*d(10)
   s(44)=c(6)*d(13)+c(13)*d(6)+c(7)*d(17)+c(17)*d(7)+c(10)*d(19)+ &
      & c(19)*d(10)
   s(45)=c(8)*d(11)+c(11)*d(8)+c(5)*d(14)+c(14)*d(5)
   s(46)=c(9)*d(11)+c(11)*d(9)+c(5)*d(15)+c(15)*d(5)
   s(47)=c(8)*d(12)+c(12)*d(8)+c(6)*d(16)+c(16)*d(6)
   s(48)=c(10)*d(12)+c(12)*d(10)+c(6)*d(17)+c(17)*d(6)
   s(49)=c(10)*d(13)+c(13)*d(10)+c(7)*d(19)+c(19)*d(7)
   s(50)=c(9)*d(13)+c(13)*d(9)+c(7)*d(18)+c(18)*d(7)
   s(51)=c(8)*d(13)+c(13)*d(8)+c(7)*d(20)+c(20)*d(7)+c(9)*d(19)+ &
      & c(19)*d(9)+c(10)*d(18)+c(18)*d(10)
   s(52)=c(10)*d(11)+c(11)*d(10)+c(5)*d(20)+c(20)*d(5)+c(9)*d(14)+ &
      & c(14)*d(9)+c(8)*d(15)+c(15)*d(8)
   s(53)=c(9)*d(12)+c(12)*d(9)+c(6)*d(20)+c(20)*d(6)+c(10)*d(16)+ &
      & c(16)*d(10)+c(8)*d(17)+c(17)*d(8)
   s(54)=c(5)*d(17)+c(17)*d(5)+c(6)*d(15)+c(15)*d(6)+c(14)*d(10)+ &
      & d(14)*c(10)+c(9)*d(16)+c(16)*d(9)+c(8)*d(20)+c(20)*d(8)
   s(55)=c(5)*d(19)+c(19)*d(5)+c(7)*d(14)+c(14)*d(7)+c(10)*d(15)+ &
      & c(15)*d(10)+c(9)*d(20)+c(20)*d(9)+c(8)*d(18)+c(18)*d(8)
   s(56)=c(6)*d(18)+c(18)*d(6)+c(7)*d(16)+c(16)*d(7)+c(10)*d(20)+ &
      & c(20)*d(10)+c(9)*d(17)+c(17)*d(9)+c(8)*d(19)+c(19)*d(8)
   !              end of d - f
   if(iff1.eq.3.or.iff2.eq.3) return
   ! f-f - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1710
   s(57)=c(11)*d(11)
   s(58)=c(12)*d(12)
   s(59)=c(13)*d(13)
   s(60)=c(11)*d(12)+c(12)*d(11)+c(14)*d(16)+c(16)*d(14)
   s(61)=c(11)*d(13)+c(13)*d(11)+c(15)*d(18)+c(18)*d(15)
   s(62)=c(12)*d(13)+c(13)*d(12)+c(17)*d(19)+c(19)*d(17)
   s(63)=c(11)*d(14)+c(14)*d(11)
   s(64)=c(11)*d(15)+c(15)*d(11)
   s(65)=c(13)*d(18)+c(18)*d(13)
   s(66)=c(13)*d(19)+c(19)*d(13)
   s(67)=c(12)*d(17)+d(12)*c(17)
   s(68)=c(12)*d(16)+c(16)*d(12)
   s(69)=c(11)*d(16)+c(16)*d(11)+c(14)*d(14)
   s(70)=c(11)*d(18)+c(18)*d(11)+c(15)*d(15)
   s(71)=c(13)*d(15)+c(15)*d(13)+c(18)*d(18)
   s(72)=c(13)*d(17)+c(17)*d(13)+c(19)*d(19)
   s(73)=c(12)*d(14)+c(14)*d(12)+c(16)*d(16)
   s(74)=c(12)*d(19)+c(19)*d(12)+c(17)*d(17)
   s(75)=c(11)*d(17)+c(17)*d(11)+c(14)*d(20)+c(20)*d(14)+ &
      & c(15)*d(16)+c(16)*d(15)
   s(76)=c(11)*d(19)+c(19)*d(11)+c(20)*d(15)+c(15)*d(20)+ &
      & c(14)*d(18)+c(18)*d(14)
   s(77)=c(12)*d(18)+c(18)*d(12)+c(16)*d(19)+c(19)*d(16)+ &
      & c(17)*d(20)+c(20)*d(17)
   s(78)=c(13)*d(14)+c(14)*d(13)+c(15)*d(19)+c(19)*d(15)+ &
      & c(18)*d(20)+c(20)*d(18)
   s(79)=c(12)*d(15)+c(15)*d(12)+c(14)*d(17)+c(17)*d(14)+ &
      & c(16)*d(20)+c(20)*d(16)
   s(80)=c(13)*d(16)+c(16)*d(13)+c(17)*d(18)+c(18)*d(17)+ &
      & c(19)*d(20)+c(20)*d(19)
   s(81)=c(11)*d(20)+c(20)*d(11)+c(14)*d(15)+c(15)*d(14)
   s(82)=c(12)*d(20)+c(20)*d(12)+c(16)*d(17)+c(17)*d(16)
   s(83)=c(13)*d(20)+c(20)*d(13)+c(18)*d(19)+c(19)*d(18)
   s(84)=c(14)*d(19)+c(19)*d(14)+c(15)*d(17)+c(17)*d(15)+ &
      & c(16)*d(18)+c(18)*d(16)+c(20)*d(20)
   !              end of f - f
   return
end subroutine prod


!! --------------------------------------------------------------[SAW1710]-
!     made pure, made explicit
pure subroutine lmnpre(l,m,n,lmnexp,lmnfak)
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(out) :: lmnfak
   integer,intent(out) :: lmnexp
   integer :: lh,mh,nh
   real(wp), parameter :: dftr(7) = &
      & (/1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   lmnfak=0
   lmnexp=0
   if (mod(l,2).ne.0) return
   lh=l/2
   if (mod(m,2).ne.0) return
   mh=m/2
   if (mod(n,2).ne.0) return
   nh=n/2
   lmnexp=lh+mh+nh
   lmnfak=dftr(lh+1)*dftr(mh+1)*dftr(nh+1)
   return
end subroutine lmnpre


!! --------------------------------------------------------------[SAW1710]-
!     made explicit, changed data in parameter, made elemental
elemental function olap2(l,m,n,arg,gama)
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(in)  :: arg,gama
   real(wp)  :: olap2
   integer :: lh,mh,nh
   real(wp), parameter :: dftr(7) = &
      & (/1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   if (mod(l,2).ne.0) goto 50
   lh=l/2
   if (mod(m,2).ne.0) goto 50
   mh=m/2
   if (mod(n,2).ne.0) goto 50
   nh=n/2
   olap2=arg*gama**(lh+mh+nh)*dftr(lh+1)*dftr(mh+1)*dftr(nh+1)
   return
   50  olap2=0.d0
   return
end function olap2

!! --------------------------------------------------------------[SAW1710]-
!     removed implicit, changed data in parameter, made elemental
elemental function olap(l,m,n,gama)
   implicit none
   integer,intent(in) :: l,m,n
   real(wp), intent(in) :: gama
   real(wp)  :: olap
   real(wp)  :: gm
   integer :: lh,mh,nh
   real(wp), parameter :: dftr(7) = &
      & (/1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   real(wp), parameter :: pi = 3.1415926535897932384626433832795029d0
   real(wp), parameter :: thehlf = 0.5d0
   gm=1.0d0/gama
   if (mod(l,2).ne.0) goto 50
   lh=l/2
   if (mod(m,2).ne.0) goto 50
   mh=m/2
   if (mod(n,2).ne.0) goto 50
   nh=n/2
   olap=(dsqrt(pi*gm))**3*(thehlf*gm)**(lh+mh+nh)* &
      & dftr(lh+1)*dftr(mh+1)*dftr(nh+1)
   return
   50  olap=0.d0
   return
end function olap

!! --------------------------------------------------------------[SAW1710]-
!     removed implicit, made pure
pure subroutine opab1(l,m,n,ga,v,d)
   !           electronic part of dipole moment
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(in)  :: ga,d(*)
   real(wp), intent(out) :: v(*)
   integer :: i
   real(wp)  :: g(4)
   g(1)=olap(l,m,n,ga)
   g(2)=olap(l+1,m,n,ga)
   g(3)=olap(l,m+1,n,ga)
   g(4)=olap(l,m,n+1,ga)
   v(1)=g(2)+d(1)*g(1)
   v(2)=g(3)+d(2)*g(1)
   v(3)=g(4)+d(3)*g(1)
   do i=1,3
      v(i)=-v(i)
   enddo
   return
end subroutine opab1

!! --------------------------------------------------------------[SAW1710]-
!     made explicit, made pure
pure subroutine opab4(l,m,n,ga,v,d)
   !           electronic part of second moment
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(in)  :: ga,d(*)
   real(wp), intent(out) :: v(*)
   real(wp)  :: g(10)
   integer :: i
   g(1)=olap(l,m,n,ga)
   g(2)=olap(l+1,m,n,ga)
   g(3)=olap(l,m+1,n,ga)
   g(4)=olap(l,m,n+1,ga)
   g(5)=olap(l+2,m,n,ga)
   g(6)=olap(l,m+2,n,ga)
   g(7)=olap(l,m,n+2,ga)
   g(8)=olap(l+1,m+1,n,ga)
   g(9)=olap(l+1,m,n+1,ga)
   g(10)=olap(l,m+1,n+1,ga)
   v(1)=g(5)+2.d0*d(1)*g(2)+d(1)**2*g(1)
   v(2)=g(6)+2.d0*d(2)*g(3)+d(2)**2*g(1)
   v(3)=g(7)+2.d0*d(3)*g(4)+d(3)**2*g(1)
   v(4)=g(8)+d(1)*g(3)+d(2)*g(2)+d(1)*d(2)*g(1)
   v(5)=g(9)+d(1)*g(4)+d(3)*g(2)+d(1)*d(3)*g(1)
   v(6)=g(10)+d(2)*g(4)+d(3)*g(3)+d(2)*d(3)*g(1)
   do i=1,6
      v(i)=-v(i)
   enddo
   return
end subroutine opab4

!! --------------------------------------------------------------[SAW1710]-
!     made explict, made pure
pure subroutine opac3(l,m,n,ga,v,d)
   !           electronic part of charge density
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(in)  :: ga,d(1)
   real(wp), intent(out) :: v(1)
   real(wp)  :: dd(3),ddsq
   integer :: i
   ddsq=0.d0
   do i=1,3
      dd(i)=-d(i)
      ddsq=ddsq+d(i)**2
   enddo
   v(1) =1.d0
   if(l.ne.0) v(1)=dd(1)**l
   if(m.ne.0) v(1)=v(1)*dd(2)**m
   if(n.ne.0) v(1)=v(1)*dd(3)**n
   v(1)=v(1)*dexp(-ga*ddsq)
   return
end subroutine opac3

!! --------------------------------------------------------------[SAW1710]-
!     made explicit, made pure
pure subroutine opad1(l,m,n,ga,v,d)
   !           electronic part of overlap
   implicit none
   integer,intent(in)  :: l,m,n
   real(wp), intent(in)  :: ga,d(1)
   real(wp), intent(out) :: v(1)
   v(1)=olap(l,m,n,ga)
   return
end subroutine opad1


!! --------------------------------------------------------------[SAW1710]-
!     changed do loops
subroutine opap4(l,m,n,gc,v,rc)
   !
   !     this subroutine calculates the expectation-values of
   !     -p**2/2 (kinetic energy) and p**4/(8*cl**2) (relativistic
   !     correction of the kinetic energy)
   !
   !     the program works correct only for s and p and the irreducible
   !     linear combinations of the primitive d,f,... functions
   !
   !     written by s. brode in april,1984
   !
   implicit real(wp) (a-h,o-z)
   logical :: lodd,modd,nodd,leven,meven,neven
   real(wp)  :: v(1),rc(3),rca(3),rcb(3),srcab(3),o(25)
   ! ----- ra,rb centres ga,gb exponents ia,ib monom indices
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
   real(wp), parameter :: cl = 137.03604d0
   real(wp), parameter :: dftr(7) = &
      & (/1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   real(wp), parameter :: pi=3.1415926535897932384626433832795029d0
   real(wp), parameter :: a0=0.0d0,a1=1.0d0,a2=2.0d0,a4=4.0d0,a8=8.0d0
   integer,parameter :: lmni(5) = (/0,3*1,6*2,10*3,15*4/)
   !
   !     define the overlap-integral (in line function definition)
   !
   !!!ola(ll,mm,nn)=gch**(ll/2+mm/2+nn/2)*dftr(ll/2+1)*dftr(mm/2+1)* &
   !!!     & dftr(nn/2+1)
   !
   !     some previous calculations
   !
   !     do 10 i=1,25
   !  10 o(i)=a0
   o(1:7)=0
   gch=a1/(a2*gc)
   ropigc=dsqrt(pi/gc)*pi/gc
   lodd=mod(l,2).eq.1
   leven=.not.lodd
   modd=mod(m,2).eq.1
   meven=.not.modd
   nodd=mod(n,2).eq.1
   neven=.not.nodd
   ga2=ga*a2
   gb2=gb*a2
   gab4=ga2*gb2
   prcaa=a0
   prcbb=a0
   do i=1,3
      rca(i)=rc(i)-ra(i)
      rcb(i)=rc(i)-rb(i)
      srcab(i)=rca(i)+rcb(i)
      prcaa=prcaa+rca(i)*rca(i)
      prcbb=prcbb+rcb(i)*rcb(i)
   enddo
   zwa=ga2*prcaa-dble(2*lmni(ia)+3)
   zwb=gb2*prcbb-dble(2*lmni(ib)+3)
   !
   !     calculation of the overlap-integrals
   !
   if(lodd.or.modd.or.nodd) go to 110
   o( 1)=ola(l   ,m  ,n  )
   o( 5)=ola(l+2,m   ,n  )
   o( 6)=ola(l   ,m+2,n  )
   o( 7)=ola(l   ,m  ,n+2)
   !     o(11)=ola(l+4,m   ,n  )
   !     o(12)=ola(l   ,m+4,n  )
   !     o(13)=ola(l   ,m  ,n+4)
   !     o(20)=ola(l+2,m+2,n   )
   !     o(21)=ola(l+2,m   ,n+2)
   !     o(22)=ola(l   ,m+2,n+2)
   110 continue
   !
   if(leven.or.modd.or.nodd) go to 120
   o( 2)=ola(l+1,m   ,n  )
   !     o( 8)=ola(l+3,m   ,n  )
   !     o(14)=ola(l+1,m+2,n   )
   !     o(15)=ola(l+1,m   ,n+2)
   120 continue
   !
   if(lodd.or.meven.or.nodd) go to 130
   o( 3)=ola(l   ,m+1,n  )
   !     o( 9)=ola(l   ,m+3,n  )
   !     o(16)=ola(l+2,m+1,n   )
   !     o(17)=ola(l   ,m+1,n+2)
   130 continue
   !
   if(lodd.or.modd.or.neven) go to 140
   o( 4)=ola(l   ,m  ,n+1)
   !     o(10)=ola(l   ,m  ,n+3)
   !     o(18)=ola(l+2,m   ,n+1)
   !     o(19)=ola(l   ,m+2,n+1)
   140 continue
   !
   !     if(leven.or.meven.or.nodd) go to 150
   !     o(23)=ola(l+1,m+1,n   )
   ! 150 continue
   !     if(leven.or.modd.or.neven) go to 160
   !     o(24)=ola(l+1,m   ,n+1)
   ! 160 continue
   !     if(lodd.or.meven.or.neven) go to 170
   !     o(25)=ola(l   ,m+1,n+1)
   ! 170 continue
   !
   !     calculation of kinetic energy
   !
   p2=-ga*(zwa*o(1) &
      & +ga2*(a2*(rca(1)*o(2)+rca(2)*o(3)+rca(3)*o(4)) &
      & +o(5)+o(6)+o(7)))


   !     calculation of relativistic correction

   !     p4=(gab4*(o(11)+o(12)+o(13)+a2*(o(20)+o(21)+o(22))
   !    .         +a2*(srcab(1)*(o( 8)+o(16)+o(18))
   !    .             +srcab(2)*(o(14)+o( 9)+o(19))
   !    .             +srcab(3)*(o(15)+o(17)+o(10)))
   !    .         +a4*(rca(1)*rcb(1)*o(5)
   !    .             +rca(2)*rcb(2)*o(6)
   !    .             +rca(3)*rcb(3)*o(7)
   !    .             +a2*((rca(1)*rcb(2)+rca(2)+rcb(1))*o(23)
   !    .                 +(rca(1)*rcb(3)+rca(3)+rcb(1))*o(24)
   !    .                 +(rca(2)*rcb(3)+rca(3)+rcb(2))*o(25))))
   !    .   +(gb2*zwa+ga2*zwb)*(o(5)+o(6)+o(7))
   !    .   +a2*(gb2*zwa*(rcb(1)*o(2)+rcb(2)*o(3)+rcb(3)*o(4))
   !    .       +ga2*zwb*(rca(1)*o(2)+rca(2)*o(3)+rca(3)*o(4)))
   !    .   +zwa*zwb*o(1))
   !    .   *(-gab4)/(a8*cl**2)

   !     store p2 and p4 to the correct storage-locations

   v(1)=p2*ropigc
   !     v(2)=p4*ropigc
   !     v(3)=v(1)+v(2)
   return

contains

   real(wp) function ola(ll,mm,nn)
      integer, intent(in) :: ll,mm,nn
      ola=gch**(ll/2+mm/2+nn/2)*dftr(ll/2+1)*dftr(mm/2+1)*dftr(nn/2+1)
   end function ola

end subroutine opap4

subroutine opam(l,m,n,gama,v,d)
   implicit real(wp)(a-h,o-z)
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
   common /gf/f(7),g(8) !to share with bip
   real(wp) pa(3),pb(3),ovl(8),d(3),v(3)

   integer,parameter :: lx(20) = &
   &(/1,2,1,1,3,1,1,2,2,1,4,1,1,3,3,2,1,2,1,2/)
   integer,parameter :: my(20) = &
   &(/1,1,2,1,1,3,1,2,1,2,1,4,1,2,1,3,3,1,2,2/)
   integer,parameter :: nz(20) = &
   &(/1,1,1,2,1,1,3,1,2,2,1,1,4,1,2,1,2,3,3,2/)
   real(wp),parameter :: pi=3.1415926535897932384626433832795029d0
   integer,parameter :: lmni(20) = &
   &(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)
   real(wp),parameter :: dftr(7) = &
   &(/1.0D0,1.0D0,3.0D0,15.0D0,105.0D0,945.0D0,10395.0D0/)
   gch=0.5D0/gama
   lmn=lmni(ia)+lmni(ib)
   pre=(pi/gama)**1.5D0
   do i=1,lmn
      ovl(i)=0.D0
   enddo
   do i=1,lmn,2
      i1=(i+1)/2
      ovl(i)=dftr(i1)*gch**(i1-1)
   enddo
   do i=1,3
      pa(i)=d(i)-ra(i)
      pb(i)=d(i)-rb(i)
   enddo
      s1=0.D0
      s2=0.D0
      s3=0.D0
      s4=0.D0
      s5=0.D0
      s6=0.D0
      s7=0.D0
      s8=0.D0
      s9=0.D0
   call bip(lx(ia),lx(ib),pa(1),pb(1))
      l=lx(ia)+lx(ib)
      l1=l-1
    do i=1,l1
      s1=s1+f(i)*ovl(i+1)
      s2=s2+f(i)*ovl(i)
    enddo
    do i=1,l
      s3=s3+g(i)*ovl(i)
    enddo
      call bip(my(ia),my(ib),pa(2),pb(2))
      m=my(ia)+my(ib)
      m1=m-1
    do i=1,m1
      s4=s4+f(i)*ovl(i+1)
      s5=s5+f(I)*ovl(i)
    enddo
    do i=1,m
      s6 =s6+g(i)*ovl(i)
    enddo
    call bip(nz(ia),nz(ib),pa(3),pb(3))
      n=nz(ia)+nz(ib)
      n1=n-1
    do i=1,n1
      s7=s7+f(i)*ovl(i+1)
      s8=s8+f(i)*ovl(i)
    enddo
    do i=1,n
      s9=s9+g(i)*ovl(i)
    enddo
      sx=s1+d(1)*s2
      sy=s4+d(2)*s5
      sz=s7+d(3)*s8
      v(1)=pre*s2*(sy*s9-sz*s6)
      v(2)=pre*s5*(sz*s3-sx*s9)
      v(3)=pre*s8*(sx*s6-sy*s3)
    return
end subroutine opam

subroutine bip(la,lb,a,b)
   implicit real(wp)(a-h,o-z)
   common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
   common /gf/f(7),g(8)
   real(wp),parameter :: c1=1.D0,c2=2.D0,c3=3.D0,c4=4.D0,c6=6.D0,c8=8.D0
   f(la+lb-1)=1.D0
   gb2=gb+gb
   g(la+lb)=-gb2
      go to(1,10,20,30),la
    1 go to (5,6,7,8),lb
!        s - s
    5 g(1)=-b*gb2
      return
!        s - p
    6 f(1)=b
      g(1)=c1-b*b*gb2
      g(2)=-c4*b*gb
      return
!       s - d
    7 f(1)=b*b
      f(2)=b+b
      g(1)=c2*b*(c1-b*b*gb)
      g(2)=c2-c6*gb*b*b
      g(3)=-c6*b*gb
      return
!       s -f
    8 f(1)=b**3
      f(2)=c3*b*b
      f(3)=c3*b
      g(1)=b*b*(c3-b*b*gb)
      g(2)=(b+b)*(c3-c2*b*b*gb2)
      g(3)=c3*(c1-c4*gb*b*b)
      g(4)=-c8*b*gb
      return
   10 go to (11,12,13,14),lb
!       p - s
   11 f(1)=a
      g(1)=-a*b*gb2
      g(2)=-gb2*(a+b)
      return
!       p - p
   12 f(1)=a*b
      f(2)=a+b
      g(1)=a*(c1-gb2*b*b)
      g(2)=c1-b*gb2*(a+a+b)
      g(3)=-gb2*(a+b+b)
      return
!       p - d
   13 f(1)=a*b*b
      f(2)=b*(a+a+b)
      f(3)=a+b+b
      g(1)=c2*a*b*(c1-gb*b*b)
      g(2)=c2*(a+b)*(c1-gb*b*b)
      g(3)=c2-c3*gb2*b*(a+b)
      g(4)=-gb2*(a+b*c3)
      return
!       p - f
   14 f(1)=a*b**3
      f(2)=b*b*(b+c3*a)
      f(3)=c3*b*(a+b)
      f(4)=a+c3*b
      g(1)=a*b*b*(c3-gb2*b*b)
      g(2)=c3*b*(c2*a+b)-gb2*(b**3)*(c4*a+b)
      g(3)=c3*(a+c2*b)-c3*gb2*b*b*(c2*b+c3*a)
      g(4)=c3-c4*b*gb*(c2*a+c3*b)
      g(5)=-gb2*(a+c4*b)
      return
   20 go to (21,22,23,24),lb
!       d - s
   21 f(1)=a*a
      f(2)=a+a
      g(1)=-gb2*a*a*b
      g(2)=-gb2*a*(a+b+b)
      g(3)=-gb2*(a+a+b)
      return
!       d - p
   22 f(1)=a*a*b
      f(2)=a*(a+b+b)
      f(3)=(a+a+b)
      g(1)=a*a*(c1-gb2*b*b)
      g(2)=c2*a*(c1-gb2*b*(a+b))
      g(3)=c1-gb2*(a*a+b*b+c4*a*b)
      g(4)=-c2*gb2*(a+b)
      return
!       d - d
   23 f(1)=(a*b)**2
      f(2)=c2*a*b*(a+b)
      f(3)=a*a+c4*a*b+b*b
      f(4)=c2*(a+b)
      g(1)=c2*a*a*b*(c1-gb2*b*b)
      g(2)=c2*a*((a+b+b)-gb*b*b*(c3*a+c2*b))
      g(3)=c2*(a+a+b)-b*gb2*(b*(b+c6*a)+c3*a*a)
      g(4)=c2-gb2*(a*(a+c6*b)+c3*b*b)
      g(5)=-gb2*(c3*b+c2*a)
      return
!       d - f
   24 f(1)=a*a*(b**3)
      f(2)=a*b*b*(c3*a+c2*b)
      f(3)=b*(b*b+c3*a*(b+b+a))
      f(4)=c3*b*(b+a+a)+a*a
      f(5)=c3*b+c2*a
      g(1)=(c3-gb2*b*b)*(a*b)**2
      g(2)=a*b*(c6*(a+b)-c4*gb*b*b*(b+a+a))
      g(3)=c3*(a*a+c4*a*b+b*b)-gb*b*b*(b*b+c8*a*b+c6*a*a)
      g(4)=c6*(a+b)-c4*b*gb2*(b*b+c3*a*b+a*a)
      g(5)=c3-gb2*(a*a+c8*a*b+b*b)
      g(6)=-gb2*c2*(a+b+b)
      return
   30 go to (31,32,33,34),lb
!       f - s
   31 f(1)=a**3
      f(2)=c3*a*a
      f(3)=c3*a
      g(1)=-gb2*(a**3)*b
      g(2)=-gb2*a*a*(a+c3*b)
      g(3)=-c6*a*gb*(a+b)
      g(4)=-gb2*(c3*a+b)
      return
!       f - p
   32 f(1)=b*a**3
      f(2)=a*a*(a+c3*b)
      f(3)=c3*a*(a+b)
      f(4)=c3*a+b
      g(1)=(c1-gb2*b*b)*a**3
      g(2)=a*a*(c3-b*gb2*(c3*b+a))
      g(3)=a*(c3-gb2*(a*a+c6*a*b+c3*b*b))
      g(4)=c1-gb2*(c3*a*a+c6*a*b+b*b)
      g(5)=-gb2*(c3*a+b)
      return
!       f - d
   33 f(1)=a*(a*b)**2
      f(2)=a*a*b*(c3*a+c2*a)
      f(3)=a*(a*a+c6*a*b+c3*b*b)
      f(4)=c3*a*a+b*(c6*a+b)
      f(5)=c3*a+c2*b
      g(1)=c2*(a**3)*b*(c1-gb*b*b)
      g(2)=c2*a*a*((a+c3*b)-c3*gb*b*b*(a+b))
      g(3)=c6*a*((a+b)-gb*b*(a*a+c3*a*b+b*b))
      g(4)=c2*((c3*a+b)-gb*(a*a*(a+9.d0*b)+b*b*(9.d0*a+b)))
      g(5)=c2-c3*gb2*(a*(a+c3*b)+b*b)
      g(6)=-c3*gb2*(a+b)
      return
!       f - f
   34 f(1)=(a*b)**3
      f(2)=c3*((a*b)**2)*(a+b)
      f(3)=c3*a*b*(a*(a+c3*b)+b*b)
      f(4)=a*a*(a+9.d0*b)+b*b*(b+9.d0*a)
      f(5)=c3*(a*(a+c3*b)+b*b)
      f(6)=c3*(a+b)
      g(1)=a*((a*b)**2)*(c3-gb2*b*b)
      g(2)=a*a*b*((c3*b+c2*a)*c3-gb2*b*b*(c4*a+c3*b))
      g(3)=c3*a*((a*(a+c6*b)+c3*b*b)-c4*gb*b*b*(c2*a*a+b*(c4*a+b)))
      g(4)=c3*(c3*a*a+b*(c6*a+b))-gb2*b*b*(b*(b+12.d0*a)+18.d0*a*a)
      g(5)=9.d0*a+c6*b-gb2*(a*a*(a+12.d0*b)+b*b*(18.d0*a+4.d0*b))
      g(6)=c3-c3*gb2*(a*(a+c4*b)+c2*b*b)
      g(7)=-gb2*(c3*a+c4*b)
      return
end subroutine bip


end module intpack
