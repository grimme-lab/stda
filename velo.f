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
      SUBROUTINE SETETA(J,MNL)
      use stdacommon
      IMPLICIT REAL*8(A-H,O-Z)                                                  

      common /carte/   lmn(0:3,0:3,0:3)
      dimension mnl(3)

      ity=lmn(mnl(1),mnl(2),mnl(3))
      eta(j,5)=float(ity)
      do k=6,25
         eta(j,k)=0.0d0
      enddo
      eta(j,5+ity)=1.00d0
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c velocity dipole integrals (antisymmetric i.e. <i|j>=-<j|i>)
c <alp1,l1,m1,n1|d/dx|alp2,l2,m2,n2>=-2*alp2*<alp1,l1,m1,n1|alp2,l2+1,m2,n2>
c                                    +    l2*<alp1,l1,m1,n1|alp2,l2-1,m2,n2>
c and so on ...
c
c s.grimme, dec. 1995
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE VELO(I,J,D)
      use stdacommon
      use intpack
      IMPLICIT REAL*8(A-H,O-Z)                                                  

      common /carte/   lmn(0:3,0:3,0:3)

      dimension v(3),point(3),nml(3),d(*)
      
      point=0
      ity=ipty(j)
      if(ity.gt.10.or.ipty(i).gt.10) then
         d(1:3)=0
         return
      endif

      do l=0,2
         do m=0,2
            do n=0,2
               if(ity.eq.lmn(l,m,n)) then
                  nml(1)=l
                  nml(2)=m
                  nml(3)=n
                  goto 100
               endif
            enddo
         enddo
      enddo
      
100   alp=exip(j)
      if(ity.gt.10) goto 99 


C first term dipole

c     point(1)=co(ipat(j),1)
c     point(2)=co(ipat(j),2)
c     point(3)=co(ipat(j),3)
c     call propa(opab1,i,j,v,point,3,dummy)
c     do k=1,3
c        d(k)=2.0*alp*v(k)
c     enddo

c     do 2 k=1,3
c        if(nml(k).eq.0) goto 2
c        nml(k)=nml(k)-1
c        call seteta(j,nml)
c        call propa(opad1,i,j,v,point,1,dummy)
c        nml(k)=nml(k)+1
c        call seteta(j,nml)
c        d(k)=d(k)+nml(k)*v(1)
c 2   continue
c     write(*,*)'typ j',ity,nml
c     write(*,'(''di in  velo'',3f20.12)') d(1),d(2),d(3)

C all terms overlap

      do k=1,3
         nml(k)=nml(k)+1
         call seteta(j,nml)
         !call propa(opad1,i,j,v,point,1,dummy)
         call propa(opad1,point,v,1,i,j)
         nml(k)=nml(k)-1
         call seteta(j,nml)
         d(k)=-2.0*alp*v(1)
      enddo

      do 3 k=1,3
         if(nml(k).eq.0) goto 3
         nml(k)=nml(k)-1
         call seteta(j,nml)
         !call propa(opad1,i,j,v,point,1,dummy)
         call propa(opad1,point,v,1,i,j)
         nml(k)=nml(k)+1
         call seteta(j,nml)
         d(k)=d(k)+nml(k)*v(1)
  3   continue

      return
       
 99   continue
c     write(*,*) 'no f-functions for velocity integrals'
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C s.grimme, april 1998
C checked against OPAM from chandra
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE LXYZ(I,J,D)
      use stdacommon
      use intpack
      IMPLICIT REAL*8(A-H,O-Z)                                                  

      common /carte/   lmn(0:3,0:3,0:3)
      
      dimension v(3),point(3),nml(3),d(*)
      
      ity=ipty(j)
      do l=0,2
         do m=0,2
            do n=0,2
               if(ity.eq.lmn(l,m,n)) then
                  nml(1)=l
                  nml(2)=m
                  nml(3)=n
                  goto 100
               endif
            enddo
         enddo
      enddo
      
100   alp=exip(j)
      if(ity.gt.10) goto 99 

C all terms overlap

      do k=1,3

      vz1=0
      vz3=0

      if(k.eq.3) then
         ii=1
         jj=2
      endif

      if(k.eq.2) then
         ii=1
         jj=3
      endif

      if(k.eq.1) then
         ii=2
         jj=3
      endif

      if(nml(jj).gt.0) then
         nml(jj)=nml(jj)-1
         call seteta(j,nml)
         !call propa(opab1,i,j,v,point,3,dummy)
         call propa(opab1,point,v,3,i,j)
         nml(jj)=nml(jj)+1
         call seteta(j,nml)
         vz1=nml(jj)*v(ii)
      endif

      nml(jj)=nml(jj)+1
      call seteta(j,nml)
      !call propa(opab1,i,j,v,point,3,dummy)
      call propa(opab1,point,v,3,i,j)
      nml(jj)=nml(jj)-1
      call seteta(j,nml)
      vz2=2.*alp*v(ii)
   
      if(nml(ii).gt.0) then
         nml(ii)=nml(ii)-1
         call seteta(j,nml)
         !call propa(opab1,i,j,v,point,3,dummy)
         call propa(opab1,point,v,3,i,j)
         nml(ii)=nml(ii)+1
         call seteta(j,nml)
         vz3=nml(ii)*v(jj)
      endif

      nml(ii)=nml(ii)+1
      call seteta(j,nml)
      !call propa(opab1,i,j,v,point,3,dummy)
      call propa(opab1,point,v,3,i,j)
      nml(ii)=nml(ii)-1
      call seteta(j,nml)
      vz4=2.*alp*v(jj)

      d(k)=vz1-vz2-vz3+vz4

      enddo

      return
       
 99   write(*,*) 'no f-functions for lxyz integrals'
      stop
      end
