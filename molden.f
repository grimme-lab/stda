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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! write an output molden file with the normalized AO basis set
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine molden_file(ncent,nprims,nmo,icdim,nbf,imethod,cc,
     .                       ccspin)
      use stdacommon
      implicit none
      integer:: i,j,prim,at,ncent,nprims,info(nprims)
      integer:: nmo,nbf,imethod,icdim,ccspin(nmo),counter_BF
      integer:: f_info(nbf),flag,flag2
      real*8:: cc(icdim),norm(nprims)

      open(unit=11,file='molden.molden')
      open(unit=12,file='fnorm')
      ! read normalization
      Do i=1,nprims
      read(12,*)norm(i)
      enddo
      close(12,status='delete')
      ! Geometry
      write(11,*) '[Molden Format]'
      write(11,*) '[Atoms] Angs'
      Do i=1, ncent
      write(11,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo
      !Basis set (sp are seperated)

      ! Counting the number of prim in each contraction
      prim=1
      info=1
      Do i=2,nprims
      if(prim/=i)then
      if(ipao(i)==ipao(i-1))then
      info(prim)=info(prim)+1
      info(i)=0
      else
      prim=i
      endif
      endif
      enddo

      write(11,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(11,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(11,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(11,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(11,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(11,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(11,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(11,23)'f',info(i),'1.000000'

      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif


      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(11,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(11,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(11,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(11,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(11,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)
      at=ipat(i)
      if(at/=ipat(i+1)) write(11,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo

      ! MOs

      write(11,*) '[MO]'
      Do i=1,nmo
      write(11,*) 'Sym= X'
      write(11,*) 'Ene=',eps(i)
      if(imethod==1)then
      write(11,*) 'Spin= Alpha'
      else
      if(ccspin(i)==1) write(11,*) 'Spin= Alpha'
      if(ccspin(i)==2) write(11,*) 'Spin= Beta'
      endif
      write(11,*) 'Occup=',int(occ(i))
      Do j=1,nbf

      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then

      if(f_info(j)==14)write(11,*) j, cc((i-1)*nbf+j+2) ! 16
      if(f_info(j)==15)write(11,*) j, cc((i-1)*nbf+j-1) ! 14
      if(f_info(j)==16)write(11,*) j, cc((i-1)*nbf+j-1) ! 15
      if(f_info(j)==17)write(11,*) j, cc((i-1)*nbf+j+1) ! 18
      if(f_info(j)==18)write(11,*) j, cc((i-1)*nbf+j+1) ! 19
      if(f_info(j)==19)write(11,*) j, cc((i-1)*nbf+j-2) ! 17

      else
      write(11,*) j, cc((i-1)*nbf+j)
      endif
      enddo
      enddo
      write(11,*)
      close(11)


 21   format(a,2i7,3f16.8)
 22   format(i7,3x,a)
 23   format(a,3x,i7,3x,a)
 24   format(2f16.8)
 25   format(3f16.8)


      end
