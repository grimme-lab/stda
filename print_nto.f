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

! written by Marc de Wegifosse 2019
! need a cleanup ...
      subroutine print_nto_rpa(uci,hci,ca,moci,nci,nroot,nao,iconf,
     .                         maxconf,no,nv)
      implicit none 
      integer:: moci,nci,nroot,nao,maxconf,no,nv
      real*4:: uci(nci,nroot),hci(nci,nroot)
      real*4:: newuci(nci,nroot)
      real*8:: ca(nao*moci)
      integer ::iconf(maxconf,2)
      
      newuci=uci+hci
      
      call print_nto(newuci,ca,moci,nci,nroot,nao,iconf,maxconf,no,nv)
      
      end
      
      subroutine print_nto(Xci,ca,moci,nci,nroot,nao,iconf,maxconf,
     .                     no,nv)
      use stdacommon
      use commonresp
      implicit none
      integer:: moci,nci,nroot,nao,maxconf,no,nv
      real*4:: Xci(nci,nroot),X(no,nv),MO(nao)
      real*8:: ca(nao*moci)
      integer:: iconf(maxconf,2)
      
      integer:: i,j,k,l
      
      real*4:: SIGMA(no),U(no,no),V(nv,nv)
      integer:: iwork(8*no), infoo, lwork
      real*4,allocatable :: work(:)
      
      real*4:: PR_nto
      
      character(len=14):: fname
      
      integer:: ncent,nprims,nmo,nbf,imethod,icdim
      real*8,allocatable:: norm(:)
      
      integer:: prim,at
      integer,allocatable:: info(:),f_info(:)
      integer:: counter_BF
      integer:: flag,flag2,counter
      
      open(unit=11,file='NTOvar')
      read(11,*)ncent,nprims,nmo,icdim,nbf
      close(11,status='delete')      
      allocate(norm(nprims))
      open(unit=12,file='fnorm')
      ! read normalization
      Do i=1,nprims
      read(12,*)norm(i)
      enddo
      close(12,status='delete')
      open(unit=11,file='NTOao')
      allocate(ipat(nprims),ipty(nprims),ipao(nprims),info(nprims))
      allocate(atnam(ncent),exip(nprims),cxip(nprims),co(ncent,4))
      allocate(f_info(nbf))
      Do i=1, nprims
      read(11,121)ipat(i),ipty(i),ipao(i),exip(i),cxip(i)
      enddo
      close(11,status='delete')
      open(unit=11,file='NTOat')
      Do i=1,ncent
      read(11,*)atnam(i),co(i,1:4)
      enddo
      close(11,status='delete')

      
 121  format(3i10,3x,2f16.9)     
      
      
      
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in print NTOs
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      open(unit=14,file='jmol.spt')
      open(unit=15,file='NTOs.html')
      write(15,'(a)')'<html>'
      write(15,'(a)')'<head>'
      write(15,'(a)')'<title>Natural transition orbitals</title>'
      write(15,'(a)')'</head>'
      write(15,'(a)')'<body>'
      write(15,'(a)')'</tr>'
      
      Do k=1, Nnto
      ! SVD of Xci(ia,k) => X(i,a) 
      ! X = U SIGMA V
      X=0.0
      Do i=1, nci
      X(iconf(i,1),iconf(i,2)-no)=Xci(i,k)
      enddo
      SIGMA=0.0
      U=0.0
      V=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',no,nv,X,no,SIGMA,U,no,V,nv,work,lwork,iwork,
     .infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',no,nv,X,no,SIGMA,U,no,V,nv,work,lwork,iwork,
     .infoo)      
      deallocate(work)
      ! Compute the participation ratio
      PR_nto=sum(SIGMA**2.0)**2.0/sum(SIGMA**4.0)
      if(k==1) write(*,*)'#state   PR_NTOs'
      write(*,*)k,PR_nto
      
      write(fname,'(a,i3.3,a)')'nto',k,'.molden'
      write(14,*)'load ',fname
      
      
      write(15,'(a,i0,a)')'<h2>NTO ',k,'</h2>'
      write(15,'(a)')'<table>'
      
      open(unit=13,file=fname)
      write(13,*) '[Molden Format]'
      write(13,*) '[Atoms] Angs'
      Do i=1, ncent
      write(13,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo      
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
      
      write(13,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(13,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(13,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(13,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(13,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(13,23)'f',info(i),'1.000000'
      
      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif
      
      
      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(13,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)      
      at=ipat(i)
      if(at/=ipat(i+1)) write(13,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo      
      counter=0
      write(13,*) '[MO]' 
      Do i=1,no
      if(SIGMA(i)**2.0>=0.05)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMA(i)**2.0,'</td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,no
      MO(j)=MO(j)+ca(j+(l-1)*nao)*U(l,i)
      enddo
      enddo
      write(13,*) 'Sym= Hole'
      write(13,*) 'Ene=',-SIGMA(i)**2.0
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMA(i)**2.0,'</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,nv
      MO(j)=MO(j)+ca(j+(l+no-1)*nao)*V(i,l)
      enddo
      enddo
      write(13,*) 'Sym= Particle'
      write(13,*) 'Ene=',-SIGMA(i)**2.0
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 2'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      write(15,*)'</table>'
      close(13)
      enddo
      write(15,'(a)')'</body>'
      write(15,'(a)')'</html>'
      write(*,*)'NTOs written'
 21   format(a,2i7,3f16.8)
 22   format(i0,3x,a)
 23   format(a,3x,i0,3x,a)
 24   format(2f16.8)
 25   format(3f16.8) 
 26   format(a,i5.5,a,i5.5,a)
 27   format(a,i5.5,a,i5.5,a,f5.3,a)
      end
      
      subroutine Uprint_nto_rpa(uci,hci,cca,mocia,nexa,nci,nroot,nao,
     .                   iconfa,
     .                   maxconfa,noa,nva,ccb,mocib,nexb,iconfb,
     .                   maxconfb,nob,nvb)
      implicit none 
      integer:: mocia,nci,nroot,nao,maxconfa,noa,nva
      integer:: nexa
      real*4:: uci(nci,nroot),hci(nci,nroot)
      real*4:: newuci(nci,nroot)
      real*8:: cca(nao*mocia)
      integer ::iconfa(maxconfa,2)
      
      integer:: mocib,maxconfb,nob,nvb
      integer:: nexb
      real*8:: ccb(nao*mocib)
      integer ::iconfb(maxconfb,2)
      
      newuci=uci+hci
      
      call Uprint_nto(uci,cca,mocia,nexa,nci,nroot,nao,iconfa,
     .                   maxconfa,noa,nva,ccb,mocib,nexb,iconfb,
     .                   maxconfb,nob,nvb)
      
      end
      
      
      subroutine Uprint_nto(Xci,cca,mocia,nexa,nci,nroot,nao,iconfa,
     .                   maxconfa,noa,nva,ccb,mocib,nexb,iconfb,
     .                   maxconfb,nob,nvb)
      use stdacommon
      use commonresp
      implicit none
      integer:: mocia,nci,nroot,nao,maxconfa,noa,nva
      real*4:: Xci(nci,nroot),Xa(noa,nva),MO(nao)
      real*8:: cca(nao*mocia)
      integer:: iconfa(maxconfa,2),nexa
      
      integer:: mocib,maxconfb,nob,nvb
      real*4:: Xb(nob,nvb)
      real*8:: ccb(nao*mocib)
      integer:: iconfb(maxconfb,2),nexb
      
      integer:: i,j,k,l
      
      real*4:: SIGMAa(noa),Ua(noa,noa),Va(nva,nva)
      real*4:: SIGMAb(nob),Ub(nob,nob),Vb(nvb,nvb)
      integer:: iworka(8*noa), infoo, lwork, iworkb(8*nob)
      real*4,allocatable :: work(:)
      
      real*4:: PR_ntoa, PR_ntob
      
      character(len=14):: fname
      
      integer:: ncent,nprims,nmo,nbf,imethod,icdim
      real*8,allocatable:: norm(:)
      
      integer:: prim,at
      integer,allocatable:: info(:),f_info(:)
      integer:: counter_BF
      integer:: flag,flag2,counter
      
      open(unit=11,file='NTOvar')
      read(11,*)ncent,nprims,nmo,icdim,nbf
      close(11,status='delete')      
      allocate(norm(nprims))
      open(unit=12,file='fnorm')
      ! read normalization
      Do i=1,nprims
      read(12,*)norm(i)
      enddo
      close(12,status='delete')
      open(unit=11,file='NTOao')
      allocate(ipat(nprims),ipty(nprims),ipao(nprims),info(nprims))
      allocate(atnam(ncent),exip(nprims),cxip(nprims),co(ncent,4))
      allocate(f_info(nbf))
      Do i=1, nprims
      read(11,121)ipat(i),ipty(i),ipao(i),exip(i),cxip(i)
      enddo
      close(11,status='delete')
      open(unit=11,file='NTOat')
      Do i=1,ncent
      read(11,*)atnam(i),co(i,1:4)
      enddo
      close(11,status='delete')

      
 121  format(3i10,3x,2f16.9)     
      
      
      
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in print NTOs
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      open(unit=14,file='jmol.spt')
      open(unit=15,file='NTOs.html')
      write(15,'(a)')'<html>'
      write(15,'(a)')'<head>'
      write(15,'(a)')'<title>Natural transition orbitals</title>'
      write(15,'(a)')'</head>'
      write(15,'(a)')'<body>'
      write(15,'(a)')'</tr>'
      
      Do k=1, Nnto
      
      ! alpha block
      
      ! SVD of Xci(ia,k) => X(i,a) 
      ! X = U SIGMA V
      Xa=0.0
      Do i=1, nexa
      Xa(iconfa(i,1),iconfa(i,2)-noa)=Xci(i,k)
      enddo
      SIGMAa=0.0
      Ua=0.0
      Va=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',noa,nva,Xa,noa,SIGMAa,Ua,noa,Va,nva,work,
     .lwork,iworka,infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',noa,nva,Xa,noa,SIGMAa,Ua,noa,Va,nva,work,
     .lwork,iworka,infoo)      
      deallocate(work)

      
      ! beta block
      
      Xb=0.0
      Do i=1, nexb
      Xb(iconfb(i,1),iconfb(i,2)-nob)=Xci(i+nexa,k)
      enddo
      SIGMAb=0.0
      Ub=0.0
      Vb=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',nob,nvb,Xb,nob,SIGMAb,Ub,nob,Vb,nvb,work,
     .lwork,iworkb,infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',nob,nvb,Xb,nob,SIGMAb,Ub,nob,Vb,nvb,work,
     .lwork,iworkb,infoo)      
      deallocate(work)
      ! Compute the participation ratio
      ! Compute the participation ratio
      PR_ntoa=(sum(SIGMAa**2.0)+sum(SIGMAb**2.0))**2.0
     ./(sum(SIGMAa**4.0)+sum(SIGMAb**4.0))
      if(k==1) write(*,*)'#state   PR_NTOs'
      write(*,*)k,PR_ntoa
      
      write(fname,'(a,i3.3,a)')'nto',k,'.molden'
      write(14,*)'load ',fname
      
      
      write(15,'(a,i0,a)')'<h2>NTO ',k,'</h2>'
      write(15,'(a)')'<table>'
      
      open(unit=13,file=fname)
      write(13,*) '[Molden Format]'
      write(13,*) '[Atoms] Angs'
      Do i=1, ncent
      write(13,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo      
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
      
      write(13,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(13,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(13,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(13,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(13,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(13,23)'f',info(i),'1.000000'
      
      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif
      
      
      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(13,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)      
      at=ipat(i)
      if(at/=ipat(i+1)) write(13,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo   
      
      counter=0
      write(13,*) '[MO]' 
      Do i=1,noa
      if(SIGMAa(i)**2.0>=0.05)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,' alpha </td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,noa
      MO(j)=MO(j)+cca(j+(l-1)*nao)*Ua(l,i)
      enddo
      enddo
      write(13,*) 'Sym= Hole',SIGMAa(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,'</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,nva
      MO(j)=MO(j)+cca(j+(l+noa-1)*nao)*Va(i,l)
      enddo
      enddo
      write(13,*) 'Sym= Particle',SIGMAa(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 1'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      
      
      
      ! beta
      

      Do i=1,noa
      if(SIGMAb(i)**2.0>=0.05)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAb(i)**2.0,' beta </td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,nob
      MO(j)=MO(j)+ccb(j+(l-1)*nao)*Ub(l,i)
      enddo
      enddo
      write(13,*) 'Sym= Hole',SIGMAb(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAb(i)**2.0,'</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,nvb
      MO(j)=MO(j)+ccb(j+(l+nob-1)*nao)*Vb(i,l)
      enddo
      enddo
      write(13,*) 'Sym= Particle',SIGMAb(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup= 1'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      
      write(15,*)'</table>'
      close(13)
      enddo
      write(15,'(a)')'</body>'
      write(15,'(a)')'</html>'
      write(*,*)'NTOs written'
 21   format(a,2i7,3f16.8)
 22   format(i0,3x,a)
 23   format(a,3x,i0,3x,a)
 24   format(2f16.8)
 25   format(3f16.8) 
 26   format(a,i5.5,a,i5.5,a)
 27   format(a,i5.5,a,i5.5,a,f5.3,a)
      end
      
      
      subroutine SFprint_nto(Xci,cca,mocia,nci,nroot,nao,
     .                   noa,nva,ccb,mocib,iconfb,
     .                   maxconfb,nob,nvb)
      use stdacommon
      use commonresp
      implicit none
      integer:: mocia,nci,nroot,nao,noa,nva
      real*4:: Xci(nci,nroot),Xa(noa,noa),MO(nao),X(noa,nvb)
      real*8:: cca(nao*mocia)
      
      integer:: mocib,maxconfb,nob,nvb
      real*4:: Xb(nvb,nvb)
      real*8:: ccb(nao*mocib)
      integer:: iconfb(maxconfb,2)
      
      integer:: i,j,k,l
      
      real*4:: SIGMAa(noa),Ua(noa,noa),Va(noa,noa)
      real*4:: SIGMAb(nvb),Ub(nvb,nvb),Vb(nvb,nvb)
      integer:: iworka(8*noa), infoo, lwork, iworkb(8*nvb)
      real*4,allocatable :: work(:)
      
      real*4:: PR_ntoa, PR_ntob, PR_nto
      
      character(len=14):: fname
      
      integer:: ncent,nprims,nmo,nbf,imethod,icdim
      real*8,allocatable:: norm(:)
      
      integer:: prim,at
      integer,allocatable:: info(:),f_info(:)
      integer:: counter_BF
      integer:: flag,flag2,counter
      
      integer::kmem,lmem,jmem,io,iv,ij,idum1,idum2,jmax,imax
      real*8::uu,pp,umax
      
      
      open(unit=11,file='NTOvar')
      read(11,*)ncent,nprims,nmo,icdim,nbf
      close(11,status='delete')      
      allocate(norm(nprims))
      open(unit=12,file='fnorm')
      ! read normalization
      Do i=1,nprims
      read(12,*)norm(i)
      enddo
      close(12,status='delete')
      open(unit=11,file='NTOao')
      allocate(ipat(nprims),ipty(nprims),ipao(nprims),info(nprims))
      allocate(atnam(ncent),exip(nprims),cxip(nprims),co(ncent,4))
      allocate(f_info(nbf))
      Do i=1, nprims
      read(11,121)ipat(i),ipty(i),ipao(i),exip(i),cxip(i)
      enddo
      close(11,status='delete')
      open(unit=11,file='NTOat')
      Do i=1,ncent
      read(11,*)atnam(i),co(i,1:4)
      enddo
      close(11,status='delete')

      
 121  format(3i10,3x,2f16.9)     
      
      
      
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in print NTOs
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      
      open(unit=14,file='jmol.spt')
      open(unit=15,file='NTOs.html')
      write(15,'(a)')'<html>'
      write(15,'(a)')'<head>'
      write(15,'(a)')'<title>Natural transition orbitals</title>'
      write(15,'(a)')'</head>'
      write(15,'(a)')'<body>'
      write(15,'(a)')'</tr>'
      
      ! Print sf determinants with largest coefficients
      
      open(unit=13,file='1st_sf_MO.molden')
      write(13,*) '[Molden Format]'
      write(13,*) '[Atoms] Angs'
      Do i=1, ncent
      write(13,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo      
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
      
      write(13,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(13,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(13,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(13,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(13,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(13,23)'f',info(i),'1.000000'
      
      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif
      
      
      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(13,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)      
      at=ipat(i)
      if(at/=ipat(i+1)) write(13,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo       
      
      umax=-1   
      kmem=1 
      ! first largest
      do j=1,nci
            io=iconfb(j,1)
            iv=iconfb(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
            uu=dble(Xci(j,1))
            pp=dble(Xci(j,1))
            if(dabs(uu*pp).gt.umax)then
               umax=uu*pp
               imax=io
               jmax=iv
               kmem=j
            endif
      enddo
      !write(*,*)imax,jmax
      
      ! print first alpha MO
      
      write(13,*) '[MO]'      
      write(13,*) 'Sym= ',imax
      write(13,*) 'Ene=',1
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup=',0
      i=imax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, cca((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, cca((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, cca((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, cca((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, cca((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, cca((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, cca((i-1)*nao+j)
      endif
      enddo
      
      write(13,*) 'Sym= ',jmax
      write(13,*) 'Ene=',2
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup=',1
      i=jmax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, ccb((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, ccb((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, ccb((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, ccb((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, ccb((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, ccb((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, ccb((i-1)*nao+j)
      endif
      enddo
      
      
      ! second largest
         umax=-1
         lmem=1
         do j=1,nci
            uu=dble(Xci(j,1))
            pp=dble(Xci(j,1))
            if((uu*pp).gt.umax.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfb(j,1)
               jmax=iconfb(j,2)
               lmem=j
            endif
         enddo
      !write(*,*)imax,jmax
      write(13,*) 'Sym= ',imax
      write(13,*) 'Ene=',3
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup=',0
      i=imax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, cca((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, cca((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, cca((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, cca((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, cca((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, cca((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, cca((i-1)*nao+j)
      endif
      enddo
      
      write(13,*) 'Sym= ',jmax
      write(13,*) 'Ene=',4
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup=',1
      i=jmax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, ccb((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, ccb((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, ccb((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, ccb((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, ccb((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, ccb((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, ccb((i-1)*nao+j)
      endif
      enddo
      
      ! third largest
         umax=-1
         jmem=-1
         do j=1,nci
            uu=dble(Xci(j,1))
            pp=dble(Xci(j,1))
            if((uu*pp).gt.umax.and.j.ne.kmem.and.j.ne.lmem)then
               umax=uu*pp
               imax=iconfb(j,1)
               jmax=iconfb(j,2)
               jmem=j
            endif
         enddo
      !write(*,*)imax,jmax
      write(13,*) 'Sym= ',imax
      write(13,*) 'Ene=',5
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup=',0
      i=imax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, cca((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, cca((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, cca((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, cca((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, cca((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, cca((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, cca((i-1)*nao+j)
      endif
      enddo
      
      write(13,*) 'Sym= ',jmax
      write(13,*) 'Ene=',6
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup=',1
      i=jmax
      Do j=1,nbf
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, ccb((i-1)*nao+j+2) ! 16
      if(f_info(j)==15)write(13,*) j, ccb((i-1)*nao+j-1) ! 14
      if(f_info(j)==16)write(13,*) j, ccb((i-1)*nao+j-1) ! 15
      if(f_info(j)==17)write(13,*) j, ccb((i-1)*nao+j+1) ! 18
      if(f_info(j)==18)write(13,*) j, ccb((i-1)*nao+j+1) ! 19
      if(f_info(j)==19)write(13,*) j, ccb((i-1)*nao+j-2) ! 17
      
      else
      write(13,*) j, ccb((i-1)*nao+j)
      endif
      enddo
      close(13)
      
      
      ! Spin-flip Xci(:,1)
      
      X=0.0
      Do i=1, nci
      X(iconfb(i,1),iconfb(i,2)-nob)=Xci(i,1)
      enddo
      SIGMAa=0.0
      Ua=0.0
      Va=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',noa,nvb,X,noa,SIGMAa,Ua,noa,Vb,nvb,work,lwork,
     .iworka,infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',noa,nvb,X,noa,SIGMAa,Ua,noa,Vb,nvb,work,lwork,
     .iworka,infoo)      
      deallocate(work)
      ! Compute the participation ratio
      PR_nto=sum(SIGMAa**2.0)**2.0/sum(SIGMAa**4.0)
      write(*,*)'PR_NTOs spin-flip 1st SF-state ',PR_nto
      
      write(fname,'(a,i3.3,a)')'nto',0,'.molden'
      write(14,*)'load ',fname
      
      
      write(15,'(a)')'<h2>NTO spin-flip 1st SF-state </h2>'
      write(15,'(a)')'<table>'
      
      open(unit=13,file=fname)
      write(13,*) '[Molden Format]'
      write(13,*) '[Atoms] Angs'
      Do i=1, ncent
      write(13,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo      
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
      
      write(13,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(13,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(13,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(13,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(13,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(13,23)'f',info(i),'1.000000'
      
      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif
      
      
      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(13,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)      
      at=ipat(i)
      if(at/=ipat(i+1)) write(13,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo      
      counter=0
      write(13,*) '[MO]' 
      Do i=1,noa
      if(SIGMAa(i)**2.0>=0.02)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',0,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',0,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,
     .' alpha </td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,noa
      MO(j)=MO(j)+cca(j+(l-1)*nao)*Ua(l,i)
      enddo
      enddo
      write(13,*) 'Sym= Hole'
      write(13,*) 'Ene=',-SIGMAa(i)**2.0
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',0,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',0,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,
     .' beta</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,nvb
      MO(j)=MO(j)+ccb(j+(l+nob-1)*nao)*Vb(i,l)
      enddo
      enddo
      write(13,*) 'Sym= Particle'
      write(13,*) 'Ene=',-SIGMAa(i)**2.0
      write(13,*) 'Spin= beta'
      write(13,*) 'Occup= 1'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      write(15,*)'</table>'
      close(13)

      ! SF states
      
      
      Do k=1, Nnto
      
      ! occ occ Xij=sum_a(Xia_1*Xja_(k+1))
      Xa=0.0
      Do i=1, nci
      Do j=1, nci
      if(iconfb(i,2)==iconfb(j,2))then
      Xa(iconfb(i,1),iconfb(j,1))=Xa(iconfb(i,1),iconfb(j,1))+
     .             Xci(i,1)*Xci(j,k+1)
      endif
      enddo
      enddo
      
      SIGMAa=0.0
      Ua=0.0
      Va=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',noa,noa,Xa,noa,SIGMAa,Ua,noa,Va,noa,work,
     .lwork,iworka,infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',noa,noa,Xa,noa,SIGMAa,Ua,noa,Va,noa,work,
     .lwork,iworka,infoo)      
      deallocate(work)
      
      
      ! unocc unocc Xab=sum_i(Xia_1*Xib_(k+1))
      Xb=0.0
      Do i=1, nci
      Do j=1, nci
      if(iconfb(i,1)==iconfb(j,1))then
      Xb(iconfb(i,2)-nob,iconfb(j,2)-nob)=
     .Xb(iconfb(i,2)-nob,iconfb(j,2)-nob)+
     .             Xci(i,1)*Xci(j,k+1)
      endif
      enddo
      enddo
      Xb=-Xb
      SIGMAb=0.0
      Ub=0.0
      Vb=0.0
      ! dimension of work
      allocate(work(1))
      lwork=-1
      call SGESDD('A',nvb,nvb,Xb,nvb,SIGMAb,Ub,nvb,Vb,nvb,work,
     .lwork,iworkb,infoo)
      lwork=work(1)
      deallocate(work)
      ! SVD
      allocate(work(lwork))
      call SGESDD('A',nvb,nvb,Xb,nvb,SIGMAb,Ub,nvb,Vb,nvb,work,
     .lwork,iworkb,infoo)      
      deallocate(work)
      ! Compute the participation ratio
      PR_ntoa=sum(SIGMAa**2.0)**2.0/sum(SIGMAa**4.0)
      if(k==1) write(*,*)'#state   PR_NTOs'
      PR_ntob=sum(SIGMAb**2.0)**2.0/sum(SIGMAb**4.0)
      PR_nto=(sum(SIGMAa**2.0)+sum(SIGMAb**2.0))**2.0
     ./(sum(SIGMAa**4.0)+sum(SIGMAb**4.0))
      write(*,*)k,PR_nto
      
      write(fname,'(a,i3.3,a)')'nto',k,'.molden'
      write(14,*)'load ',fname
      
      
      write(15,'(a,i0,a)')'<h2>NTO ',k,'</h2>'
      write(15,'(a)')'<table>'
      
      open(unit=13,file=fname)
      write(13,*) '[Molden Format]'
      write(13,*) '[Atoms] Angs'
      Do i=1, ncent
      write(13,21) atnam(i),i,int(co(i,4)),co(i,1:3)*0.52917721092
      enddo      
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
      
      write(13,*) '[GTO]'
      at=0
      counter_BF=0
      f_info=0
      flag=0
      flag2=0
      Do i=1,nprims
      if(at/=ipat(i)) write(13,22) ipat(i),'0'
      if(info(i)/=0)then
      counter_BF=counter_BF+1
      if(ipty(i)==1.and.exip(i)==exip(i+1).and.info(i)==1.and.
     .   ipty(i+1)==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.exip(i)==exip(i+info(i)).and.info(i)>1.and.
     .   ipty(i+info(i))==2)then
      flag=1
      flag2=i
      write(13,23)'sp',info(i),'1.000000'
      endif
      if(ipty(i)==1.and.flag==0)  write(13,23)'s',info(i),'1.000000'
      if(ipty(i)==2.and.flag==0)  write(13,23)'p',info(i),'1.000000'
      if(ipty(i)==5)  write(13,23)'d',info(i),'1.000000'
      if(ipty(i)==11) write(13,23)'f',info(i),'1.000000'
      
      if(ipty(i)==14.or.ipty(i)==15.or.ipty(i)==16.or.
     . ipty(i)==17.or.ipty(i)==18.or.ipty(i)==19) then
      f_info(counter_BF)=ipty(i)

      endif
      
      
      endif
      !
      !   Contractactions are normalized
      !
      if(ipty(i)==1.and.flag==1)then
      write(13,25)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i),
     .      cxip(i+info(flag2))*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i+info(flag2)))**2.0
     .      *dsqrt(2.0*exip(i+info(flag2)))))/norm(i+info(flag2))
      endif
      if(ipty(i)==1.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317
     .      /(2.0*exip(i)*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==2.and.flag==0) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.5
     .      /((2.0*exip(i))**2.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==5) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*0.75
     .      /((2.0*exip(i))**3.0*dsqrt(2.0*exip(i))))/norm(i)
      if(ipty(i)==11) write(13,24)exip(i),
     .      cxip(i)*dsqrt(5.5683279968317*1.875
     .      /((2.0*exip(i))**4.0*dsqrt(2.0*exip(i))))/norm(i)      
      at=ipat(i)
      if(at/=ipat(i+1)) write(13,*)
      if(flag==1.and.ipty(i)==4)then
      flag=0
      flag2=0
      endif
      enddo   
      
      counter=0
      write(13,*) '[MO]' 
      Do i=1,noa
      if(SIGMAa(i)**2.0>=0.05)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,
     .' occ-occ alpha </td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,noa
      MO(j)=MO(j)+cca(j+(l-1)*nao)*Va(i,l)
      enddo
      enddo
      write(13,'(a,f5.3)') 'Sym= Hole',SIGMAa(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAa(i)**2.0,'</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,noa
      MO(j)=MO(j)+cca(j+(l-1)*nao)*Ua(l,i)
      enddo
      enddo
      write(13,'(a,f5.3)') 'Sym= Particle',SIGMAa(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Alpha'
      write(13,*) 'Occup= 1'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      
      
      
      ! unocc-unocc
      

      Do i=1,nvb
      if(SIGMAb(i)**2.0>=0.05)then ! Print significative NTOs
      counter=counter+1
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,'(a)')'<tr>'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAb(i)**2.0,
     .' unocc-unocc beta</td>'
      ! hole MO
      MO=0.0
      Do j=1,nao
      Do l=1,nvb
      MO(j)=MO(j)+ccb(j+(l-1+nob)*nao)*Ub(l,i)
      enddo
      enddo
      write(13,'(a,f5.3)') 'Sym= Hole',SIGMAb(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup= 0'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo      
      counter=counter+1
      ! electron MO
      write(14,*)'mo ',counter
      write(14,*)'background white'
      write(14,*)'mo titleformat ""'
      write(14,*)'mo resolution 10'
      write(14,*)'mo nomesh'
      write(14,*)'mo fill'
      write(14,*)'mo cutoff 0.04'
      write(14,26)'write image png "NTO_',k,'_',counter,'.png"'
      write(15,27)'<td><img src="NTO_',k,'_',counter,
     .'.png" border="1" width="400"><br>',SIGMAb(i)**2.0,'</td>'
      write(15,'(a)')'<tr>'
      MO=0.0
      Do j=1,nao
      Do l=1,nvb
      MO(j)=MO(j)+ccb(j+(l+nob-1)*nao)*Vb(i,l)
      enddo
      enddo
      write(13,'(a,f5.3)') 'Sym= Particle',SIGMAb(i)**2.0
      write(13,*) 'Ene=',counter
      write(13,*) 'Spin= Beta'
      write(13,*) 'Occup= 1'
      Do j=1,nao
      
      ! f molden order
      if(f_info(j)==14.or.f_info(j)==15.or.f_info(j)==16.or.
     . f_info(j)==17.or.f_info(j)==18.or.f_info(j)==19) then
     
      if(f_info(j)==14)write(13,*) j, MO(j+2) ! 16
      if(f_info(j)==15)write(13,*) j, MO(j-1) ! 14
      if(f_info(j)==16)write(13,*) j, MO(j-1) ! 15
      if(f_info(j)==17)write(13,*) j, MO(j+1) ! 18
      if(f_info(j)==18)write(13,*) j, MO(j+1) ! 19
      if(f_info(j)==19)write(13,*) j, MO(j-2) ! 17
      
      else
      write(13,*) j, MO(j)
      endif
      enddo        
      
      endif
      enddo
      
      write(15,*)'</table>'
      close(13)
      enddo
      write(15,'(a)')'</body>'
      write(15,'(a)')'</html>'
      write(*,*)'NTOs written'
 21   format(a,2i7,3f16.8)
 22   format(i0,3x,a)
 23   format(a,3x,i0,3x,a)
 24   format(2f16.8)
 25   format(3f16.8) 
 26   format(a,i5.5,a,i5.5,a)
 27   format(a,i5.5,a,i5.5,a,f5.3,a)
      end      
      
      
      
      
      
      
      
      
      
      
