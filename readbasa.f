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
c reads the speical tm2molden binary file 

      subroutine readbas0a(mode,ncent,nmo,nbf,nprims,wfn) 
      use stdacommon
      implicit double precision (a-h,o-z)  

      character*(*)wfn
      character*80 out
      character*128 a128    
      character*20 a20      
      dimension xx(10)
      logical ex
      integer i,j,maxlen

      write(*,*)   
      write(*,*)'reading: ',wfn
      call header('M O / A O    I N P U T ',0)
      inquire(file=wfn,exist=ex)
      if(.not.ex)then
          write(*,*)'file:',wfn,' not found'
          stop
      endif
      open(unit=iwfn,file=wfn,form='unformatted')
      read(iwfn) nmo,nbf,nprims,ncent
      close(iwfn)

      ! determine length of ncent integer (for fitting printout with next routine to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      write(*,'(a)',advance='no')'atom '
      do i=1,maxlen-1
         write(*,'(a)',advance='no')' '
      enddo

      write(*,'(''#'',10x,''x'',13x,''y'',
     .                      13x,''z'',12x,''charge'')')

      end

      subroutine readbasa(mode,imethod,ncent,nmo,nbf,nprims,cc,
     .icdim,wfn,iaobas) 
      use stdacommon
      implicit double precision (a-h,o-z)  

      dimension cc(icdim)  
      integer imethod

      character*(*) wfn                               
      character*80 out                               
      character*128 a128    
      character*20 a20      
      logical ex,mosgen
      dimension xx(10)
      character*79 prntfrmt
      integer maxlen

      iaobas=0

      ! determine length of ncent integer (for printout to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      prntfrmt=' '
      write(prntfrmt,'(a,i0,a)')'(2x,a2,x,i',maxlen,
     .               ',2x,3f14.8,3x,f10.2)'

      iwfn=42
      open(unit=iwfn,file=wfn,form='unformatted')
      read(iwfn) nmo,nbf,nprims,ncent
      if(imethod.eq.2) nmo = 2*nmo
      do 100 i = 1,ncent                        
        read (iwfn) atnam(i),co(i,1),co(i,2),co(i,3),co(i,4)
        if(co(i,4).lt.1.0d0) atnam(i)='xx'
        write(*,prntfrmt) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
100   continue                                                 
      read(iwfn) (ipat(i),i=1,nprims)
c     ipat - primitive to atom                
      read(iwfn) (ipty(i),i=1,nprims)
c     ipty - angular momemtum type of primitive
      read(iwfn) (ipao(i),i=1,nprims)
c     ipao - primitive to contracted              
      read(iwfn) (exip(i),i=1,nprims) 
c     exip - exponents of primitives               
      read(iwfn) (cxip(i),i=1,nprims) 

! for debugging purposes
!      do i=1,nprims
!        write(*,*) i,ipty(i)
!        write(*,*) exip(i),cxip(i)
!          write(*,*)k,jprimao,jprimtyp,cxip(k),cxip(k)**2
!      enddo
          
      do i=1,nmo
      read(iwfn) occ(i),eps(i) 
!       write(*,*) occ(i),eps(i)       
      enddo
      do i=1,nmo
      read(iwfn) (cc(j+(i-1)*nbf),j=1,nbf)          
      enddo
!      do i=1,nmo
!       write(*,*) (cc(j+(i-1)*nbf),j=1,nbf)
!      enddo
      read(iwfn) tote,gamma                              
      close(iwfn)

      iaobas=idint(gamma)

      write(*,95)ncent,nmo,nprims,nbf
 95   format (/,1x,'# atoms          =',i5,/,
     .          1x,'# mos            =',i5,/,
     .          1x,'# primitive  aos =',i5,/,
     .          1x,'# contracted aos =',i5,/)   

      if(iaobas.eq.0)then
         write(*,*) 'spherical AO basis'
      else
         write(*,*) 'cartesian AO basis'
      endif

      call etafill(nprims)

!203   format(2x,a2,i3,2x,3f14.8,3x,f10.2)  

      end     
             
      subroutine readbasb(mode,imethod,ncent,nmo,nbf,nprims,cc,ccspin,
     .icdim,wfn,iaobas)
      use stdacommon
      implicit double precision (a-h,o-z)

      dimension cc(icdim)
      integer ccspin(nmo)
      integer imethod

      character*(*) wfn
      character*80 out
      character*128 a128
      character*20 a20
      logical ex,mosgen
      dimension xx(10)
      character*100 line
      integer iostatus
      character*5 spin,sym
      character*79 prntfrmt
      integer maxlen

      ! determine length of ncent integer (for printout to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      prntfrmt=' '
      write(prntfrmt,'(a,i0,a)')'(2x,a2,x,i',maxlen,
     .               ',2x,3f14.8,3x,f10.2)'


      iaobas=0

      iwfn=42
      open(unit=iwfn,file=wfn,form='unformatted')
      read(iwfn) nmo,nbf,nprims,ncent
      if(imethod.eq.2) nmo = 2*nmo
      do 100 i = 1,ncent
        read (iwfn) atnam(i),co(i,1),co(i,2),co(i,3),co(i,4)
        if(co(i,4).lt.1.0d0) atnam(i)='xx'
        write(*,prntfrmt) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
100   continue
      read(iwfn) (ipat(i),i=1,nprims)
c     ipat - primitive to atom                
      read(iwfn) (ipty(i),i=1,nprims)
c     ipty - angular momemtum type of primitive
      read(iwfn) (ipao(i),i=1,nprims)
c     ipao - primitive to contracted              
      read(iwfn) (exip(i),i=1,nprims)
c     exip - exponents of primitives               
      read(iwfn) (cxip(i),i=1,nprims)

      do i=1,nmo
      read(iwfn) occ(i),eps(i)
      enddo
      do i=1,nmo
      read(iwfn) (cc(j+(i-1)*nbf),j=1,nbf)
      enddo
      read(iwfn) tote,gamma
      close(iwfn)
   
      if(imethod.eq.2) then

      write(*,'(/,A,/)') 'Reading orbitals data from molden.input file '

      open(unit=iwfn,file='molden.input',status='OLD')
      do
       read(iwfn,'(A)',IOSTAT=iostatus) line
       if(line.eq.'[MO]'.or.iostatus.lt.0) exit
      enddo
      do i = 1, nmo
       read(iwfn,*) line, sym
       read(iwfn,*) line, eps(i)
       read(iwfn,*) line, spin
       if(spin.eq.'Alpha') then
        ccspin(i) = 1
       else
        ccspin(i) = 2
       endif
       read(iwfn,*) line, occ(i)
       do j = 1, nbf
        read(iwfn,*) ibf, ccmolden
       enddo
      enddo
      close(iwfn)

!      call header('Orbitals',0)
!      write(*,'(/,A,/)') '  Occupancy,  Energy (eV), Orbital Spin'
!      do i = 1, nmo
!       write(*,'(F8.2,F12.4,I4)') occ(i),eps(i)*27.21139,ccspin(i)
!      enddo

      endif

      iaobas=idint(gamma)

      write(*,95) ncent,nmo,nprims,nbf
 95   format (/,1x,'# atoms          =',i5,/,
     .          1x,'# mos            =',i5,/,
     .          1x,'# primitive  aos =',i5,/,
     .          1x,'# contracted aos =',i5,/)

      if(iaobas.eq.0)then
         write(*,*) 'spherical AO basis'
      else
         write(*,*) 'cartesian AO basis'
      endif

      call etafill(nprims)

!203   format(2x,a2,i3,2x,3f14.8,3x,f10.2)

      end 
