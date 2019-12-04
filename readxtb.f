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
ccccccccccccccccccccccccccccccccc
!    read out xTB     input     c
ccccccccccccccccccccccccccccccccc           
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! co(ncent,1:3) : Cartesian coordinates 
! co(ncent,4)   : nuclear charge
! cxip(nprims) : contraction coefficients of primitives
! exip(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients 
! eps(nmo)     : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO 
! ipat(ncent)  : # of atom, the primitive is located on


      subroutine readxtb0(imethod,ncent,nmo,nbf,nprims)
      implicit double precision (a-h,o-z)

      integer, intent( out ) :: imethod,ncent,nmo,nbf,nprims
      ! temporary variables
      integer ii,i,j,k,maxlen
      logical ex

      write(*,*)
      write(*,*)'reading: wfn.xtb'
      call header('M O / A O    I N P U T ',0)
      inquire(file='wfn.xtb',exist=ex)
      if(.not.ex)then
          write(*,*)'file: wfn.xtb not found'
          stop 'input file not found'
      endif
 
      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .      status='old')
      rewind(iwfn)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! read rhf/uhf flag                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      read(iwfn)imethod
! read dimensions
      read(iwfn)ncent,nbf,nmo,nprims
      close(29)

! determine length of ncent integer (for fitting printout with next routine to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      write(*,'(a)',advance='no')'atom '
      do i=1,maxlen-1
         write(*,'(a)',advance='no')' '
      enddo

      write(*,'(''#'',10x,''x'',13x,''y'',
     .                      13x,''z'',12x,''charge'')')

      return
      end


      subroutine readxtb(imethod,ncent,nmo,nbf,nprims,cc)
      use stdacommon
      implicit double precision (a-h,o-z)

      integer, intent( in ) :: imethod,ncent,nmo,nbf,nprims
      real*8, intent ( out ) :: cc(imethod*nbf*nmo)
      ! temporary variables
      integer ii,i,j,k,maxlen
      real*8 dum
      character*79 prntfrmt
! determine length of ncent integer (for printout to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      prntfrmt=' '
      write(prntfrmt,'(a,i0,a)')'(2x,a2,x,i',maxlen,
     .               ',2x,3f14.8,3x,f10.2)'

      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .      status='old')
      rewind(iwfn)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! read rhf/uhf flag                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      read(iwfn)ii
! read dimensions
      read(iwfn)ii,i,j,k
! now read coordinates
      do i = 1,ncent
        read(iwfn) atnam(i)
      enddo
      do i = 1,ncent
         do j=1,3
           read(iwfn) dum
           co(i,j)=dum
         enddo 
         read(iwfn) k
         co(i,4)=dble(k)
         if(co(i,4).lt.1.0d0) atnam(i)='xx'
      enddo
*************************
* print out coordinates *
*************************
      do i=1,ncent
        write(*,prntfrmt) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
      enddo
!303   format(2x,a2,i3,2x,3f14.8,3x,f10.2)

**************************
! Now read basis set data                       
**************************
! ipty
      do i=1,nprims
        read(iwfn) k
        ipty(i)=k 
      enddo
! ipat
      do i=1,nprims          
        read(iwfn) k 
        ipat(i) = k
      enddo
! ipao
      do i=1,nprims
        read(iwfn) k 
        ipao(i) = k
      enddo       

! first exponents, then contraction coefficients
      read(iwfn) exip(1:nprims)
      read(iwfn) cxip(1:nprims)
*********************
! now the mo data   *
*********************
      k=0
      if(imethod.eq.2) then 
        !uks case: nmo = nmo_a + nmo_b
        ! alpha first, beta second
        ! occs + energies
        k=nmo/2 
        read(iwfn) occ(1:k)     
        read(iwfn) eps(1:k)
        k=k+1
        read(iwfn) occ(k:nmo)
        read(iwfn) eps(k:nmo)
        ! read MO coefficients
        i=nmo*nbf/2
        read(iwfn) cc(1:i)
        i=i+1
        k=nmo*nbf
        read(iwfn) cc(i:k)       
      else
      !rks case
      ! occs + energies      
        read(iwfn) occ(1:nmo)
        read(iwfn) eps(1:nmo)
      ! read MO coefficients 
        read(iwfn) cc
      endif

      close(iwfn)

      write(*,95) ncent,nmo,nprims,nbf
 95   format (/,1x,'# atoms          =',i5,/,
     .          1x,'# mos            =',i5,/,
     .          1x,'# primitive  aos =',i5,/,
     .          1x,'# contracted aos =',i5,/)

      if(imethod*nbf.gt.nmo)then
         write(*,*) 'spherical AO basis'
      else
         write(*,*) 'cartesian AO basis'
      endif

      call etafill(nprims)
      
      return
      end


