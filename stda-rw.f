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

! adapted by Marc de Wegifosse 2018-2019

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ncent: number of atoms
C nmo  : number of MOs
C nao  : number of contracted AOs
C xyz  : array with atomic coordinates (1:3) and nuclear charge (4) in
C        Bohr
C c    : MO coefficients (nao*nmo)
C eps  : MO energies (au)
C occ  : occupation numbers
C iaoat: index array (1:nao) indicating on which atom the AO is centered
C thr  : energy threshold in eV up to which energy the excited states
C        are computed = spectral range (input in eV)
C thrp : threshold for perturbation selection of CSF (input)
C ax   : Fock exchange mixing parameter in DF used
C othr and vthr NOT used (computed)
C othr : occ. orbitals with up to <othr> lower than Fermi level are
C        included (input in eV)
C vthr : virt. orbitals with up to <vthr> above Fermi level are
C        included (input in eV)
C fthr : threshold for CSF consideration in PT2
C nvec : integer, # roots for which eigenvectors are wanted
!
!!! used logicals from commonlogicals !!!
C triplet: logical=.true. if triplet states are to be calculated
C rpachk : logical=.true. if sTD-DFT is performed
C eigvec : logical=.true. print eigenvectors, ggavec is needed if eigvec and GGAs are used together
C nvec : integer, # roots for which eigenvectors are wanted
C screen : prescreen in pt selection and for CSFs with small transition strengths
C dokshift : shift A(ia,ia) elements if K(ia,ia) is small
C
C writes file <tda.dat> for spectrum plotting
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE stda_rw(ncent,nmo,nao,xyz,c,eps,occ,iaoat,thr,thrp,
     .            ax,alphak,betaj,fthr,nvec)
      use commonlogicals
      use commonresp
      use omp_lib
      IMPLICIT NONE
c input:
      integer ncent,nmo,nao
      integer iaoat(*)
      real*8  thr,thrp,ax,othr,vthr
      real*8  c(*),eps(*),occ(*),xyz(4,*)
      logical ggavec,ex
c local varibles:

c l=dipole lengths, v=dipole velocity (d/dxyz), m=angular momentum
      real*8 ::dipole(3),dipole_mo(3)
      integer::iwrk,jwrk
      real*8, allocatable ::xl(:),yl(:),zl(:)
      real*8, allocatable ::xv(:),yv(:),zv(:)
      real*8, allocatable ::xm(:),ym(:),zm(:)
      real*8, allocatable ::help(:),scr(:),dum(:),x(:,:)
c MOs, orbital energies and CSF printout stuff
      real*8, allocatable ::ca(:),epsi(:)
      real*8, allocatable ::umerk(:,:),umrkx(:,:),umrky(:,:),umrkz(:,:)
      real*8, allocatable ::rvp(:)

c stuff for diag of TDA matrix (or rpa)
c critical regarding memory
      integer info,lwork,liwork,il,iu,nfound
      real*4, allocatable ::uci  (:,:)
      real*4, allocatable ::eci  (:)
      real*4, allocatable ::hci  (:,:)
      real*4, allocatable ::work (:)
c RPA stuff
      real*4, allocatable ::apb(:)
      real*4, allocatable ::ambsqr(:)
      real*4, allocatable ::amb(:)
c Linear response
      real*4 :: start_time, end_time, stda_time
      integer :: STATUS

ccccccccccccc
      real*4  vu,vl
      integer,allocatable ::iwork(:)
      integer,allocatable::isuppz(:)

c Löwdin MOs, repulsion terms, charges and half-transformed stuff
c critical regarding memory
      real*8, allocatable ::clow(:)
      real*4, allocatable ::gamj(:,:)
      real*4, allocatable ::gamk(:,:)
      real*4, allocatable ::qia(:,:),pia(:,:)
      real*4, allocatable ::pij(:,:),qab(:,:),qij(:,:)
      real*4, allocatable ::q1(:),q2(:),q3(:)
      real*4 sdot, integral

!c the maximum size of the TDA expansion space
      integer maxconf
      real*8, allocatable ::ed(:),edpt(:)
      integer, allocatable :: iconf(:,:),kconf(:,:)



! check for vector printout
      integer, allocatable ::  vecchk(:)

c intermediates
      real*8 omax,vmin,pert,de,ek,ej,ak,xc,rabx,ef,fthr,aksqrt
      real*8 beta1,alpha2,pp,hilf,uu,sss,rl,rv,time,coc(3)
      real*8 fl,fv,ec,p23,xp,umax,xvu,yvu,zvu,xmu,ymu,zmu,xlu,ylu,zlu
      real*8 xj,amat,xmolw,xk,betaj,alphak,beta2,alpha1,deps,loc,jii,jaa
      real*8 alp_real(6),sumf,xms,yms,zms
      integer moci,i,j,k,l,ii,ihomo,io,iv,ihilf,nex,new,jo,jv,m,jmem
      integer nci,jj,idum1,idum2,kmem,imax,jmax,lmem,ij,nroot,lin,ierr
      integer no,nv,n,nexpt,jhomo,nvec
      integer*8 imem1,imem2,imem3
c atomic Hubbard parameters
      real*8 eta(94)
c atomic masses
      common /amass/ ams(107)
      real*8 ams
      character*79 dummy

      call cpu_time(stda_time)

c just a printout
      call header('s T D A',0)

      thr =thr /27.211385050d0
c estimate the orbital energy window which corresponds to the desired
c spectra range thr
      deps=(1.+0.8*ax)*thr
c make it save
      deps=deps*2.0

      moci=0
      omax=-1d+42
      vmin= 1d+42

      do i=1,nmo
         if(occ(i).gt.1.990d0.and.eps(i).gt.omax) omax=eps(i)
         if(occ(i).lt.0.010d0.and.eps(i).lt.vmin)vmin=eps(i)
      enddo

! if eigenvectors are wanted in TM format,check now how many occupied there are in general
      if(eigvec) then
        jhomo=0
        do i=1,nmo
          if(occ(i).gt.1.990d0) jhomo=jhomo+1
        enddo
      endif

      othr=vmin-deps
      vthr=deps+omax

      write(*,*)'spectral range up to (eV)     : ', thr*27.211385050d0
      write(*,*)'occ MO cut-off (eV)           : ', othr*27.211385050d0
      write(*,*)'virtMO cut-off (eV)           : ', vthr*27.211385050d0
      write(*,*)'perturbation thr              : ', thrp
      if(fthr.lt.1.79d308) then
      write(*,*)'max. CSF selection range (eV) : ', fthr
      fthr = fthr /27.211385050d0
      endif
      write(*,*)'triplet                       : ', triplet

      do i=1,nmo
         if(occ(i).gt.1.990d0.and.eps(i).gt.othr)moci=moci+1
         if(occ(i).lt.0.010d0.and.eps(i).lt.vthr)moci=moci+1
      enddo

      allocate(
     .         xl(moci*(moci+1)/2),yl(moci*(moci+1)/2),
     .         zl(moci*(moci+1)/2),
     .         xv(moci*(moci+1)/2),yv(moci*(moci+1)/2),
     .         zv(moci*(moci+1)/2),
     .         xm(moci*(moci+1)/2),ym(moci*(moci+1)/2),
     .         zm(moci*(moci+1)/2),
     .         help(nao*(nao+1)/2),clow(nao*moci),
     .         scr(nao*nao),dum(nao*nao),x(nao,nao),
     .         gamj(ncent,ncent),
     .         gamk(ncent,ncent),
     .         ca(nao*moci),epsi(moci)
     .        )

      write(*,*)'MOs in TDA : ', moci

! make two cases: 1st one) eigenvectors are needed, 2) eigenvectors are not needed
      if(eigvec.or.nto) then ! we want eigenvectors to be printed out
      allocate(vecchk(nmo), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for vecchk'
        vecchk=0
        moci=0
        do i=1,nmo
           if(occ(i).gt.1.990d0.and.eps(i).gt.othr)then
              moci=moci+1
              do j=1,nao
                 ca(j+(moci-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsi(moci)=eps(i)
              vecchk(i)=moci
           endif
        enddo
        ihomo=moci
        do i=1,nmo
           if(occ(i).lt.0.010d0.and.eps(i).lt.vthr)then
              moci=moci+1
              do j=1,nao
                 ca(j+(moci-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsi(moci)=eps(i)
              vecchk(i)=moci
           endif
        enddo
      else ! no eigenvectors needed
        moci=0
        do i=1,nmo
           if(occ(i).gt.1.990d0.and.eps(i).gt.othr)then
              moci=moci+1
              do j=1,nao
                 ca(j+(moci-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsi(moci)=eps(i)
           endif
        enddo
        ihomo=moci
        do i=1,nmo
           if(occ(i).lt.0.010d0.and.eps(i).lt.vthr)then
              moci=moci+1
              do j=1,nao
                 ca(j+(moci-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsi(moci)=eps(i)
           endif
        enddo
      endif

      no=ihomo
      nv=moci-no
      write(*,*)'oMOs in TDA: ', no
      write(*,*)'vMOs in TDA: ', nv
      if(no.eq.0.or.nv.eq.0) then
        stop 'no CSF, increase energy threshold (-e option)'
      endif

      maxconf=no*nv
      allocate(ed(maxconf),edpt(maxconf), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for ed and edpt'
      allocate(iconf(maxconf,2), kconf(maxconf,2), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for iconf and kconf'
      ed=0.0d0
      edpt=0.0d0
      iconf=0
      kconf=0
c we arrange MOS according to energy from 1:HOMO to LUMO:MOCI
c (in TM they come in irreps)
      write(*,*)'sorting MOs ...'
c sort for E diag
      call sort_vec(moci,nao,c,epsi)

      write(*,*)'reading and transforming R..V..L AO ints ...'

c read L,V,M with xyz components each and transform
c to MO basis (original but sorted MOs in array ca)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            dipole lengths
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=31,file='xlint',form='unformatted',status='old')
      read(31) help
      call onetri(1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,xl)
      close(31,status='delete')
      open(unit=32,file='ylint',form='unformatted',status='old')
      read(32) help
      call onetri(1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,yl)
      close(32,status='delete')
      open(unit=33,file='zlint',form='unformatted',status='old')
      read(33) help
      call onetri(1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,zl)

      !dipole_mo(3)=0.0
      !Do i=1,moci
      !dipole_mo(3)=dipole_mo(3)-zl(lin(i,i))
      !enddo
      !write(*,*)dipole_mo(3)

      close(33,status='delete')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            magnetic dipole
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=34,file='xmint',form='unformatted',status='old')
      read(34) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,xm)
      close(34,status='delete')
      open(unit=35,file='ymint',form='unformatted',status='old')
      read(35) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,ym)
      close(35,status='delete')
      open(unit=36,file='zmint',form='unformatted',status='old')
      read(36) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,zm)
      close(36,status='delete')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            velocity dipole
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=37,file='xvint',form='unformatted',status='old')
      read(37) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,xv)
      close(37,status='delete')
      open(unit=38,file='yvint',form='unformatted',status='old')
      read(38) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,yv)
      close(38,status='delete')
      open(unit=39,file='zvint',form='unformatted',status='old')
      read(39) help
      call onetri(-1,help,dum,scr,ca,nao,moci)
      call shrink(moci,dum,zv)
      close(39,status='delete')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c            calc S^1/2 and q(GS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=40,file='sint',form='unformatted',status='old')
      read(40) help
      write(*,*) 'ints done.'
      close(40,status='delete')


      write(*,*) 'S^1/2 ...'
      call makel(nao,help,x)

      call dgemm('n','n',nao,moci,nao,1.d0,X,nao,CA,nao,0.d0,SCR,nao)
      write(*,*) 'S^1/2 orthogonalized MO coefficients done.'

c check and copy to real*4 array
      do i=1,moci
      sss=0.0d0
      do j=1,nao
         sss=sss+scr(j+(i-1)*nao)**2
         clow(j+(i-1)*nao)=scr(j+(i-1)*nao)
      enddo
      if(abs(sss-1.0d0).gt.1.d-2)then
         write(*,*) 'MO norm ',i,sss
         stop 'internal MO norm error'
      endif
      enddo

      deallocate(scr,dum,help,x)

!      call cpu_time(time)

c these are standard (CIS) factors for singlet
      ak=2.0d0
      if(triplet) ak=0.0d0

c     open(unit=1,file='~/.param')
c     read(1,*)beta1,beta2,alpha1,alpha2,fk
c     close(1)

c the global parameters of the method:
      beta1=0.20d0
      beta2=1.830d0
      alpha1=1.420d0
      alpha2=0.480d0

      if(betaj.lt.-99.0d0) then ! if no beta parameter was read in
      betaj=beta1+beta2*ax
      endif
      if(alphak.lt.-99.0d0) then ! if no alpha parameter was read in
      alphak=alpha1+alpha2*ax
      endif

      write(*,*)
      write(*,*) 'ax(DF)   : ',ax
      write(*,*) 's^K      : ',ak
      write(*,*) 'beta  (J): ',betaj
      write(*,*) 'alpha (K): ',alphak
      write(*,*)

c set gamma's
      call setrep(eta)
      write(*,*) 'hardness table read.'

      write(*,*) 'setting up gammas ...'
c distances for gamma calc
      xmolw=0
      do i=1,ncent
         ii=idint(xyz(4,i))
c ams is the atomic mass (-> mol weight for output file)
         xmolw=xmolw+ams(ii)
         do j=1,i
            jj=idint(xyz(4,j))
            xj  =0.50d0*(eta(ii)+eta(jj)) * ax
            xk  =0.50d0*(eta(ii)+eta(jj))
            rabx=sqrt((xyz(1,i)-xyz(1,j))**2
     .               +(xyz(2,i)-xyz(2,j))**2
     .               +(xyz(3,i)-xyz(3,j))**2)
            gamj(j,i)=1./(rabx**betaj+1./xj**betaj)**(1.0d0/betaj)
            gamk(j,i)=1./(rabx**alphak+1./xk**alphak)**(1.0d0/alphak)
            gamj(i,j)=gamj(j,i)
            gamk(i,j)=gamk(j,i)
         enddo
      enddo

      !
! write transition charges on disc
!
      write(*,*)'write transition charges on disc'
       open(unit=70,file='qii',form='unformatted',status='replace')
!       open(unit=71,file='qij',form='unformatted',status='replace')
      open(unit=710,file='pij',form='unformatted',status='replace')
      allocate(q1(ncent),q2(ncent),qij(ncent,no*(no+1)/2))
      q1=0.0
      q2=0.0
      Do i=1, no
      Do j=1, i-1
      call lo12pop(i,j,ncent,nao,iaoat,clow,q1)
!       write(71)q1
      qij(1:ncent,lin(i,j))=q1(1:ncent)
      enddo
      call lo12pop(i,i,ncent,nao,iaoat,clow,q1)
       write(70)q1
      qij(1:ncent,lin(i,i))=q1(1:ncent)
      q2(1:ncent)=q2(1:ncent)+q1(1:ncent)
      enddo
       close(70)
!       close(71)
      allocate(pij(ncent,no*(no+1)/2))
      call ssymm('l','l',ncent,no*(no+1)/2,1.0,gamj,ncent,qij,ncent,0.0
     .             ,pij,ncent)
      deallocate(qij)
      Do i=1, no*(no+1)/2
      write(710)pij(1:ncent,i)
      enddo
      deallocate(pij)
      close(710)
      open(unit=72,file='qaa',form='unformatted',status='replace')
      open(unit=73,file='qab',form='unformatted',status='replace')
      Do i=no+1, moci
      Do j=no+1, i-1
      call lo12pop(i,j,ncent,nao,iaoat,clow,q1)
      write(73)q1
      enddo
      call lo12pop(i,i,ncent,nao,iaoat,clow,q1)
      write(72)q1
      enddo
      close(72)
      close(73)
      open(unit=74,file='qia',form='unformatted',status='replace')
      open(unit=740,file='pia',form='unformatted',status='replace')
      allocate(qia(ncent,no*nv))
      Do i=1, no
      Do j=no+1, moci
      call lo12pop(i,j,ncent,nao,iaoat,clow,q1)
      write(74)q1
      ij=(i-1)*nv+j-no
      qia(1:ncent,ij)=q1(1:ncent)
      enddo
      enddo
      close(74)
      allocate(pia(ncent,no*nv))
      pia=0.0
      call ssymm('l','l',ncent,no*nv,1.0,gamk,ncent,qia,ncent,0.0
     .             ,pia,ncent)
      Do i=1, no*nv
      write(740)pia(1:ncent,i)
      enddo
      close(740)
      deallocate(clow)

      write(*,'(/'' SCF atom population (using active MOs):'')')
      write(*,'(10F7.3)')q2(1:ncent)*2.0
      write(*,*)
      write(*,'('' # electrons in TDA:'',F8.3)') 2.0*sum(q2(1:ncent))
      write(*,*)

      open(unit=70,file='qii',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      allocate(qij(ncent,no),qab(ncent,nv))

      Do i=1, no
      read(70)qij(1:ncent,i)
      enddo
      Do i=1, nv
      read(72)qab(1:ncent,i)
      enddo
      close(70)
      close(72)

      allocate(pij(ncent,no), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii| intermediate'
      pij=0.0
      call ssymm('l','l',ncent,no,1.0,gamj,ncent,qij,ncent,0.0,pij
     .           ,ncent)
      allocate(uci(nv,no), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii|aa) matrix'
      uci=0.0
      ! now calc (ii|aa)^J
      call sgemm('t','n',nv,no,ncent,1.0,qab,ncent,pij,ncent,0.0,
     .           uci,nv)
      deallocate(pij)

      allocate(q3(1:ncent))
c determine singles which are lower than thr
      k=0
      j=0
      do io=1,no ! occ loop
         q1(1:ncent)=qij(1:ncent,io) !qii
         do iv=no+1,moci ! virt loop
            de=epsi(iv)-epsi(io)
            q2(1:ncent)=qab(1:ncent,iv-no) !qaa
            ej=dble(uci(iv-no,io))
            de=de-ej
            i=nv*(io-1)+iv-no
            q2(1:ncent)=pia(1:ncent,i)
            q3(1:ncent)=qia(1:ncent,i) !qia
            ek=sdot(ncent,q2,1,q3,1)
            de=de+ak*ek

! optional: perform K(ia,ia) dependent shift
            if(dokshift) then
              call kshift_to_ediag(de,ek)
            endif

c the primary ones
            if(de.le.thr)then
               k=k+1
               iconf(k,1)=io
               iconf(k,2)=iv
               ed(k)=de
            endif
c for PT
            if(de.gt.thr.and.de.lt.fthr)then ! needs to be on if fthr is specified
               j=j+1
               kconf(j,1)=io
               kconf(j,2)=iv
               edpt(j)=de
            endif

         enddo ! virt loop
      enddo ! occ loop
      deallocate(qij,qab,qia,pia,uci)

      nci=nex
      nex=k
      nexpt=j

      write(*,*)
      write(*,*)nex,'CSF included by energy.'
      write(*,*)
      write(*,*)nexpt,'considered in PT2.'

c errors and warning
      if(nex.lt.1) stop 'no CSF, increase energy threshold (-e option)'
!      if(nex.eq.maxconf)
!     .   stop 'primary CSF space exceeded. use -e option!'
!      if(nexpt.eq.maxconf)
!     .   write(*,*)'CSF PT2 space exceeded. try -p option!'

c sort for E diag
      do 141 ii = 2,nex
         i = ii - 1
         k = i
         pp= ed(i)
         do 121 j = ii, nex
            if (ed(j) .ge. pp) go to 121
            k = j
            pp=ed(j)
  121    continue
         if (k .eq. i) go to 141
         ed(k) = ed(i)
         ed(i) = pp
         do m=1,2
         ihilf=iconf(i,m)
         iconf(i,m)=iconf(k,m)
         iconf(k,m)=ihilf
         enddo
  141 continue

c just printout
      write(*,*)'ordered frontier orbitals'
      write(*,*)'        eV'
      j=max(1,no-10)
      do i=j,no
         write(*,'(i4,F10.3,F8.1)') i,epsi(i)*27.21139
      enddo
      write(*,*)
      j=min(moci,no+11)
      do i=no+1,j
         write(*,'(i4,F10.3,F8.1)') i,epsi(i)*27.21139
      enddo

      open(unit=70,file='qii',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      open(unit=74,file='qia',form='unformatted',status='old')
      allocate(qij(ncent,no),qab(ncent,nv),qia(ncent,no*nv))

      Do i=1, no
      read(70)qij(1:ncent,i)
      enddo
      Do i=1, nv
      read(72)qab(1:ncent,i)
      enddo
      Do i=1, no
      Do j=no+1, moci
      ij=(i-1)*nv+j-no
      read(74)qia(1:ncent,ij)
      enddo
      enddo
      close(70,status='delete')
      close(72)
      close(74)

      write(*,*)
      write(*,*)'            lowest CSF states'
      write(*,*)'      eV     nm      excitation i->a               eV'
      do i=1,min(nex,25)
         io=iconf(i,1)
         iv=iconf(i,2)
         q1(1:ncent)=qij(1:ncent,io)
         jii=integral(q1,q1,gamj,ncent)
         k=iv-no
         q2(1:ncent)=qab(1:ncent,k)
         ej=integral(q1,q2,gamj,ncent)
         jaa=integral(q2,q2,gamj,ncent)
         l=nv*(io-1)+k
         q1(1:ncent)=qia(1:ncent,l)
         ek=integral(q1,q1,gamk,ncent)
! de is now the Kia shift
         de=0
         loc=ej/sqrt(jii*jaa) ! locality
         if(dokshift) call kshift_to_ediag(de,ek)
         write(*,14) i,27.211*ed(i),
     .   1.d+7/(ed(i)*2.19474625d+5),iconf(i,1:2),
     .   27.211*(epsi(iv)-epsi(io)),27.211*ej,27.211*ek,27.211*de,loc
      enddo

      deallocate(q1,q2,q3,qij,qab,qia)
 14   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3,
     .       3x,'Kshft:',f8.3,2x,'locality:',f6.3,E12.5)
!       call cpu_time(hilf)
!       write(*,*) 'time elapsed:',hilf-time

      write(*,*)
      write(*,*)'selecting CSF ...'
      if(pt_off)then
      deallocate(edpt,kconf)
      nroot=nex
      nci=nex
      write(*,*)nci,'CSF in total.'
      else
      call ptselect_rw(nex,ncent,no,nv,nexpt,maxconf,iconf,kconf,
     .              ak,ax,ed,edpt,gamj,gamk,thrp,new,moci)
      deallocate(edpt,kconf)

c nroot at this point is the number of primary CSF. The
c number of roots in the range 0-thr (as desired) is not
c known but will be determined by the diag routine.
      nroot=nex
      write(*,*)new,'CSF included by PT.'
      nci=nex+new
      write(*,*)nci,'CSF in total.'
      endif
!      call cpu_time(hilf)
!      write(*,*) 'time elapsed:',hilf-time
      if(rpachk) then
**************************
c Obtional RPA procedure *
**************************
      write(*,*) 'sTD-DFT procedure...'
      write(*,*) 'setting up A+B and A-B matrices'

c allocate A+B and A-B in packed form
      allocate( apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2),
     .          stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A+B or A-B'

      call rrpamat_rw(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed,gamj
     .             ,gamk,apb,ambsqr,moci)

*****************************
c Linear Response functions *
*****************************
      if(optrota)then
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53,status='delete')
      if(velo_OR==.false.)call optrot(nci,apb,amb,iconf,maxconf,
     .xl,yl,zl,moci,no,nv,xm,ym,zm,xmolw)
      if(velo_OR==.true.)call optrot_velo(nci,apb,amb,iconf,maxconf,
     .xv,yv,zv,moci,no,nv,xm,ym,zm,xmolw)
      call cpu_time(end_time)
      print '("Opt. Rot.   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      CALL EXIT(STATUS)
      endif
      if(resp) then
      if(triplet) stop 'not available'
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53,status='delete')

      call lresp1(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
     .                     no,nv)

      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

!       call lresp1_noinv(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv)
!
!       call cpu_time(end_time)
!       print '("Lresp   Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0
      print '("sTD-DFT-rw Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      CALL EXIT(STATUS)
      endif

      if(aresp) then
      if(triplet) stop 'not available'
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53,status='delete')

      call lresp(nci,apb,ambsqr,iconf,maxconf,xl,yl,zl,moci,
     .                     no,nv)

      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

!       call lresp_noinv(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv)
!       call lresp_noinv1(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv)
!
!       call cpu_time(end_time)
!       print '("Lresp   Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0
      print '("sTD-DFT-rw Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      CALL EXIT(STATUS)
      endif



c big arrays not needed anymore
      deallocate(ed)
!      call prmat4(6,ambsqr,nci,0,'A-B^0.5')

************************************************************************************
! if eigenvectors in TM format wanted for GGAs, this is done here (not so nice)
************************************************************************************
      ggavec=.false.
      if (ax.eq.0.0d0.and.eigvec) then
        open(unit=39,file='TmPvEcInFo',status='replace')
        ! print to temporary file
        if(triplet) then
          write(39,*) 1
        else
         write(39,*) 0
        endif
        write(39,*) nvec,nmo,jhomo
        do i=1,nmo
          write(39,*) vecchk(i)
        enddo
        do i=1,nci
           write(39,*) iconf(i,1),iconf(i,2)
        enddo
        close(39)
        ggavec=.true.
        eigvec=.false.
      endif
************************************************************************************

      allocate( eci(nci) ,stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for eigenvalue vector'
      allocate(hci(nci,nci),uci(nci,nci),stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for eigenvector matrix'


c call sRPA routine (solve RPA problem)
c uci is now X, hci is Y
      call srpapack(nci,thr,ambsqr,apb,eci,uci,hci,info,ggavec)
      if(info.lt.1) stop 'internal error in diag'
      nroot=info
      info=0


!       if(resp)then
!
!       call pol_sos(nroot,nci,eci,uci,hci,xl,yl,zl,moci,
!      .maxconf,iconf,ak)
!
!
!       endif

*****************************
c Lin. Response func. 2PA   *
*****************************
      if(TPA)then
      if(triplet) stop 'not available'
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53)

      call lresp_2PA(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
     .                     no,nv,eci,uci,hci,nroot)

      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

!       call lresp_2PA_noinv(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv,eci,uci,hci,nroot)
!
!       call cpu_time(end_time)
!       print '("Lresp   Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0

      print '("sTD-DFT-rw Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      write(*,*)
      !CALL EXIT([STATUS])
      endif

      deallocate(apb,ambsqr)
!********************************************************************************
      open(unit=53,file='amb',form='unformatted',status='old')
      close(53,status='delete')
      else

!********************************************************************************
cccccccccccccccccccccccccc
c standard TDA procedure c
cccccccccccccccccccccccccc

!********************************************************************************
! construct ( 0.5 * B ) for X trafo (velocity correction) and print to file *
!********************************************************************************
      if(velcorr) then
        call rtdacorr_rw(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed
     .                ,gamj,gamk,moci)
      endif
!********************************************************************************

      allocate( hci(nci,nci), stat=ierr  )
      if(ierr.ne.0)stop 'allocation failed for TDA matrix'
      write(*,*)'calculating TDA matrix ...'
      call rtdamat_rw(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed,gamj
     .             ,gamk,hci,moci)
      deallocate(gamj,gamk)
!      call prmat4(6,hci,nci,nci,'A-Matrix')
!********************************************************************************

c big arrays not needed anymore
      deallocate(ed)

      write(*,*)'diagonalizing ...'
      write(*,'('' estimated time (min) '',f8.2)')
     .            float(nci)**2*float(nroot)/8.d+8/60.

c if LAPACK does not work
c     allocate(eci(nci),uci(nci,nroot),
c    .         stat=ierr)
c     call sHQRII(hci,nci,nroot,eci,uci)

c faster by a factor of 2-3
      lwork =26*nci
      liwork=10*nci
c we allocate uci with nci (and not with nroot) as save
c choice (other values gave segfaults)
      nroot=min(nci,int(1.5*nroot))
      allocate(eci(nci),uci(nci,nci),work(lwork),
     .         iwork(liwork),isuppz(nci))
      if(ierr.ne.0)stop 'allocation failed for TDA matrix diag'

      vl=0
      vu=thr
      call ssyevr('V','V','U',nci,hci,nci,vl,vu,il,iu,1.e-6,
     .            nfound,eci,uci,nci,isuppz,
     .            work,lwork,iwork,liwork,info)

      nroot=nfound
      if(nfound.lt.1) stop 'internal error in diag'

c     call prmat4(6,uci,nci,nci,'uci')

!!      internal check for orthonormality
!         hci=0.0
!         call sgemm('T','n',nci,nroot,nci,1.0,uci,nci,uci,nci,0.0,
!     .              hci,nci)
!         do i=1,min(12,nroot)
!          write(*,'(12f10.6)') (hci(j,i),j=1,min(12,nroot))
!         enddo
!        write(*,*)
!        do i=max(1,nroot-11),nroot
!          write(*,'(12f10.6)') (hci(j,i),j=max(1,nroot-11),nroot)
!        enddo
!!!!

      endif

      write(*,'(i5,'' roots found, lowest/highest eigenvalue : '',
     .2F8.3,i4)')nroot,eci(1)*27.21139,eci(nroot)*27.21139,info
      if(info.gt.0) stop 'diag error (ssyevr)'

      allocate(umerk(14,nroot))
c contract WF with MO integrals for intensities
      p23= 2.0d0/3.0d0
      aksqrt=sqrt(ak)

! if aniso, store and print the x,y,z-resolved data
      if(.not.velcorr.and.aniso) then
        allocate(umrkx(4,nroot),umrky(4,nroot),umrkz(4,nroot),stat=ierr)
        if(ierr.ne.0) stop 'error in xyz-umerk alloc'
        umrkx=0.0d0
        umrky=0.0d0
        umrkz=0.0d0
      endif

! if desired, print out data for exciton coupling, i.e., transition moments
      if(printexciton) then
        open(unit=27,file='tda.exc')
        write(27,*)ncent,nroot
        write(27,'(a)',advance='yes')'$coord'
        ihilf=0
        dummy=' '
        write(dummy,'(a)')'(x,f14.8,x,f14.8,x,f14.8,x,i0)'
        call cofc(ncent,xyz,coc)
        do i=1,ncent
          write(27,dummy) xyz(1:3,i)-coc(1:3),nint(xyz(4,i))
        enddo
        write(27,'(a)',advance='yes')'$states'
      endif
      allocate(q1(ncent))
      alp_real=0.0d0
      sumf=0.0d0
      allocate(rvp(nroot))
      rvp=0.0d0
c loop over excited states
      do i=1,nroot
         de=dble(eci(i))
         ef=1.0d0/(de+1.0d-8)
         xlu=0.0d0
         ylu=0.0d0
         zlu=0.0d0
         xmu=0.0d0
         ymu=0.0d0
         zmu=0.0d0
         xvu=0.0d0
         yvu=0.0d0
         zvu=0.0d0
         xms=0.0d0
         yms=0.0d0
         zms=0.0d0
         l=0
         umax=-1
         kmem=1
         q1=0.0
c loop over  CSF
         do j=1,nci
            io=iconf(j,1)
            iv=iconf(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if(dabs(uu*pp).gt.umax)then
               umax=uu*pp
               imax=io
               jmax=iv
               kmem=j
            endif
            xlu=xlu+xl(ij)*uu
            ylu=ylu+yl(ij)*uu
            zlu=zlu+zl(ij)*uu
            xvu=xvu+xv(ij)*pp
            yvu=yvu+yv(ij)*pp
            zvu=zvu+zv(ij)*pp
            xmu=xmu+xm(ij)*pp
            ymu=ymu+ym(ij)*pp
            zmu=zmu+zm(ij)*pp
            if(.not.velcorr.and.printexciton) then
               l=(io-1)*nv+(iv-no)
               q1(1:ncent)=q1(1:ncent)+qia(1:ncent,l)*real(uu)
            endif
         enddo

         ! multiply with factor from spin-integration
         xlu=xlu*aksqrt
         ylu=ylu*aksqrt
         zlu=zlu*aksqrt
         xvu=xvu*aksqrt
         yvu=yvu*aksqrt
         zvu=zvu*aksqrt
         xmu=xmu*aksqrt
         ymu=ymu*aksqrt
         zmu=zmu*aksqrt

c polarizability
         alp_real(1)=alp_real(1)+xlu*xlu*ef
         alp_real(2)=alp_real(2)+xlu*ylu*ef
         alp_real(3)=alp_real(3)+ylu*ylu*ef
         alp_real(4)=alp_real(4)+xlu*zlu*ef
         alp_real(5)=alp_real(5)+ylu*zlu*ef
         alp_real(6)=alp_real(6)+zlu*zlu*ef

         if(.not.velcorr.and.printexciton) then
           xms=(yvu*coc(3)-zvu*coc(2))
           yms=(zvu*coc(1)-xvu*coc(3))
           zms=(xvu*coc(2)-yvu*coc(1))
           write(27,'(i5,a,x,F10.4,x,a2)') i,':',de*27.21139,'eV'
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'l:',xlu,ylu,zlu
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'p:',xvu,yvu,zvu
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'m:',xmu-xms,ymu-yms,
     .                                                  zmu-zms
           do k=1,ncent
             write(27,'(f12.8)') real(aksqrt)*q1(k) ! factor of 2**0.5 arises from spin integration
           enddo
!           write(27,*)
         endif

         xp=xlu*xlu+ylu*ylu+zlu*zlu
         fl=de*p23*xp
         xp=xvu*xvu+yvu*yvu+zvu*zvu
         fv=ef*p23*xp
         xp=xlu*xmu+ylu*ymu+zlu*zmu
         rl=-235.7220d0*xp
! multiplying with ak (2=singlet excitation -> restricted case) is already included in moments: see book by Harada & Nakanishi
! recalculated value for Rot with recent CODATA values and changed it (used to be 235.730)
! the factor is the conversion from e*a0*mu_b in atomic units to 10^40 cgs (m_b is the Bohr magneton and a0 the Bohr radius)
         xp=xvu*xmu+yvu*ymu+zvu*zmu
         rv=-235.7220d0*xp*ef
         rvp(i)=rv
c sum rule
         sumf=sumf+fl

         umerk(1,i)=de
         umerk(2,i)=fl
         umerk(3,i)=fv
         umerk(4,i)=rl
         umerk(5,i)=rv
         umerk(6,i)=dble(uci(kmem,i))
         if(rpachk) umerk(6,i)=dble(uci(kmem,i))**2-dble(hci(kmem,i))**2
         umerk(7,i)=imax
         umerk(8,i)=jmax
! now find 2nd largest contribution
         umax=-1
         lmem=1
         do j=1,nci
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconf(j,1)
               jmax=iconf(j,2)
               lmem=j
            endif
         enddo
         umerk(9,i) =dble(uci(lmem,i))
         if(rpachk) umerk(9,i)=dble(uci(lmem,i))**2-dble(hci(lmem,i))**2
         umerk(10,i)=imax
         umerk(11,i)=jmax
         umax=-1
! now find 3rd largest contribution
         do j=1,nci
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.kmem.and.j.ne.lmem)then
               umax=uu*pp
               imax=iconf(j,1)
               jmax=iconf(j,2)
               jmem=j
            endif
         enddo
         umerk(12,i)=dble(uci(jmem,i))
         if(rpachk)umerk(12,i)=dble(uci(jmem,i))**2-dble(hci(jmem,i))**2
         umerk(13,i)=imax
         umerk(14,i)=jmax

         if(.not.velcorr.and.aniso) then
           sss=de*2.0d0
           umrkx(1,i)=xlu*xlu*sss
           umrky(1,i)=ylu*ylu*sss
           umrkz(1,i)=zlu*zlu*sss
           sss=ef*2.0d0
           umrkx(2,i)=xvu*xvu*sss
           umrky(2,i)=yvu*yvu*sss
           umrkz(2,i)=zvu*zvu*sss
           sss=-235.7220d0
           umrkx(3,i)=xlu*xmu*sss
           umrky(3,i)=ylu*ymu*sss
           umrkz(3,i)=zlu*zmu*sss
           sss=-235.7220d0*ef
           umrkx(4,i)=xvu*xmu*sss
           umrky(4,i)=yvu*ymu*sss
           umrkz(4,i)=zvu*zmu*sss
         endif
      enddo


c sort (output)
c     do 140 ii = 2,nroot
c        i = ii - 1
c        k = i
c        pp= umerk(1,i)
c        do 120 j = ii, nroot
c           if (umerk(1,j) .ge. pp) go to 120
c           k = j
c           pp=umerk(1,j)
c 120    continue
c        if (k .eq. i) go to 140
c        umerk(1,k) = umerk(1,i)
c        umerk(1,i) = pp
c        do m=2,11
c           hilf=umerk(m,i)
c           umerk(m,i)=umerk(m,k)
c           umerk(m,k)=hilf
c        enddo
c 140 continue

 11   format(i5,f9.3,f8.1,2f11.4,3x,3(F6.2,'(',i4,'->',i4,')'))
 12   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3,
     .       3x,'Kshft:',f8.3)
 13   format(                   41x,f8.3,13x,f8.3)

      if(.not.rpachk) deallocate(hci)


      if(.not.velcorr) then

        write(*,*)'writing spectral data to tda.dat ...'
        call print_tdadat(nroot,xmolw,eci,umerk(2:5,:),thr,'tda.dat')

        if(aniso) then
          write(*,*)'writing anisotropic data ...'
          write(*,*)'...tdax.dat'
          call print_tdadat(nroot,xmolw,eci,umrkx(1:4,:),thr,'tdax.dat')
          write(*,*)'...tday.dat'
          call print_tdadat(nroot,xmolw,eci,umrky(1:4,:),thr,'tday.dat')
          write(*,*)'...tdaz.dat'
          call print_tdadat(nroot,xmolw,eci,umrkz(1:4,:),thr,'tdaz.dat')
          deallocate(umrkx,umrky,umrkz)
          write(*,*)'data given such that f = (f_x + f_y + f_z)/3'
          write(*,*)'                 and R = R_x + R_y + R_z'
        endif
      else
         ! if velocity correction is on, get improved rotatory strengths
         if (printexciton) then
           call apbtrafoexc(nci,nroot,uci,eci,xl,yl,zl,xv,yv,zv,xm,ym,zm
     .                  ,xmolw,no,nv,coc,ncent,qia,maxconf,iconf,ak,rvp)
         else
           call apbtrafo(nci,nroot,uci,eci,xl,yl,zl,xv,yv,zv,xm,ym,zm
     .                   ,xmolw,maxconf,iconf,ak,rvp)
         endif
         do i=1,nroot
            umerk(5,i)=rvp(i)
         enddo
      endif

      if(printexciton) then
        deallocate(qia)
        close(27)
      endif


c output
      write(*,*)
      write(*,*)
     .'excitation energies, transition moments and TDA amplitudes'
      if(velcorr) then
        write(*,*) 'state    eV      nm       fL        Rv(corr)'
      else
        write(*,*) 'state    eV      nm       fL         Rv'
      endif
      do i=1,nroot
         ec=umerk(1,i)
         write(*,11)
     .   i,ec*27.21139,
     .   1.d+7/(ec*2.19474625d+5),umerk(2,i),umerk(5,i)
     .   ,umerk(6,i),idint(umerk(7,i)),idint(umerk(8,i))
     .   ,umerk(9,i),idint(umerk(10,i)),idint(umerk(11,i))
     .   ,umerk(12,i),idint(umerk(13,i)),idint(umerk(14,i))
      enddo

      alp_real=alp_real*2.0d0
      call prmat(6,alp_real,3,0,'alpha tensor')
      write(*,*) 'trace alpha_L[0] / au:',
     .(alp_real(1)+alp_real(3)+alp_real(6))/3.d0
      write(*,*) 'sum rule f_L        ',sumf

      call sosor(nroot,xmolw,eci,rvp)
      deallocate(rvp,q1)

*****************************
c Lin. Response func. ESA   *
*****************************
      if(ESA)then
      if(triplet) stop 'not available'
      ! Unrelaxed state-to-state transition dipole moments
      call cpu_time(start_time)
      if(rpachk) then
      call lresp_ESA(nci,iconf,maxconf,xl,yl,zl,moci,
     .              no,nv,eci,uci,hci,nroot,xmolw,thr)
      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      else
      call lresp_ESA_tda(nci,iconf,maxconf,xl,yl,zl,moci,
     .              no,nv,eci,uci,nroot,xmolw,thr)
      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTDA Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      endif
      write(*,*)
      !CALL EXIT([STATUS])
      endif

!       call hyperpol_sos(nroot,nci,eci,uci,hci,xl,yl,zl,moci,
!      .maxconf,iconf,no,nv) !to test the moments
!       call tpa_sos(nroot,nci,eci,uci,hci,xl,yl,zl,moci,
!      .maxconf,iconf,no,nv) !idem, it gives also information about how lresp relaxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optional: print eigenvectors in TM format
! Write $eigenpairs to scratch file
      if (eigvec) then
        ihilf=0 ! use as switch for output
        if(triplet) ihilf=1

        open(unit=39,file='TmPvEcInFo',status='replace')
        ! print to temporary file
        write(39,*) ihilf
        write(39,*) nvec,nmo,jhomo
        do i=1,nmo
          write(39,*) vecchk(i)
        enddo
        do i=1,nci
           write(39,*) iconf(i,1),iconf(i,2)
        enddo
        close(39)
        deallocate(vecchk)
        if (rpachk) then
          call printvecrpa(nci,nroot,uci,hci,eci) ! RPA vectors
          deallocate(hci)
        else
          call printvectda(rpachk,nci,nroot,uci,eci) ! TDA vectors
        endif
      endif

      if(nto)then
      if(rpachk)then
      call print_nto_rpa(uci,hci,ca,moci,nci,nroot,nao,iconf,
     .                   maxconf,no,nv)
      else

      call print_nto(uci,ca,moci,nci,nroot,nao,iconf,maxconf,no,nv)
      endif
      endif

      deallocate(uci)

      if(rpachk) then
      write(*,*) 'sTD-DFT-rw done.'
      else
      write(*,*) 'sTDA-rw done.'
      endif

!!! old fitting part
!c used for fitting
       inquire(file='.REF',exist=ex)
       if(.not.ex) return
       open(unit=27,file='.REF')
       open(unit=28,file='.OUT')
       read(27,*) i,k

       if(k.eq.1)then
c take only bright ones
      read(27,*) hilf
      do j=1,nroot
         if(umerk(2,j).gt.0.03)then
         write(28,'(4f12.6)')
     .   (umerk(1,j)*27.211-hilf),hilf,umerk(1,j)*27.211,umerk(2,j)
         exit
         endif
      enddo
      elseif(k.eq.2)then
c take only very bright ones
      read(27,*) hilf
      do j=1,nroot
         if(umerk(2,j).gt.0.2)then
         write(28,'(4f12.6)')
     .   (umerk(1,j)*27.211-hilf),hilf,umerk(1,j)*27.211,umerk(2,j)
         exit
         endif
      enddo
      elseif(k.eq.3)then
c take only very bright CD ones
      read(27,*) hilf
      do j=1,nroot
         if(abs(umerk(5,j)).gt.100)then
         write(28,'(4f12.6)')
     .   (umerk(1,j)*27.211-hilf),hilf,umerk(1,j)*27.211,umerk(5,j)
         exit
         endif
      enddo
      else
      do j=1,i
         read(27,*) hilf
         if(hilf.gt.0.01)
     .   write(28,'(4f12.6)')
     .   (umerk(1,j)*27.211-hilf),hilf,umerk(1,j)*27.211,umerk(2,j)
      enddo
      endif

      deallocate(iconf,umerk)

      close(27)
      close(28)

      end





***********************************************************************
* select important CSFs beyond energy threshold
***********************************************************************
      subroutine ptselect_rw(nex,ncent,no,nv,nexpt,mxcnf,iconf,kconf,
     .                    dak,dax,ed,edpt,gamj,gamk,thrp,new,moci)
      use omp_lib
      implicit none
      integer, intent(in) :: nex,ncent,nexpt,mxcnf
      integer, intent(in) :: kconf(mxcnf,2),no,nv
      integer, intent (inout) :: iconf(mxcnf,2)
      integer, intent (out) :: new
      real*8,  intent(in) :: gamj(1:ncent,1:ncent),gamk(1:ncent,1:ncent)
      real*8, intent(in)  :: dak,dax,thrp,edpt(mxcnf)
      real*8, intent(inout) :: ed(mxcnf)
      real*4, allocatable :: qia(:,:),pia(:,:)
      real*4, allocatable :: qab(:,:),qij(:,:),pij(:,:)
      real*4, allocatable :: q1(:),q2(:),q3(:)
      real*8, allocatable :: pt(:),pt2(:)
      integer i,j,k,io,iv,jo,jv,ierr,lin,iwrk,jwrk,iiv,jjv,ij,moci,l
      real*4 ek,ej,sdot,tmpi,tmpj,ak,ax,integral
      real*8 de,pert,amat
      logical, allocatable :: incl_conf(:)
      ak=real(dak)
      ax=real(dax)
      allocate(q1(ncent),q2(ncent),q3(ncent),pt(nex),pt2(nex),
     .         incl_conf(nexpt),stat=ierr)
      if(ierr.ne.0)stop 'allocation for PT intermediates crashed'
      incl_conf=.false.


      allocate(qia(ncent,mxcnf),qab(ncent,nv*(nv+1)/2),
     .         pij(ncent,no*(no+1)/2),pia(ncent,mxcnf))

      open(unit=710,file='pij',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      open(unit=73,file='qab',form='unformatted',status='old')
      open(unit=74,file='qia',form='unformatted',status='old')
      open(unit=740,file='pia',form='unformatted',status='old')
      Do i=1, no
      Do j=1, i
      ij=lin(i,j)
      read(710)pij(1:ncent,ij)
      enddo
      enddo
      close(710)
      Do i=no+1, moci
      k=i-no
      Do j=no+1, i-1
      l=j-no
      ij=lin(k,l)
      read(73)qab(1:ncent,ij)
      enddo
      ij=lin(k,k)
      read(72)qab(1:ncent,ij)
      enddo
      close(72)
      close(73)
      Do i=1, no
      Do j=no+1, moci
      ij=(i-1)*nv+j-no
      read(74)qia(1:ncent,ij)
      read(740)pia(1:ncent,ij)
      enddo
      enddo
      close(74)
      close(740)

      new=0
      pt2=0.0d0


! loop over secondary/neglected CSF, this is done omp-parallel
!$omp parallel private(i,k,io,iv,iiv,iwrk,q1,de,j,jo,jv,jjv,jwrk,ek,ej,
!$omp&                 q2,q3,pt,amat,pert) reduction (+:pt2)
!$omp do
      do k=1,nexpt
         io=kconf(k,1)
         iv=kconf(k,2)
         iiv=iv-no
         iwrk=(io-1)*nv + iiv
         q1(1:ncent)=pia(1:ncent,iwrk)
         de=edpt(k)
c loop over primary CSF
         do j=1,nex
            jo=iconf(j,1)
            jv=iconf(j,2)
            jjv=jv-no
            jwrk=(jo-1)*nv + jjv
            q2(1:ncent)=qia(1:ncent,jwrk)
            ek=sdot(ncent,q1,1,q2,1)
c coupling
            q2(1:ncent)=pij(1:ncent,lin(io,jo))
            q3(1:ncent)=qab(1:ncent,lin(iiv,jjv))
            ej=sdot(ncent,q2,1,q3,1)
            pt(j)=0.0d0
            amat=ak*ek-ej
            pt(j)=amat**2/(de-ed(j)+1.d-10)
         enddo
         pert=sum(pt(1:nex))
c if sum > threshold include it
         if(pert.gt.thrp)then
            incl_conf(k)=.true.
         else
            do i=1,nex
               pt2(i)=pt2(i)+pt(i)
            enddo
         endif
      enddo
!$omp end do
!$omp end parallel
      deallocate(qia,qab,pij,pia)
      do i=1,nexpt
         if(.not.incl_conf(i)) cycle
            io=kconf(i,1)
            iv=kconf(i,2)
            de=edpt(i)
            new=new+1
            iconf(nex+new,1)=io
            iconf(nex+new,2)=iv
            ed   (nex+new  )=de
      enddo
      pert=0.0d0
      do i=1,nex
         pert=pert+pt2(i)
         ed(i)=ed(i)-pt2(i)
      enddo

      write(*,'('' average/max PT2 energy lowering (eV):'',2F10.3)')
     .             27.21139*pert/float(nex),maxval(pt2)*27.21139


      deallocate(q1,q2,q3,pt,pt2,incl_conf)
      return

      end subroutine ptselect_rw
***********************************************************************


***********************************************************************
* set up A+B and A-B (packed form) in RKS case
***********************************************************************
      subroutine rrpamat_rw(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,
     .                  gamj,gamk,apb,ambsqr,moci)
      use commonlogicals
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*8,  intent(in) :: gamj(1:ncent,1:ncent),gamk(1:ncent,1:ncent)
      real*4, allocatable :: qia(:,:),pia(:,:)
      real*4, allocatable :: qab(:,:),qij(:,:),pij(:,:)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
      real*4, allocatable :: q1(:),q2(:),q3(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,l,moci,k
      real*4 ek,ej,sdot,ak,ax,de,integral
      allocate(q1(ncent),q2(ncent),q3(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for q1 crashed'
      ak=real(dak)
      ax=real(dax)
! calculate A+B and A-B
      apb=0.0e0
      ambsqr=0.0e0
      allocate(qab(ncent,nv*(nv+1)/2),
     .         pij(ncent,no*(no+1)/2))

      open(unit=710,file='pij',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      open(unit=73,file='qab',form='unformatted',status='old')
      open(unit=74,file='qia',form='unformatted',status='old')
      open(unit=740,file='pia',form='unformatted',status='old')
      Do i=1, no
      Do j=1, i
      ij=lin(i,j)
      read(710)pij(1:ncent,ij)
      enddo
      enddo
      close(710,status='delete')
      Do i=no+1, moci
      k=i-no
      Do j=no+1, i-1
      l=j-no
      ij=lin(k,l)
      read(73)qab(1:ncent,ij)
      enddo
      ij=lin(k,k)
      read(72)qab(1:ncent,ij)
      enddo
      close(72,status='delete')
      close(73,status='delete')
      if(abs(dax).ge.1.0d-6) then
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,jjv,
!$omp&  q2,q3,ej)
!$omp do
        do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              q2(1:ncent)=pij(1:ncent,lin(io,jo))
              q3(1:ncent)=qab(1:ncent,lin(iiv,jjv))
              ej=sdot(ncent,q2,1,q3,1) !  ej = (ij|ab)
              ambsqr(ij)=-ej
              apb(ij)=-ej
           enddo ! j
        enddo ! i
!$omp end do
!$omp end parallel
      endif
      deallocate(pij,qab)
      allocate(qia(ncent,mxcnf),
     .         pia(ncent,mxcnf))
      Do i=1, no
      Do j=no+1, moci
      ij=(i-1)*nv+j-no
      read(74)qia(1:ncent,ij)
      read(740)pia(1:ncent,ij)
      enddo
      enddo
      close(74,status='delete')
      close(740,status='delete')


! if ax=0, A-B is diagonal and its off-diagonal elements do not need to be calculated
      if(abs(dax).lt.1.0d-6) then
        ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,
!$omp&  jwrk,q1,q2,ek,de)
!$omp do
        do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           q1(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv+jjv
              q2(1:ncent)=qia(1:ncent,jwrk)
              ek=sdot(ncent,q1,1,q2,1)! ek = (ia|jb)
              apb(ij)=2.0*ak*ek
              ambsqr(ij)=0.0
           enddo ! j
           de=real(ed(i))
           ij=lin(i,i)
           q2(1:ncent)=qia(1:ncent,iwrk)
           ek=sdot(ncent,q1,1,q2,1)
           if(aresp.or.resp.or.optrota) then
           ambsqr(ij)=de-ak*ek ! diagonal element of (A-B)
           else
           ambsqr(ij)=sqrt(de-ak*ek) ! diagonal element of (A-B)^0.5
           endif
           apb(ij)=de+ak*ek       ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel
        open(unit=53,file='amb',form='unformatted',status='replace')
        write(53) ambsqr
        close(53)
      else
        ij=0
        ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,
!$omp&  jwrk,q1,q2,q3,ek)
!$omp do
        do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           q1(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv+jjv
              q2(1:ncent)=qia(1:ncent,jwrk)
              ek=sdot(ncent,q1,1,q2,1) ! ek = (ia|jb)
              ambsqr(ij)=2.0*ak*ek+ambsqr(ij)
              jwrk=(io-1)*nv+jjv
              q2(1:ncent)=pia(1:ncent,jwrk)
              jwrk=(jo-1)*nv+iiv
              q3(1:ncent)=qia(1:ncent,jwrk)
              ek=sdot(ncent,q2,1,q3,1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ambsqr(ij)=ambsqr(ij)-ax*ek
              apb(ij)=ax*ek+apb(ij)
           enddo ! j
           ij=lin(i,i)
           q2(1:ncent)=qia(1:ncent,iwrk)
           ek=sdot(ncent,q1,1,q2,1)
           apb(ij)=real(ed(i))-ak*ek+ax*ek    ! diagonal element of A-B
           ambsqr(ij)=real(ed(i))-ax*ek+ak*ek ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel
      deallocate(qia,pia)
      deallocate(q1,q2,q3)


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

      end subroutine rrpamat_rw
***********************************************************************

***********************************************************************
* set up A matrix in RKS case
***********************************************************************
      subroutine rtdamat_rw(nci,ncent,no,nv,mxcnf,iconf,dak,dax
     .                   ,ed,gamj,gamk,hci,moci)
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*8,  intent(in) :: gamj(1:ncent,1:ncent),gamk(1:ncent,1:ncent)
      real*4, allocatable :: qia(:,:),pia(:,:)
      real*4, allocatable :: qab(:,:),qij(:,:),pij(:,:)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: hci(nci,nci)
      real*4, allocatable :: q1(:),q2(:),q3(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk,l,k,moci
      real*4 ek,ej,sdot,ak,ax,de,integral
      allocate(q1(ncent),q2(ncent),q3(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj crashed'
      ak=real(dak)
      ax=real(dax)

      allocate(qia(ncent,mxcnf),qab(ncent,nv*(nv+1)/2),
     .         pij(ncent,no*(no+1)/2),pia(ncent,mxcnf))

      open(unit=710,file='pij',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      open(unit=73,file='qab',form='unformatted',status='old')
      open(unit=74,file='qia',form='unformatted',status='old')
      open(unit=740,file='pia',form='unformatted',status='old')
      Do i=1, no
      Do j=1, i
      ij=lin(i,j)
      read(710)pij(1:ncent,ij)
      enddo
      enddo
      close(710,status='delete')
      Do i=no+1, moci
      k=i-no
      Do j=no+1, i-1
      l=j-no
      ij=lin(k,l)
      read(73)qab(1:ncent,ij)
      enddo
      ij=lin(k,k)
      read(72)qab(1:ncent,ij)
      enddo
      close(72,status='delete')
      close(73,status='delete')
      Do i=1, no
      Do j=no+1, moci
      ij=(i-1)*nv+j-no
      read(74)qia(1:ncent,ij)
      read(740)pia(1:ncent,ij)
      enddo
      enddo
      close(74,status='delete')
      close(740,status='delete')

! calculate CIS matrix A
      hci=0.0e0
!$omp parallel private(i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,q1,
!$omp& q2,q3,ek,ej)
!$omp do
      do i=1,nci
         io=iconf(i,1)
         iv=iconf(i,2)
         iiv=iv-no
         iwrk=(io-1)*nv + iiv
         q1(1:ncent)=pia(1:ncent,iwrk)
         do j=1,i-1
            jo=iconf(j,1)
            jv=iconf(j,2)
            jjv=jv-no
            jwrk=(jo-1)*nv + jjv
            q2(1:ncent)=qia(1:ncent,jwrk)
            ek=sdot(ncent,q1,1,q2,1)
            q2(1:ncent)=pij(1:ncent,lin(io,jo))
            q3(1:ncent)=qab(1:ncent,lin(iiv,jjv))
            ej=sdot(ncent,q2,1,q3,1)
            hci(j,i)=ak*ek-ej
            hci(i,j)=hci(j,i)
         enddo
         hci(i,i)=real(ed(i))
      enddo
!$omp end do
!$omp end parallel
      deallocate(qia,qab,pij,pia)
      deallocate(q1,q2,q3)
      return
      end subroutine rtdamat_rw
***********************************************************************

      real*4 function integral(q_ij,q_ab,gamma,ncent)
      use omp_lib
      implicit none
      ! Compute semi-empirical integrals directly
      integer, intent(in) :: ncent
      real*8,  intent(in) :: gamma(1:ncent,1:ncent)
      real*4 :: q_ij(1:ncent), q_ab(1:ncent), intermediate(1:ncent)
      real*4 :: sdot
      call ssymv('u',ncent,1.0,gamma,ncent,q_ab,1,0.0,intermediate,1)
      integral = sdot(ncent,q_ij,1,intermediate,1)
      return
      end

***********************************************************************
* set up 0.5*B  (packed form) in RKS case
***********************************************************************
      subroutine rtdacorr_rw(nci,ncent,no,nv,mxcnf,iconf,dak,dax
     .                    ,ed,gamj,gamk,moci)
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*8,  intent(in) :: gamj(1:ncent,1:ncent),gamk(1:ncent,1:ncent)
      real*4, allocatable :: qia(:,:),pia(:,:)
      real*4, allocatable :: qab(:,:),qij(:,:),pij(:,:)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, allocatable :: q1(:),q2(:),q3(:),bmat(:)
      integer i,j,io,iv,jo,jv,ierr,iiv,jjv,iwrk,jwrk,l,k,moci
      integer*8 ij,lin8
      real*4 ek,ej,sdot,ak,ax,de,fact,integral
      ij=nci
      ij=ij*(ij+1)/2
      allocate(q1(ncent),q2(ncent),q3(ncent),bmat(ij), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj/bmat crashed'
      ak=real(dak)
      ax=real(dax)

      allocate(qia(ncent,mxcnf),qab(ncent,nv*(nv+1)/2),
     .         pij(ncent,no*(no+1)/2),pia(ncent,mxcnf))

      open(unit=710,file='pij',form='unformatted',status='old')
      open(unit=72,file='qaa',form='unformatted',status='old')
      open(unit=73,file='qab',form='unformatted',status='old')
      open(unit=74,file='qia',form='unformatted',status='old')
      open(unit=740,file='pia',form='unformatted',status='old')
      Do i=1, no
      Do j=1, i
      ij=lin8(i,j)
      read(710)pij(1:ncent,ij)
      enddo
      enddo
      close(710)
      Do i=no+1, moci
      k=i-no
      Do j=no+1, i-1
      l=j-no
      ij=lin8(k,l)
      read(73)qab(1:ncent,ij)
      enddo
      ij=lin8(k,k)
      read(72)qab(1:ncent,ij)
      enddo
      close(72)
      close(73)
      Do i=1, no
      Do j=no+1, moci
      ij=(i-1)*nv+j-no
      read(74)qia(1:ncent,ij)
      read(740)pia(1:ncent,ij)
      enddo
      enddo
      close(74)
      close(740)


! calculate 0.5*B
      bmat=0.0e0
      fact=0.50d0 ! this is the scaling of the B-contribution
      open(unit=52,file='bmat',form='unformatted',status='replace')
      ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,
!$omp&                 jwrk,q1,q2,q3,ek,ej)
!$omp do
      do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           q1(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin8(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv + jjv
              q2(1:ncent)=qia(1:ncent,jwrk)
              ek=sdot(ncent,q1,1,q2,1) ! ek = (ia|bj)
              bmat(ij)=(fact)*ak*ek
              jwrk=(io-1)*nv+jjv
              q2(1:ncent)=pia(1:ncent,jwrk)
              jwrk=(jo-1)*nv+iiv
              q3(1:ncent)=qia(1:ncent,jwrk)
              ek=sdot(ncent,q2,1,q3,1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              bmat(ij)=bmat(ij)-fact*ax*ek ! scaled by ax
           enddo
           ij=lin8(i,i)
           q2(1:ncent)=qia(1:ncent,iwrk)
           ek=sdot(ncent,q1,1,q2,1)
           bmat(ij)=fact*(ak*ek-ax*ek) ! diagonal element of 0.5*B
      enddo
!$omp end do
!$omp end parallel
      write(52)bmat
      close(52)
      deallocate(qia,qab,pij,pia)
      deallocate(bmat,q1,q2,q3)
      return

      end subroutine rtdacorr_rw
***********************************************************************
