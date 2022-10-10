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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ncent: number of atoms
C nmo  : number of MOs   
C nao  : number of contracted AOs
C xyz  : array with atomic coordinates (1:3) and nuclear charge (4) in
C        Bohr
C c    : MO coefficients (nao*nmo)
C eps  : MO energies (au)             
C occ  : occupation numbers           
C csp  : molecular orbital spin (alpha - 1, beta - 2)
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

      SUBROUTINE sutda(ncent,nmo,nao,xyz,c,eps,occ,csp,iaoat,thr,thrp,
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
      integer csp(*)
      logical ggavec,ex

c local varibles:

c l=dipole lengths, v=dipole velocity (d/dxyz), m=angular momentum      
      real*8, allocatable ::xla(:),yla(:),zla(:)
      real*8, allocatable ::xva(:),yva(:),zva(:)
      real*8, allocatable ::xma(:),yma(:),zma(:)
      real*8, allocatable ::xlb(:),ylb(:),zlb(:)
      real*8, allocatable ::xvb(:),yvb(:),zvb(:)
      real*8, allocatable ::xmb(:),ymb(:),zmb(:)
      real*8, allocatable ::help(:),scr(:),dum(:),x(:,:)
c MOs, orbital energies and CSF printout stuff      
      real*8, allocatable ::cca(:),epsia(:)
      real*8, allocatable ::ccb(:),epsib(:)
      real*8, allocatable ::umerk(:,:)
      real*8, allocatable ::umrkx(:,:),umrky(:,:),umrkz(:,:)
      real*8, allocatable ::rvp(:)

c stuff for diag of TDA matrix      
c critical regarding memory     
      integer info,lwork,liwork,il,iu,nfound
      real*4, allocatable ::uci  (:,:)
      real*4, allocatable ::eci  (:)
      real*4, allocatable ::hci  (:,:)
      real*4, allocatable ::work (:)
      real*4  vu,vl
      integer,allocatable ::iwork(:)
      integer,allocatable::isuppz(:)

c Löwdin MOs, repulsion terms, charges and half-transformed stuff
c critical regarding memory     
      real*8, allocatable ::clowa(:)
      real*8, allocatable ::clowb(:)
      real*4, allocatable ::gamj(:,:)
      real*4, allocatable ::gamk(:,:)
      real*4, allocatable ::qiaa(:,:),qiab(:,:),piaa(:,:),piab(:,:)
      real*4, allocatable ::pija(:,:),qaba(:,:),qija(:,:)
      real*4, allocatable ::pijb(:,:),qabb(:,:),qijb(:,:)
      real*4, allocatable ::q1a(:),q1b(:),q2(:)
!      real*4 q1a(ncent),q1b(ncent),q2(ncent)
! prescreening vectors K_ia,ia
      real*4 sdot

c the maximum size of the TDA expansion space
!      integer maxconf1,maxconf2                                 
!      parameter (maxconf1=500000)
!c the maximum size of the TDA pt2 space
!      parameter (maxconf2=maxconf1*10)
!      integer kconfa(maxconf2,2)
!      integer kconfb(maxconf2,2)
!      real*8  edpta(maxconf2)
!      real*8  edptb(maxconf2)
! maximum sizes of P-CSF and S-CSF expansion space
      integer maxconfa,maxconfb
      real*8, allocatable ::eda(:),edb(:)
      real*8, allocatable ::edpta(:),edptb(:)
      integer, allocatable :: iconfa(:,:),iconfb(:,:)
      integer, allocatable :: kconfa(:,:),kconfb(:,:)

c RPA stuff
      real*4, allocatable ::apb(:)
      real*4, allocatable ::ambsqr(:)

c intermediates       
      real*8 omax,vmin,pert,de,ek,ej,ak,xc,rabx,ef
      real*8 pp,hilf,uu,sss,rl,rv,rm,time,coc(3)
      real*8 fl,fv,ec,p23,xp,umax,xvu,yvu,zvu,xmu,ymu,zmu,xlu,ylu,zlu
      real*8 xj,amat,xmolw,xk,loc,jii,jaa,xms,yms,zms
      real*8 betaj,alphak,beta1,beta2,alpha1,alpha2,deps,fthr
      integer moci,i,j,k,ii,io,iv,ihilf,jo,jv,m,l,idum1,idum2
      integer*8 imem1,imem2,imem3
      integer jj,kmem,imax,jmax,lmem,jmem,ij,lin,ierr
! state/orbital indices and dimensions
      integer n,no,nv,noa,nob,nva,nvb,nuhf,nsomoa,nsomob
      integer nexpt,nexpta,nexptb
      integer nex,nexa,nexb
      integer new,newa,newb
      integer mocia,mocib
      integer ispin,nci,nroot
      ! SOS polarizability and sum rule check
      real*8 alp_real(6),sumf

! variables for vector printout
      integer, allocatable :: vecchka(:),vecchkb(:)
      integer nvec,jhomo,jhomoa,jhomob
c atomic Hubbard parameters      
      real*8 eta(94)
c atomic masses      
      common /amass/ ams(107)
      real*8 ams
      character*79 dummy
c Linear response
      real*4 :: start_time, end_time, stda_time
      integer :: STATUS

c just a printout      
      call header('s T D A',0)
      
      thr =thr /27.211385050d0
c estimate the orbital energy window which corresponds to the desired
c spectra range thr
      deps=(1.+0.8*ax)*thr
c make it safe
      deps=deps*2.0

      omax=-1d+42
      vmin= 1d+42
      do i=1,nmo
         if(occ(i).gt.0.99.and.eps(i).gt.omax) omax=eps(i)
         if(occ(i).lt.0.01.and.eps(i).lt.vmin) vmin=eps(i)
      enddo


! optional: if eigenvectors are wanted in TM format, check now how many occupied there are in general
! this is needed to get the CSF sorting of TM
      jhomo=0
      jhomoa=0
      jhomob=0
      do i=1,nmo
        if(occ(i).gt.0.990d0) then 
          jhomo=jhomo+1
          if(csp(i).eq.1)  jhomoa = jhomoa + 1
          if(csp(i).eq.2)  jhomob = jhomob + 1
        endif
      enddo

      nuhf=jhomoa-jhomob

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

      mocia = 0
      mocib = 0
      do i = 1,nmo
         if(eps(i).gt.othr.and.eps(i).lt.vthr) then
          if(csp(i).eq.1)  mocia = mocia + 1
          if(csp(i).eq.2)  mocib = mocib + 1
         endif
      enddo

      moci = mocia + mocib

      allocate(
     .         xla(mocia*(mocia+1)/2),yla(mocia*(mocia+1)/2),
     .         zla(mocia*(mocia+1)/2),
     .         xlb(mocib*(mocib+1)/2),ylb(mocib*(mocib+1)/2),
     .         zlb(mocib*(mocib+1)/2),
     .         xva(mocia*(mocia+1)/2),yva(mocia*(mocia+1)/2),
     .         zva(mocia*(mocia+1)/2),
     .         xvb(mocib*(mocib+1)/2),yvb(mocib*(mocib+1)/2),
     .         zvb(mocib*(mocib+1)/2),
     .         xma(mocia*(mocia+1)/2),yma(mocia*(mocia+1)/2),
     .         zma(mocia*(mocia+1)/2),
     .         xmb(mocib*(mocib+1)/2),ymb(mocib*(mocib+1)/2),
     .         zmb(mocib*(mocib+1)/2),
     .         help(nao*(nao+1)/2),
     .         clowa(nao*mocia),clowb(nao*mocib),
     .         scr(nao*nao),dum(nao*nao),x(nao,nao),
     .         gamj(ncent,ncent),gamk(ncent,ncent),
     .         cca(nao*mocia),epsia(mocia),
     .         ccb(nao*mocib),epsib(mocib)
     .        )

      write(*,*) 'Active space MOs in TDA : ',moci,mocia,mocib        

      mocia=0
      mocib=0
      moci =0
! optional: if eigenvectors are wanted, make set a checking variable (vecchk) that maps the 
! TM orbital space onto the sTDA orbital space (which is reduced due to thresholds)
      if(eigvec.or.nto) then
        allocate (vecchka(nmo/2),vecchkb(nmo/2))
        vecchka=0 ! for vector printout
        vecchkb=0
        l=0
        k=0
        do i = 1,nmo
           if(occ(i).gt.0.99) then 
            if(csp(i).eq.1) then
              l=l+1 ! increase alpha counter
              if(eps(i).gt.othr)then
                mocia = mocia + 1
                moci = moci +1
                do j=1,nao
                  cca(j+(mocia-1)*nao)=c(j+(i-1)*nao)
                enddo
                epsia(mocia)=eps(i)
                vecchka(l)=mocia
              endif
            endif
            if(csp(i).eq.2) then
              k = k+1 ! increase beta counter
              if(eps(i).gt.othr)then
                mocib = mocib + 1
                moci = moci + 1
                do j=1,nao
                  ccb(j+(mocib-1)*nao)=c(j+(i-1)*nao)
                enddo
                epsib(mocib)=eps(i)
                vecchkb(k)=mocib
              endif
            endif
           endif
        enddo

        no = moci
        noa = mocia
        nob = mocib

        do i = 1,nmo
           if(occ(i).lt.0.01) then 
            if(csp(i).eq.1) then
              l=l+1 ! increase alpha counter
              if(eps(i).lt.vthr)then
                mocia = mocia + 1
                moci = moci +1
                do j=1,nao
                  cca(j+(mocia-1)*nao)=c(j+(i-1)*nao)
                enddo
                epsia(mocia)=eps(i)
                vecchka(l)=mocia
              endif
            endif
            if(csp(i).eq.2) then
              k = k+1 ! increase beta counter
              if(eps(i).lt.vthr)then
                mocib = mocib + 1
                moci = moci + 1
                do j=1,nao
                  ccb(j+(mocib-1)*nao)=c(j+(i-1)*nao)
                enddo
                epsib(mocib)=eps(i)
                vecchkb(k)=mocib
              endif
            endif
           endif
        enddo

      else
! do not print eigenvectors case
        do i = 1,nmo
           if(occ(i).gt.0.99.and.eps(i).gt.othr)then
            moci = moci + 1
            if(csp(i).eq.1) then
              mocia = mocia + 1
              do j=1,nao
                cca(j+(mocia-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsia(mocia)=eps(i)
            endif
            if(csp(i).eq.2) then
              mocib = mocib + 1
              do j=1,nao
                ccb(j+(mocib-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsib(mocib)=eps(i)
            endif
           endif
        enddo
 
        no = moci
        noa = mocia
        nob = mocib
  
        do i = 1,nmo
           if(occ(i).lt.0.01.and.eps(i).lt.vthr)then          
              moci=moci+1
              if(csp(i).eq.1) then
              mocia = mocia + 1
              do j=1,nao
                cca(j+(mocia-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsia(mocia)=eps(i)
            endif
            if(csp(i).eq.2) then
              mocib = mocib + 1
              do j=1,nao
                ccb(j+(mocib-1)*nao)=c(j+(i-1)*nao)
              enddo
              epsib(mocib)=eps(i)
            endif
           endif
        enddo

      endif
      nv  = moci - no
      nva = mocia - noa
      nvb = mocib - nob

      write(*,*) 'Occupied active MOs in TDA: ', no,noa,nob           
      write(*,*) 'Virtual active MOs in TDA: ', nv,nva,nvb             
      if((noa.eq.0.or.nva.eq.0).and.(nob.eq.0.or.nvb.eq.0)) then 
        stop 'no CSF, increase energy threshold (-e option)'
      endif

      maxconfa=noa*nva
      maxconfb=nob*nvb
      allocate(eda(maxconfa),edpta(maxconfa),
     .         edb(maxconfb),edptb(maxconfb), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for ed and edpt'
      allocate(iconfa(maxconfa,2),kconfa(maxconfa,2), 
     .         iconfb(maxconfb,2),kconfb(maxconfb,2), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for iconf and kconf'
      eda=0.0d0      
      edpta=0.0d0
      iconfa=0
      kconfa=0
      edb=0.0d0
      edptb=0.0d0
      iconfb=0
      kconfb=0


c we arrange MOS according to energy from 1:HOMO to LUMO:MOCI      
c (in TM they come in irreps)
      write(*,*) 'Sorting MOs ...'
c sort for E diag

      call sort_vec(mocia,nao,cca,epsia)
      call sort_vec(mocib,nao,ccb,epsib)

      write(*,*) 'Reading and transforming R..V..L AO ints ...'

c read L,V,M with xyz components each and transform
c to MO basis (original but sorted MOs in array ca)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            dipole lengths
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=31,file='xlint',form='unformatted',status='old')
      read(31) help
      call onetri(1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,xla)
      call onetri(1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,xlb)
      close(31,status='delete')

      open(unit=32,file='ylint',form='unformatted',status='old')
      read(32) help
      call onetri(1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,yla)
      call onetri(1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,ylb)
      close(32,status='delete')

      open(unit=33,file='zlint',form='unformatted',status='old')
      read(33) help
      call onetri(1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,zla)
      call onetri(1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,zlb)
      close(33,status='delete')


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            magnetic dipole
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=34,file='xmint',form='unformatted',status='old')
      read(34) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,xma)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,xmb)
      close(34,status='delete')

      open(unit=35,file='ymint',form='unformatted',status='old')
      read(35) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,yma)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,ymb)
      close(35,status='delete')

      open(unit=36,file='zmint',form='unformatted',status='old')
      read(36) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,zma)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,zmb)
      close(36,status='delete')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c            velocity dipole
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      open(unit=37,file='xvint',form='unformatted',status='old')
      read(37) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,xva)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,xvb)
      close(37,status='delete')

      open(unit=38,file='yvint',form='unformatted',status='old')
      read(38) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,yva)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,yvb)
      close(38,status='delete')     

      open(unit=39,file='zvint',form='unformatted',status='old')
      read(39) help
      call onetri(-1,help,dum,scr,cca,nao,mocia)
      call shrink(mocia,dum,zva)
      call onetri(-1,help,dum,scr,ccb,nao,mocib)
      call shrink(mocib,dum,zvb)
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
 
      call dgemm('n','n',nao,mocia,nao,1.d0,x,nao,cca,nao,0.d0,scr,nao)

      do i=1,mocia
      sss=0
       do j=1,nao
         sss = sss+scr(j+(i-1)*nao)**2
         clowa(j+(i-1)*nao) = scr(j+(i-1)*nao)
       enddo
       if(abs(sss-1.).gt.1.d-2)then
         write(*,*) 'alpha MO norm ',i,sss
         stop 'internal MO norm error'
       endif
      enddo

      call dgemm('n','n',nao,mocib,nao,1.d0,x,nao,ccb,nao,0.d0,scr,nao)

      do i=1,mocib
      sss=0
       do j=1,nao
         sss = sss+scr(j+(i-1)*nao)**2
         clowb(j+(i-1)*nao) = scr(j+(i-1)*nao)
       enddo
       if(abs(sss-1.0d0).gt.1.d-2)then
         write(*,*) 'beta MO norm ',i,sss
         stop 'internal MO norm error'
       endif
      enddo

      write(*,*) 'S^1/2 orthogonalized MO coefficients done.'

      deallocate(scr,dum,help,x)

!      q1a = 0
!      do i=1,noa
!         call lo12pop(i,i,ncent,nao,iaoat,clowa,q2)
!         q1a(1:ncent)=q1a(1:ncent)+q2(1:ncent)
!      enddo
!      q1b = 0
!      do i=1,nob
!         call lo12pop(i,i,ncent,nao,iaoat,clowb,q2)
!         q1b(1:ncent)=q1b(1:ncent)+q2(1:ncent)
!      enddo
      call cpu_time(time)
      allocate(q1a(ncent),q1b(ncent),q2(ncent))
      q1a=0.0
      q1b=0.0
      q2=0.0
      allocate(qija(ncent,noa),qaba(ncent,nva),
     .         qijb(ncent,nob),qabb(ncent,nvb),stat=ierr)
      if(ierr.ne.0) stop 'error in diag. J charges allocation'
      qija=0.0
      qaba=0.0
      qijb=0.0
      qabb=0.0

      do i=1,noa
         call lo12pop(i,i,ncent,nao,iaoat,clowa,q2)
         q1a(1:ncent)=q1a(1:ncent)+q2(1:ncent)
         qija(1:ncent,i)=q2(1:ncent)
      enddo
      do i=1,nob
         call lo12pop(i,i,ncent,nao,iaoat,clowb,q2)
         q1b(1:ncent)=q1b(1:ncent)+q2(1:ncent)
         qijb(1:ncent,i)=q2(1:ncent)
      enddo
      write(*,'(/'' SCF atom population (using active MOs):'')')
      write(*,'(10F7.3)') q1a(1:ncent)
      write(*,'(10F7.3)') q1b(1:ncent)
      write(*,'(10F7.3)') (q1a(1:ncent)+q1b(1:ncent))
      write(*,*)
      write(*,'('' # of alpha and beta electrons in TDA:'',2F8.3)') 
     .sum(q1a(1:ncent)),sum(q1b(1:ncent))
      write(*,*)
      do i=noa+1,mocia
         j=i-noa
         call lo12pop(i,i,ncent,nao,iaoat,clowa,q2)
         qaba(1:ncent,j)=q2(1:ncent)
      enddo
      do i=nob+1,mocib
         j=i-nob
         call lo12pop(i,i,ncent,nao,iaoat,clowb,q2)
         qabb(1:ncent,j)=q2(1:ncent)
      enddo

c  UCIS scaling factor for K(ia|jb):
      ak = 1.0d0

c the global parameters of the method:
      beta1 = 0.20d0              
      beta2 = 1.830d0              
      alpha1 = 1.420d0              
      alpha2 = 0.480d0    

      if(betaj.lt.-99.0d0) then ! if no beta parameter was read in
      betaj=beta1+beta2*ax
      endif
      if(alphak.lt.-99.0d0) then ! if no alpha parameter was read in
      alphak=alpha1+alpha2*ax
      endif
 
      write(*,*)
      write(*,*) 'ax(DF)  : ',ax
      write(*,*) 's^K      : ',ak
      write(*,*) 'beta  (J): ',betaj
      write(*,*) 'alpha (K): ',alphak
      write(*,*)

c set the hardness table
      call setrep(eta)
      write(*,*) 'hardness table read.'

      write(*,*) 'setting up gammas ...'
c compute gamma(j/k)
      xmolw=0
      do i = 1,ncent
         ii = idint(xyz(4,i))
c ams is the atomic mass (mol weight for output file)        
         xmolw = xmolw+ams(ii)
         do j=1,i      
            jj = idint(xyz(4,j))
            xj  = 0.50d0*(eta(ii)+eta(jj)) * ax 
            xk  = 0.50d0*(eta(ii)+eta(jj)) 
            rabx = dsqrt((xyz(1,i)-xyz(1,j))**2
     .                 +(xyz(2,i)-xyz(2,j))**2
     .                 +(xyz(3,i)-xyz(3,j))**2)
            gamj(j,i)=real(1./(rabx**betaj+1./xj**betaj)**(1.0d0/betaj))
            gamk(j,i)=real(1./(rabx**alphak+1./xk**alphak)
     .                **(1.0d0/alphak))
            gamj(i,j) = gamj(j,i)
            gamk(i,j) = gamk(j,i)
         enddo
      enddo


! compute intermediates q's refer to charges, p's to gam*q (i.e., contracted)
      imem1=noa
      imem2=nva
      imem3=(noa*nva+nob*nvb)*2
      imem1=noa*(noa+1)/2+nob*(nob+1)/2
      imem2=nva*(nva+1)/2+nvb*(nvb+1)/2
      imem3=imem3+imem1+imem2
      hilf=dble(imem3)/1024.0**2
      hilf=dble(4*ncent)*hilf
      imem1=idint(hilf)
      write(*,*)'memory needed for qa and qb data (Mb) ',imem1
      write(*,*)'computing q(ij,n) ...'

! first alpha part

!     Coulomb type terms: calc qij*gam^J
      allocate(pija(ncent,noa), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii| intermediate'
      pija=0.0
      call ssymm('l','l',ncent,noa,1.0,gamj,ncent,qija,ncent,0.0,pija
     .           ,ncent)
      allocate(uci(nva,noa), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii|aa) matrix'
      uci=0.0
      ! now calc (ii|aa)^J
      call sgemm('t','n',nva,noa,ncent,1.0,qaba,ncent,pija,ncent,0.0,
     .           uci,nva)
      deallocate(pija)

      nex=noa*nva
      allocate(piaa(ncent,nex),qiaa(ncent,nex), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for intermediates'
      qiaa=0.0
      piaa=0.0
!       K-type terms      
      k=0
!$omp parallel private(i,j,k,q2)
!$omp do
      do i=1,noa
         do j=noa+1,mocia
!            k=k+1           
            k=(i-1)*nva+j-noa
            call lo12pop(i,j,ncent,nao,iaoat,clowa,q2)
            qiaa(1:ncent,k)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      piaa=0.0
      call ssymm('l','l',ncent,nex,1.0,gamk,ncent,qiaa,ncent,0.0
     .             ,piaa,ncent)



c determine singles which are lower than thr using A(ia,ia)  
c alpha excitations first:
 
      k=0
      j=0
      l=0
      nsomoa=noa-nuhf
      do io=1,noa
         do iv=noa+1,mocia
c compute A(ia,ia) = eps(a)-eps(i) + (ia|ia) - (ii|aa)
            de=epsia(iv)-epsia(io)
            l=iv-noa
            ej=0.0d0
            ej=dble(uci(l,io))
            de=de-ej
            ek=0.0d0
            i=nva*(io-1)+l
            q1a(1:ncent)=piaa(1:ncent,i)
            ek=sdot(ncent,q1a,1,qiaa(1,i),1)
            de=de+ak*ek

            if(dokshift) then
! perform K(ia,ia) dependent shift
              call kshift_to_ediag(de,ek)
              call somo_shift(nsomoa,io,iv,de,ek)
            endif

c the primary ones
            if(de.le.thr)then
               k=k+1
               iconfa(k,1)=io
               iconfa(k,2)=iv
               eda(k)=de
            endif
c for PT            
            if(de.gt.thr.and.de.lt.fthr)then ! needs to be on if fthr is specified
               j=j+1
               kconfa(j,1)=io
               kconfa(j,2)=iv
               edpta(j)=de
            endif

         enddo
      enddo
      nexa=k
      nexpta=j
      deallocate(uci)

! now beta part

!     Coulomb type terms: calc qij*gam^J
      allocate(pijb(ncent,nob), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii| intermediate'
      pijb=0.0
      call ssymm('l','l',ncent,nob,1.0,gamj,ncent,qijb,ncent,0.0,pijb
     .           ,ncent)
      allocate(uci(nvb,nob), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii|aa) matrix'
      uci=0.0
      ! now calc (ii|aa)^J
      call sgemm('t','n',nvb,nob,ncent,1.0,qabb,ncent,pijb,ncent,0.0,
     .           uci,nvb)
      deallocate(pijb)

      nex=nob*nvb
      allocate(piab(ncent,nex),qiab(ncent,nex), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for intermediates'
      qiab=0.0
      piab=0.0
!       K-type terms      
      k=0
!$omp parallel private(i,j,k,q2)
!$omp do
      do i=1,nob
         do j=nob+1,mocib
!            k=k+1           
            k=(i-1)*nvb+j-nob
            call lo12pop(i,j,ncent,nao,iaoat,clowb,q2)
            qiab(1:ncent,k)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      piab=0.0
      call ssymm('l','l',ncent,nex,1.0,gamk,ncent,qiab,ncent,0.0
     .             ,piab,ncent)
      deallocate(gamk)

      k=0
      j=0
      l=0
      nsomob=nob+nuhf
      do io=1,nob
         do iv=nob+1,mocib
c compute A(ia,ia) = eps(a)-eps(i) + (ia|ia) - (ii|aa)
            de=epsib(iv)-epsib(io)
            l=iv-nob
            ej=0.0d0
            ej=dble(uci(l,io))
            de=de-ej
            ek=0.0d0
            i=nvb*(io-1)+l
            q1b(1:ncent)=piab(1:ncent,i)
            ek=sdot(ncent,q1b,1,qiab(1,i),1)
            de=de+ak*ek

            if(dokshift) then
! perform K(ia,ia) dependent shift
              call kshift_to_ediag(de,ek)
              call somo_shift(nsomob,io,iv,de,ek)
            endif

c the primary ones
            if(de.le.thr)then
               k=k+1
               iconfb(k,1)=io
               iconfb(k,2)=iv
               edb(k)=de
            endif
c for PT           
            if(de.gt.thr.and.de.lt.fthr)then ! needs to be on if fthr is specified 
               j=j+1
               kconfb(j,1)=io
               kconfb(j,2)=iv
               edptb(j)=de
            endif

         enddo
      enddo
      nexb=k
      nexptb=j
      deallocate(uci) 
      nex = nexa + nexb
      nexpt = nexpta + nexptb

      write(*,*)
      write(*,*) '# CSF included by energy:',nexa,nexb,nex
      write(*,*)
      write(*,*) '# CSF considered in PT2:',nexpta,nexptb,nexpt

c errors and warning
      if(nex.lt.1) stop 'No CSF, increase energy threshold (-e option)'
      if(nexa.lt.1.or.nexb.lt.1) 
     .   write(*,*) 'WARNING: No alpha or beta excitations included!'
!      if(nexa.eq.maxconf1.or.nexb.eq.maxconf1) 
!     .   stop 'Primary CSF space exceeded. use -e option!'
!      if(nexpta.eq.maxconf2.or.nexptb.eq.maxconf2) 
!     .   write(*,*)'CSF PT2 space exceeded. try -p option!'
  
c sort for E diag in each spin manifold
      do 141 ii = 2,nexa                                      
         i = ii - 1                                         
         k = i                                             
         pp= eda(i) 
         do 121 j = ii, nexa                               
            if (eda(j) .ge. pp) go to 121    
            k = j                                        
            pp=eda(j)
  121    continue                            
         if (k .eq. i) go to 141            
         eda(k) = eda(i)   
         eda(i) = pp       
         do m=1,2 
          ihilf=iconfa(i,m)        
          iconfa(i,m)=iconfa(k,m)         
          iconfa(k,m)=ihilf   
         enddo
  141 continue 

      do 142 ii = 2,nexb
         i = ii - 1
         k = i
         pp= edb(i)
         do 122 j = ii, nexb
            if (edb(j) .ge. pp) go to 122
            k = j
            pp=edb(j)
  122    continue
         if (k .eq. i) go to 142
         edb(k) = edb(i)
         edb(i) = pp
         do m=1,2
          ihilf=iconfb(i,m)
          iconfb(i,m)=iconfb(k,m)
          iconfb(k,m)=ihilf
         enddo
  142 continue         

c just printout  
      write(*,*) ' '       
      write(*,*)'Ordered frontier alpha orbitals:'
      write(*,*)'        eV     # centers'
      j=max(1,noa-10)
      do i=j,noa
         xc=0
         do k=1,ncent
            xc=xc+qija(k,i)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsia(i)*27.21139,1./(xc+1.d-8) 
      enddo
!      write(*,*)'Unoccupied:'
      write(*,*) ' '
      j=min(mocia,noa+11)
      do i=noa+1,j
         xc=0
         do k=1,ncent
            xc=xc+qaba(k,i-noa)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsia(i)*27.21139,1./(xc+1.d-8) 
      enddo

      write(*,*) ' ' 
      write(*,*)'Ordered frontier beta orbitals:'
      write(*,*)'        eV     # centers'
      j=max(1,nob-10)
      do i=j,nob
         xc=0
         do k=1,ncent
            xc=xc+qijb(k,i)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsib(i)*27.21139,1./(xc+1.d-8)
      enddo
!      write(*,*)'Unoccupied:'
      write(*,*) ' '
      j=min(mocib,nob+11)
      do i=nob+1,j
         xc=0
         do k=1,ncent
            xc=xc+qabb(k,i-nob)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsib(i)*27.21139,1./(xc+1.d-8) 
      enddo

      write(*,*)
      write(*,*)'            Lowest CSF states:'
      write(*,*)'            Alpha Excitations:'
      write(*,*)'      eV     nm      excitation i->a               eV'
      do i=1,min(nexa,25)
         io=iconfa(i,1)
         iv=iconfa(i,2)
         q1a(1:ncent)=qija(1:ncent,io)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1a,1,0.0,q2,1)
         jii=sdot(ncent,q1a,1,q2,1)
         k=iv-noa
         q1a(1:ncent)=qaba(1:ncent,k)
         ej=sdot(ncent,q1a,1,q2,1)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1a,1,0.0,q2,1)
         jaa=sdot(ncent,q1a,1,q2,1)
         l=nva*(io-1)+k
         q1a(1:ncent)=piaa(1:ncent,l)
         q2(1:ncent)=qiaa(1:ncent,l)
         ek=sdot(ncent,q1a,1,q2,1)
! de is now the Kia shift
         de=0
         loc=ej/sqrt(jii*jaa) ! locality
         if(dokshift) call kshift_to_ediag(de,ek)
         if(dokshift) call somo_shift(nsomoa,io,iv,de,ek)
         write(*,16) i,27.211*eda(i),
     .   1.d+7/(eda(i)*2.19474625d+5),iconfa(i,1:2),
     .   27.211*(epsia(iv)-epsia(io)),27.211*ej,27.211*ek,27.211*de,loc
      enddo
 16   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3, 
     .       3x,'Kshft:',f8.3,2x,'locality:',f6.3)
      write(*,*)
      write(*,*)'            Beta Excitations:'
      write(*,*)'      eV     nm      excitation i->a               eV'
      do i=1,min(nexb,25)
         io=iconfb(i,1)
         iv=iconfb(i,2)
         q1b(1:ncent)=qijb(1:ncent,io)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1b,1,0.0,q2,1)
         jii=sdot(ncent,q1b,1,q2,1)
         k=iv-nob
         q1b(1:ncent)=qabb(1:ncent,k)
         ej=sdot(ncent,q1b,1,q2,1)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1b,1,0.0,q2,1)
         jaa=sdot(ncent,q1b,1,q2,1)
         l=nvb*(io-1)+k
         q1b(1:ncent)=piab(1:ncent,l)
         q2(1:ncent)=qiab(1:ncent,l)
         ek=sdot(ncent,q1b,1,q2,1)
! de is now the Kia shift
         de=0
         loc=ej/sqrt(jii*jaa) ! locality
         if(dokshift) call kshift_to_ediag(de,ek)
         if(dokshift) call somo_shift(nsomob,io,iv,de,ek)
          write(*,16) i,27.211*edb(i),
     .   1.d+7/(edb(i)*2.19474625d+5),iconfb(i,1:2),
     .   27.211*(epsib(iv)-epsib(io)),27.211*ej,27.211*ek,27.211*de,loc
      enddo

      deallocate(qija,qaba,qijb,qabb)
!!! now set up pij an qab -alpha first
      ihilf=noa*(noa+1)/2
      allocate(qija(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qij allocation'
      qija=0.0
      ij=0
!$omp parallel private(i,j,ij,q2)
!$omp do
      do i=1,noa
         do j=1,i
            ij=lin(i,j)
            call lo12pop(i,j,ncent,nao,iaoat,clowa,q2)
            qija(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel 
      allocate(pija(ncent,ihilf), stat=ierr)
        if(ierr.ne.0)stop 'allocation failed for (ij| intermediate'
        pija=0.0
        call ssymm('l','l',ncent,ihilf,1.0,gamj,ncent,qija,ncent,0.0
     .             ,pija,ncent)

      deallocate(qija)

      q2=0.0
      ihilf=nva*(nva+1)/2
      allocate(qaba(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qab allocation'
      ij=0
!$omp parallel private(i,j,l,k,ij,q2)
!$omp do
      do i=noa+1,mocia
         do j=noa+1,i
           k=i-noa
           l=j-noa
           ij=lin(k,l)
           call lo12pop(i,j,ncent,nao,iaoat,clowa,q2)
           qaba(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      deallocate(clowa)
!!! now beta part of pij and qab !!!
      ihilf=nob*(nob+1)/2
      allocate(qijb(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qij allocation'
      qijb=0.0
      ij=0
!$omp parallel private(i,j,ij,q2)
!$omp do
      do i=1,nob
         do j=1,i
            ij=lin(i,j)
            call lo12pop(i,j,ncent,nao,iaoat,clowb,q2)
            qijb(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel 
      allocate(pijb(ncent,ihilf), stat=ierr)
        if(ierr.ne.0)stop 'allocation failed for (ij| intermediate'
        pijb=0.0
        call ssymm('l','l',ncent,ihilf,1.0,gamj,ncent,qijb,ncent,0.0
     .             ,pijb,ncent)

      deallocate(gamj,qijb)

      q2=0.0
      ihilf=nvb*(nvb+1)/2
      allocate(qabb(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qab allocation'
      ij=0
!$omp parallel private(i,j,l,k,ij,q2)
!$omp do
      do i=nob+1,mocib
         do j=nob+1,i
           k=i-nob
           l=j-nob
           ij=lin(k,l)
           call lo12pop(i,j,ncent,nao,iaoat,clowb,q2)
           qabb(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel

      deallocate(clowb)

      write(*,*)
      write(*,*)'selecting CSF ...'
c determine singles which interact most strongly with all basis
c CSF which have been selected by the  diagonal element A(ia,ia)
      deallocate(q1a,q1b,q2)
      call ptselect_uks(nexa,nexb,ncent,noa,nva,nob,nvb,nexpta,nexptb,
     .                  maxconfa,maxconfb,iconfa,iconfb,kconfa,kconfb,
     .                  ak,ax,eda,edb,edpta,edptb,piaa,qiaa,piab,qiab,
     .                  pija,qaba,pijb,qabb,thrp,newa,newb)

      deallocate(edpta,kconfa,edptb,kconfb)

      new = newa + newb


c nroot at this point is the number of primary CSF. The
c number of roots in the range 0-thr (as desired) is not
c known but will be determined by the diag routine. 
      nroot = nex
      write(*,*) 'CSF included by PT:',new,newa,newb
      nexa = nexa + newa
      nexb = nexb + newb
      nex= nexa + nexb
      nci = nex
      write(*,*) 'CSF in total:',nex,nexa,nexb

      if(rpachk) then
!**************************
!                         *
!    sTD-DFT   procedure  *
!                         *
!**************************
********************************************************************************
      write(*,*) 'sTD-DFT procedure...'
      write(*,*) 'setting up A+B and A-B matrices'
********************************************************************************
      allocate(ambsqr(nci*(nci+1)/2),stat=ierr)
      call urpamat(nci,nexa,nexb,ncent,noa,nva,nob,nvb,maxconfa,
     .             maxconfb,iconfa,iconfb,ax,eda,edb,piaa,qiaa,piab,
     .             qiab,pija,qaba,pijb,qabb,ambsqr)
! get rid of large arrays
      deallocate(pija,qaba,piaa,eda,pijb,qabb,piab,edb)
      if(.not.printexciton)deallocate(qiaa,qiab)

! A+B was computed in routine and written to file 
! read it now
      allocate(apb(nci*(nci+1)/2),stat=ierr)
      open(unit=52,file='apbmat',form='unformatted',status='old') 
      read(52) apb
      close(52,status='delete')
********************************************************************************

! allocate uci and hci as X and Y vectors respectively 
      allocate(hci(nci,nci),uci(nci,nci),stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for X and Y matrix'
! allocate eigenvalue vector
      allocate( eci(nci) ,stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for eci vector'

**************************************************************************
! if eigenvectors in TM format wanted for GGAs, this needs to be done here
      ggavec=.false.
      if (ax.eq.0.0d0.and.eigvec) then
        open(unit=39,file='TmPvEcInFo',status='replace')
        ! print to temporary file
        write(39,*) 2
        write(39,*) nvec,nmo,jhomo
        write(39,*) nexa,nexb,jhomoa,jhomob
        do i=1,nmo/2
          write(39,*) vecchka(i)
        enddo
        do i=1,nmo/2
          write(39,*) vecchkb(i)
        enddo
        do i=1,nexa
           write(39,*) iconfa(i,1),iconfa(i,2)
        enddo
        do i=1,nexb
           write(39,*) iconfb(i,1),iconfb(i,2)
        enddo
        close(39)
        ggavec=.true.
        eigvec=.false.
        deallocate(vecchka,vecchkb)
      endif
************************************************************************

c call sRPA routine
c uci is now X, hci is Y

      call srpapack(nci,thr,ambsqr,apb,eci,uci,hci,info,ggavec)
      if(info.lt.1) stop 'internal error in diag'
      nroot=info
      info=0
c      call prmat4(6,hci,nci,nci,'hci')
c      call prmat4(6,uci,nci,nci,'uci')
      deallocate(ambsqr,apb)

      else
cccccccccccccccccccccccccccc
c                          c
c Standard TDA procedure   c
c                          c
cccccccccccccccccccccccccccc

c allocate the TDA matrix
      allocate(hci(nci,nci),stat=ierr)
      if(ierr.ne.0) stop 'allocation failed for TDA matrix'

********************************************************************************
c set the TDA matrix up                                                        *
      write(*,*)'calculating TDA matrix ...'              
********************************************************************************
      call utdamat(nci,nexa,nexb,ncent,noa,nva,nob,nvb,maxconfa,
     .             maxconfb,iconfa,iconfb,ax,eda,edb,piaa,qiaa,piab,
     .             qiab,pija,qaba,pijb,qabb,hci)
********************************************************************************

!********************************************************************************
! construct ( 0.5 * B ) for X trafo (velocity correction) and print to file *
!********************************************************************************
      if(velcorr) then
       call utdacorr(nexa,nexb,ncent,noa,nva,nob,nvb,maxconfa,maxconfb,
     .               iconfa,iconfb,ax,piaa,qiaa,piab,qiab,pija,qaba,
     .               pijb,qabb)
      endif
!********************************************************************************

c big arrays not needed anymore
      deallocate(pija,qaba,piaa,eda,pijb,qabb,piab,edb)
      if(.not.printexciton)deallocate(qiaa,qiab)

c Diagonalize hci (A-matrix)

      write(*,*) 'Diagonalizing ...'
      write(*,'('' estimated time (min) '',f8.2)')
     .            float(nci)**2*float(nroot)/8.d+8/60.

c if LAPACK does not work      
c     allocate(eci(nci),uci(nci,nroot),
c    .         stat=ierr)
c     call sHQRII(hci,nci,nroot,eci,uci)

c faster by a factor of 2-3
      lwork =26*nci 
      liwork=10*nci
c we allocate uci with nci (and not with nroot) as safe 
c choice (other values gave segfaults)
      nroot=min(nci,int(1.5*nroot))
      allocate(eci(nci),uci(nci,nci),work(lwork),
     .         iwork(liwork),isuppz(nci),stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for TDA matrix diag'

      vl=0  
      vu=thr
      call ssyevr('V','V','U',nci,hci,nci,vl,vu,il,iu,1.e-6,
     .            nfound,eci,uci,nci,isuppz,
     .            work,lwork,iwork,liwork,info)
      nroot=nfound
      if(nfound.lt.1) stop 'internal error in diag'

      endif ! if statement to differentiate between TDA and RPA

      write(*,'(i5,'' roots found, lowest/highest eigenvalue : '',
     .2F8.3,i4)') nroot,eci(1)*27.21139,eci(nroot)*27.21139,info
      if(info.gt.0) stop 'diag error (ssyevr)'

      allocate(umerk(17,nroot))

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
          write(27,dummy) xyz(1:3,i)-coc(1:3),nint(xyz(4,i)) ! print molecular geometry with origin shifted to c.o.c.
        enddo
        write(27,'(a)',advance='yes')'$states'
      endif

      allocate(q2(ncent))
      allocate(rvp(nroot),stat=ierr)
      if(ierr.ne.0) stop 'error in rvp allocation'

      alp_real=0.0d0
      sumf=0.0d0

c contract WF with MO integrals for intensities
       p23=2.0d0/3.0d0

! if aniso, store and print the x,y,z-resolved data
      if(.not.velcorr.and.aniso) then
        allocate(umrkx(4,nroot),umrky(4,nroot),umrkz(4,nroot),stat=ierr)
        if(ierr.ne.0) stop 'error in xyz-umerk alloc'
        umrkx=0.0d0
        umrky=0.0d0
        umrkz=0.0d0
      endif

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
         umax=-1   
         kmem=1
         ispin=1
         q2=0.0
c loop over basis P-CSF
c alpha P-CSF 
         do j=1,nexa
            io=iconfa(j,1)
            iv=iconfa(j,2)
            ij=lin(io,iv)
            uu=dble(uci(j,i))
            pp=dble(uci(j,i)) 
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax)then
               umax=uu*pp
               imax=io
               jmax=iv
               kmem=j
               ispin=1
            endif
            xlu=xlu+xla(ij)*uu
            ylu=ylu+yla(ij)*uu
            zlu=zlu+zla(ij)*uu
            xvu=xvu+xva(ij)*pp
            yvu=yvu+yva(ij)*pp
            zvu=zvu+zva(ij)*pp
            xmu=xmu+xma(ij)*pp
            ymu=ymu+yma(ij)*pp
            zmu=zmu+zma(ij)*pp
            if(.not.velcorr.and.printexciton) then ! collect contributions to charges on each atom
               l=(io-1)*nva+(iv-noa)
               q2(1:ncent)=q2(1:ncent)+qiaa(1:ncent,l)*real(uu)
            endif
         enddo
c beta P-CSF
         do j=nexa+1,nex
            io=iconfb(j-nexa,1)
            iv=iconfb(j-nexa,2)
            ij=lin(io,iv)
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax)then
               umax=uu*pp
               imax=io
               jmax=iv
               kmem=j
               ispin=2
            endif
            xlu=xlu+xlb(ij)*uu
            ylu=ylu+ylb(ij)*uu
            zlu=zlu+zlb(ij)*uu
            xvu=xvu+xvb(ij)*pp
            yvu=yvu+yvb(ij)*pp
            zvu=zvu+zvb(ij)*pp
            xmu=xmu+xmb(ij)*pp
            ymu=ymu+ymb(ij)*pp
            zmu=zmu+zmb(ij)*pp
            if(.not.velcorr.and.printexciton) then ! collect contributions to charges on each atom
               l=(io-1)*nvb+(iv-nob)
               q2(1:ncent)=q2(1:ncent)+qiab(1:ncent,l)*real(uu)
            endif
         enddo
c polarizability         
         alp_real(1)=alp_real(1)+xlu*xlu*ef 
         alp_real(2)=alp_real(2)+xlu*ylu*ef
         alp_real(3)=alp_real(3)+ylu*ylu*ef
         alp_real(4)=alp_real(4)+xlu*zlu*ef
         alp_real(5)=alp_real(5)+ylu*zlu*ef
         alp_real(6)=alp_real(6)+zlu*zlu*ef

c end of loop over configurations within an excited state

         if(.not.velcorr.and.printexciton) then
           ! this is the shift of the magnetic moment
           xms=(yvu*coc(3)-zvu*coc(2))
           yms=(zvu*coc(1)-xvu*coc(3))
           zms=(xvu*coc(2)-yvu*coc(1))
           write(27,'(i5,a,x,F10.4,x,a2)') i,':',de*27.21139,'eV'
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'l:',xlu,ylu,zlu
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'p:',xvu,yvu,zvu
           write(27,'(a2,x,f12.8,x,f12.8,x,f12.8)')'m:',xmu-xms,ymu-yms,
     .                                                  zmu-zms
           do k=1,ncent
             write(27,'(f12.8)') q2(k)
           enddo
         endif

         xp=xlu*xlu+ylu*ylu+zlu*zlu       
         fl=de*p23*xp 
         xp=xvu*xvu+yvu*yvu+zvu*zvu       
         fv=ef*p23*xp              
         xp=xlu*xmu+ylu*ymu+zlu*zmu       
         rl=-235.7220d0*xp          
! recalculated value for Rot with recent CODATA values and changed it (used to be 235.730) 
         xp=xvu*xmu+yvu*ymu+zvu*zmu       
         rv=-235.7220d0*xp*ef
         rvp(i)=rv
c sum rule
         sumf=sumf+fl
! memorize for printout
         umerk(1,i)=de
         umerk(2,i)=fl
         umerk(3,i)=fv
         umerk(4,i)=rl
         umerk(5,i)=rv
         umerk(6,i)=dble(uci(kmem,i))
         if(rpachk) umerk(6,i)=dble(uci(kmem,i))**2-dble(hci(kmem,i))**2
         umerk(7,i)=imax
         umerk(8,i)=jmax
         umerk(15,i)=ispin

c second largest transition
         umax=-1
         lmem=1
         ispin=1
         do j=1,nexa
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfa(j,1)
               jmax=iconfa(j,2)
               lmem=j
               ispin=1
            endif
         enddo
         do j=nexa+1,nex
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfb(j-nexa,1)
               jmax=iconfb(j-nexa,2)
               lmem=j
               ispin=2
            endif
         enddo
         umerk(9,i) =dble(uci(lmem,i))
         if(rpachk) umerk(9,i)=dble(uci(lmem,i))**2-dble(hci(lmem,i))**2
         umerk(10,i)=imax
         umerk(11,i)=jmax
         umerk(16,i)=ispin

c third largest transition
         umax=-1
         ispin=1
         do j=1,nexa
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.lmem.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfa(j,1)
               jmax=iconfa(j,2)
               jmem=j
               ispin=1
            endif
         enddo
         do j=nexa+1,nex
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(rpachk)then
              uu=uu+dble(hci(j,i))
              pp=pp-dble(hci(j,i))
            endif
            if((uu*pp).gt.umax.and.j.ne.lmem.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfb(j-nexa,1)
               jmax=iconfb(j-nexa,2)
               jmem=j
               ispin=2
            endif
         enddo
         umerk(12,i)=dble(uci(jmem,i))
         if(rpachk)umerk(12,i)=dble(uci(jmem,i))**2-dble(hci(jmem,i))**2
         umerk(13,i)=imax
         umerk(14,i)=jmax
         umerk(17,i)=ispin

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
c end loop over excited states

      if(.not.rpachk) deallocate(hci)



 11   format(i5,f9.3,f8.1,2f11.4)
 12   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3,
     .       3x,'Kshft:',f8.3)
 13   format(                   41x,f8.3,13x,f8.3)
 14   format(3x,F6.2,'(',i4,'a','->',i4,'a',')')
 15   format(3x,F6.2,'(',i4,'b','->',i4,'b',')')


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

! perform TDA eigenvector correction
      else
         ! if velocity correction is on, get improved rotatory strengths
         if (printexciton) then
           call apbtrafoexc_uks(nci,nexa,nexb,nroot,uci,eci,xla,yla,zla,
     . xva,yva,zva,xma,yma,zma,xlb,ylb,zlb,xvb,yvb,zvb,xmb,ymb,zmb,xmolw
     . ,noa,nva,nob,nvb,coc,ncent,qiaa,qiab,maxconfa,maxconfb,iconfa,
     . iconfb,rvp)
         else 
           call apbtrafo_uks(nci,nexa,nexb,nroot,uci,eci,xla,yla,zla,xva
     . ,yva,zva,xma,yma,zma,xlb,ylb,zlb,xvb,yvb,zvb,xmb,ymb,zmb,xmolw,
     . maxconfa,maxconfb,iconfa,iconfb,rvp) 
         endif
         do i=1,nroot
            umerk(5,i)=rvp(i)
         enddo
      endif

      if(printexciton) then
        deallocate(qiaa,qiab)
        close(27)
      endif

c output
      write(*,'(/,A,/)')
     .'excitation energy, transition moments and TDA amplitudes'
      if(velcorr) then
        write(*,*) 'state    eV      nm       fL        Rv(corr)'
      else
        write(*,*) 'state    eV      nm       fL         Rv'
      endif
      do i=1,nroot
         ec=umerk(1,i)
         write(*,11,advance='no') i,ec*27.21139,
     .   1.d+7/(ec*2.19474625d+5),umerk(2,i),umerk(5,i)
         if(umerk(15,i).eq.1.) then
          write(*,14,advance='no') umerk(6,i),idint(umerk(7,i)),
     .           idint(umerk(8,i))
         else
          write(*,15,advance='no') umerk(6,i),idint(umerk(7,i)),
     .           idint(umerk(8,i))
         endif
         if(umerk(16,i).eq.1.) then
          write(*,14,advance='no') umerk(9,i),idint(umerk(10,i)),
     .           idint(umerk(11,i))
         else
          write(*,15,advance='no') umerk(9,i),idint(umerk(10,i)),
     .           idint(umerk(11,i))
         endif
         if(umerk(17,i).eq.1.) then
          write(*,14,advance='no') umerk(12,i),idint(umerk(13,i)),
     .           idint(umerk(14,i))
         else
          write(*,15,advance='no') umerk(12,i),idint(umerk(13,i)),
     .           idint(umerk(14,i))
         endif
         write(*,'(a)')' '
      enddo


      alp_real=alp_real*2.0d0
      call prmat(6,alp_real,3,0,'alpha tensor')
      write(*,*) 'trace alpha_L[0] / au:',
     .(alp_real(1)+alp_real(3)+alp_real(6))/3.d0
      write(*,*) 'sum rule f_L        ',sumf

      call sosor(nroot,xmolw,eci,rvp)
      deallocate(rvp)
*****************************
c Lin. Response func. ESA   *
*****************************    
      if(ESA)then
      ! Unrelaxed state-to-state transition dipole moments
      call cpu_time(start_time)
      if(rpachk) then
      call ulresp_ESA(nexa,nexb,nex,iconfa,iconfb,maxconfa,maxconfb,
     .              xla,yla,zla,xlb,ylb,zlb,mocia,mocib,
     .              noa,nob,nva,nvb,eci,uci,hci,nroot,xmolw,thr,ak)
      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      else
      call ulresp_ESA_tda(nexa,nexb,nex,iconfa,iconfb,maxconfa,
     .              maxconfb,xla,yla,zla,xlb,ylb,zlb,mocia,mocib,
     .              noa,nob,nva,nvb,eci,uci,nroot,xmolw,thr,ak)
      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      endif
      write(*,*)
      !CALL EXIT([STATUS]) 
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optional: print eigenvectors in TM format
! Write $eigenpairs to scratch file
      if (eigvec) then
        ihilf=2 ! use as switch for output 

        open(unit=39,file='TmPvEcInFo',status='replace')
        ! print to temporary file
        write(39,*) 2
        write(39,*) nvec,nmo,jhomo
        write(39,*) nexa,nexb,jhomoa,jhomob
        do i=1,nmo/2
          write(39,*) vecchka(i)
        enddo
        do i=1,nmo/2
          write(39,*) vecchkb(i)
        enddo
        do i=1,nexa
           write(39,*) iconfa(i,1),iconfa(i,2)
        enddo
        do i=1,nexb
           write(39,*) iconfb(i,1),iconfb(i,2)
        enddo
        close(39)
        deallocate(vecchka,vecchkb)
!        call printeigenvec(rpachk,ihilf,nci,nroot,uci,hci,eci,nmo, 
!     .                    nvec,jhomo,maxconf1,vecchk,iconf)
        if(rpachk)then
          call printvecrpa(nci,nroot,uci,hci,eci)
          deallocate(hci)
        else
          call printvectda(rpachk,nci,nroot,uci,eci)
        endif
      endif
      
      if(nto)then
      if(rpachk)then
      call Uprint_nto_rpa(uci,hci,cca,mocia,nexa,nci,nroot,nao,iconfa,
     .                   maxconfa,noa,nva,ccb,mocib,nexb,iconfb,
     .                   maxconfb,nob,nvb)
      else
      
      call Uprint_nto(uci,cca,mocia,nexa,nci,nroot,nao,iconfa,
     .                   maxconfa,noa,nva,ccb,mocib,nexb,iconfb,
     .                   maxconfb,nob,nvb)
      endif
      endif      

      deallocate(uci)

      if(rpachk) then
      write(*,*) 'sTD-DFT done.'
      else
      write(*,*) 'sTDA done.'
      endif


c used for fitting      
      inquire(file='.REF',exist=ex)
      if(.not.ex) return
      open(unit=27,file='.REF')
      open(unit=28,file='.OUT')
      read(27,*) i,k

c     do j=1,i
c        read(27,*) hilf
c        if(hilf.gt.0.01)
c    .   write(28,'(4f12.6)') 
c    .   (umerk(1,j)*27.211-hilf),hilf,umerk(1,j)*27.211,umerk(2,j)
c     enddo

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
      elseif(k.eq.4)then
c take only very very bright ones      
      read(27,*) hilf
      do j=1,nroot
         if(umerk(2,j).gt.2.0)then
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

      close(27)
      close(28)

      return
      end


**********************************************************************

      subroutine sort_vec(nmo,nao,c,epsi)
      implicit none
      integer, intent ( in ) :: nmo,nao
      real*8, intent ( inout ) :: c(*),epsi(*)
      integer ii,i,j,k,m
      real*8 pp,hilf

      do 1141 ii = 2, nmo
         i = ii - 1
         k = i
         pp= epsi(i)
         do 1121 j = ii, nmo
            if (epsi(j) .ge. pp) go to 1121
            k = j
            pp=epsi(j)
 1121    continue
         if (k .eq. i) go to 1141
         epsi(k) = epsi(i)
         epsi(i) = pp
         do m=1,nao
          hilf=c(m+(i-1)*nao)
          c(m+(i-1)*nao)=c(m+(k-1)*nao)
          c(m+(k-1)*nao)=hilf
         enddo
 1141 continue

      return
      end


***********************************************************************
* select important CSFs beyond energy threshold (UKS case)
***********************************************************************
      subroutine ptselect_uks(nexa,nexb,ncent,noa,nva,nob,nvb,nexpta,
     .                  nexptb,mxcnfa,mxcnfb,iconfa,iconfb,kconfa,
     .                  kconfb,dak,dax,eda,edb,edpta,edptb,piaa,qiaa,
     .                  piab,qiab,pija,qaba,pijb,qabb,thrp,newa,newb)
      use omp_lib
      implicit none
      integer, intent(in) :: nexa,nexb,ncent,nva,nvb,nexpta,nexptb
      integer, intent(in) :: noa,nob,mxcnfa,mxcnfb
      integer, intent(in) :: kconfa(mxcnfa,2),kconfb(mxcnfb,2)
      integer, intent (inout) :: iconfa(mxcnfa,2),iconfb(mxcnfb,2)
      integer, intent (out) :: newa,newb
      real*4, intent(in)  :: qiaa(ncent,mxcnfa),piaa(ncent,mxcnfa)
      real*4, intent(in)  :: qiab(ncent,mxcnfb),piab(ncent,mxcnfb)
      real*4, intent(in)  :: qaba(ncent,nva*(nva+1)/2)
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*4, intent(in)  :: pijb(ncent,nob*(nob+1)/2)
      real*8, intent(in)  :: dak,dax,thrp
      real*8, intent(in)  :: edpta(mxcnfa),edptb(mxcnfb)
      real*8, intent(inout) :: eda(mxcnfa),edb(mxcnfb)
      real*4, allocatable :: qk(:),qj(:)
      real*8, allocatable :: pta(:),pt2a(:),ptb(:),pt2b(:)
      integer i,j,k,l,io,iv,jo,jv,ierr,lin,iwrk,jwrk,iiv,jjv,nex
      real*4 ek,ej,sdot,tmpi,tmpj
      real*8 de,pert,amat
      logical, allocatable :: incl_conf(:)

      allocate(qk(ncent), pta(nexa),pt2a(nexa),ptb(nexb),pt2b(nexb), 
     .         qj(ncent),incl_conf(nexpta),stat=ierr)
      if(ierr.ne.0)stop 'allocation for PT intermediates crashed'
      incl_conf=.false.
      qk=0.0
      qj=0.0
      nex=nexa+nexb
      newa = 0
      newb = 0
      ptb = 0.0d0
      pta= 0.0d0
      pt2b = 0.0d0
      pt2a = 0.0d0


! loop over secondary/neglected CSF, this is done omp-parallel
!$omp parallel private(i,k,io,iv,iiv,iwrk,qk,de,j,jo,jv,jjv,jwrk,ek,ej,
!$omp&                 l,qj,pta,ptb,amat,pert) reduction (+:pt2a,pt2b) 
!$omp do 
c outer loop over alpha S-CSF
      do k=1,nexpta
         io=kconfa(k,1)
         iv=kconfa(k,2)
         iiv=iv-noa
         iwrk=(io-1)*nva + iiv
         qk(1:ncent)=piaa(1:ncent,iwrk)
         de=edpta(k)
c loop over alpha P-CSF         
         do j=1,nexa
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
c coupling (K-J)(ia(alpha),jb(alpha))
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
c coupling
            qj(1:ncent)=pija(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qaba(1,lin(iiv,jjv)),1)
            pta(j)=0.0d0
            amat=ek-ej
c PT2 - e(2) = -((K-J)(ia,jb)**2)/(A(jb,jb)-A(ia,ia))
            pta(j)=amat**2/(de-eda(j)+1.d-10)
         enddo
c loop over beta P-CSF         
         do j=1,nexb
            jo=iconfb(j,1)
            jv=iconfb(j,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
c coupling
            ptb(j)=0.0d0
c coupling (K)(ia(alpha),jb(beta))
            ek=sdot(ncent,qk,1,qiab(1,jwrk),1)
            amat=ek
c PT2 - e(2) = -((K)(ia,jb)**2)/(A(jb,jb)-A(ia,ia))
            ptb(j)=amat**2/(de-edb(j)+1.d-10)
         enddo
         pert=sum(pta(1:nexa)) + sum(ptb(1:nexb))
c if sum > threshold include the alpha S-CSF in the alpha P-CSF
         if(pert.gt.thrp)then
            incl_conf(k)=.true.
         else
c else accumulate in the E(PT2) contribution to the alpha and beta P-CSF
            do i=1,nexa
               pt2a(i)=pt2a(i) + pta(i)              
            enddo
            do l=1,nexb
               pt2b(l)=pt2b(l) + ptb(l)
            enddo
         endif
      enddo
!$omp end do
!$omp end parallel

      do i=1,nexpta
         if(.not.incl_conf(i)) cycle
            io=kconfa(i,1)
            iv=kconfa(i,2)
            de=edpta(i)
            newa=newa+1
            iconfa(nexa+newa,1)=io
            iconfa(nexa+newa,2)=iv
            eda  (nexa+newa  )=de
      enddo

! reset incl_conf for beta space
      deallocate(incl_conf)
      allocate(incl_conf(nexptb),stat=ierr)
      if(ierr.ne.0)stop 'allocation for incl_conf_b crashed'
      incl_conf=.false.

!$omp parallel private(i,k,io,iv,iiv,iwrk,qk,de,j,jo,jv,jjv,jwrk,ek,ej,
!$omp&                l,qj,pta,ptb,amat,pert) reduction (+:pt2a,pt2b) 
!$omp do 
c outer loop over beta S-CSF
      do k=1,nexptb
         io=kconfb(k,1)
         iv=kconfb(k,2)
         iiv=iv-nob
         iwrk=(io-1)*nvb + iiv
         qk(1:ncent)=piab(1:ncent,iwrk)
         de=edptb(k)
c loop over beta P-CSF         
         do j=1,nexb
            jo=iconfb(j,1)
            jv=iconfb(j,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
c coupling (K-J)(ia(beta),jb(beta))
            ek=sdot(ncent,qk,1,qiab(1,jwrk),1)
c J coupling
            qj(1:ncent)=pijb(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1)
            ptb(j)=0.0d0
            amat=ek-ej
c PT2 - e(2) = -((K-J)(ia,jb)**2)/(A(jb,jb)-A(ia,ia))
            ptb(j)=amat**2/(de-edb(j)+1.d-10)
         enddo
c loop over alpha P-CSF         
         do j=1,nexa
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
c coupling
            pta(j)=0.0d0
c coupling (K)(ia(alpha),jb(beta))
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
            amat=ek
c PT2 - e(2) = -((K)(ia,jb)**2)/(A(jb,jb)-A(ia,ia))
            pta(j)=amat**2/(de-eda(j)+1.d-10)
         enddo
         pert=sum(pta(1:nexa)) + sum(ptb(1:nexb))
c if sum > threshold include the beta S-CSF in the beta P-CSF
         if(pert.gt.thrp)then
            incl_conf(k)=.true.
         else
c else accumulate in the E(PT2) contribution to the alpha and beta P-CSF 
            do i=1,nexa
               pt2a(i)=pt2a(i) + pta(i)
            enddo
            do l=1,nexb
               pt2b(l)=pt2b(l) + ptb(l)
            enddo
         endif
      enddo
!$omp end do
!$omp end parallel

      do i=1,nexptb
         if(.not.incl_conf(i)) cycle
            io=kconfb(i,1)
            iv=kconfb(i,2)
            de=edptb(i)
            newb=newb+1
            iconfb(nexb+newb,1)=io
            iconfb(nexb+newb,2)=iv
            edb  (nexb+newb  )=de
      enddo

      deallocate(incl_conf)

      pert=0.0d0
      do i=1,nexa
         pert=pert+pt2a(i)
         eda(i)=eda(i)-pt2a(i)
      enddo
      do i=1,nexb
         pert=pert+pt2b(i)
         edb(i)=edb(i)-pt2b(i)
      enddo

      write(*,'('' average/max(a)/max(b) PT2 energy lowering (eV):'',
     .3F10.3)')
     .             27.21139*pert/float(nex),maxval(pt2a)*27.21139,
     .maxval(pt2b)*27.21139


      deallocate(qk,qj,pta,ptb,pt2a,pt2b)
      return

      end subroutine ptselect_uks
***********************************************************************

***********************************************************************
* set up A+B and A-B (packed form) in UKS case 
***********************************************************************
      subroutine urpamat(nex,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,
     .             mxcnfb,iconfa,iconfb,dax,eda,edb,piaa,qiaa,piab,
     .             qiab,pija,qaba,pijb,qabb,ambsqr)
      use omp_lib
      implicit none
      integer, intent(in) :: nex,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa
      integer, intent(in) :: mxcnfb,iconfa(mxcnfa,2),iconfb(mxcnfb,2)
      real*4, intent(in)  :: qiaa(ncent,mxcnfa),piaa(ncent,mxcnfa)
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qaba(ncent,nva*(nva+1)/2)
      real*4, intent(in)  :: qiab(ncent,mxcnfb),piab(ncent,mxcnfb)
      real*4, intent(in)  :: pijb(ncent,nob*(nob+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*4, intent(out) :: ambsqr(nex*(nex+1)/2)
      real*8, intent(in)  :: dax,eda(mxcnfa),edb(mxcnfb)
      real*4, allocatable :: qj(:),qk(:),apb(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,k,iiv,jjv,iwrk,jwrk
      real*4 ek,ej,sdot,ax,de
      allocate(qk(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qk crashed'
      ax=real(dax)

      open(unit=52,file='apbmat',form='unformatted',status='replace')
      ambsqr=0.0e0

      if(abs(dax).lt.1.0d-6) then 
c calculate A+B and (A-B)^0.5
        allocate(apb(nex*(nex+1)/2),stat=ierr)
        if(ierr.ne.0)stop 'allocation for A+B crashed'
        apb=0.0e0
        ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek)
!$omp do
c alpha-alpha block
        do i = 1,nexa
           io=iconfa(i,1)
           iv=iconfa(i,2)
           iiv=iv-noa
           iwrk=(io-1)*nva + iiv
           qk(1:ncent)=piaa(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconfa(j,1)
              jv=iconfa(j,2)
              jjv=jv-noa
              jwrk=(jo-1)*nva+jjv
              ek=sdot(ncent,qk,1,qiaa(1,jwrk),1) ! ek = (ia|jb)
              apb(ij)=2.0*ek
              ambsqr(ij)=0.0
           enddo
           ij=lin(i,i)
           ek=sdot(ncent,qk,1,qiaa(1,iwrk),1)
           ambsqr(ij)=sqrt(real(eda(i))-ek) ! diagonal element of (A-B)^0.5
           apb(ij)=real(eda(i))+ek ! diagonal element of A+B
        enddo
!$omp end do
!$omp end parallel
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek)
!$omp do
c beta block
        do i = nexa+1,nex
           io=iconfb(i-nexa,1)
           iv=iconfb(i-nexa,2)
           iiv=iv-nob
           iwrk=(io-1)*nvb + iiv
           qk(1:ncent)=piab(1:ncent,iwrk)
!... with alpha
           do j = 1,nexa
              ij=lin(i,j)
              jo=iconfa(j,1)
              jv=iconfa(j,2)
              jjv=jv-noa
              jwrk=(jo-1)*nva+jjv
              ek=sdot(ncent,qk,1,qiaa(1,jwrk),1) ! ek = (ia|jb)
              apb(ij)=2.0*ek
           enddo   
!... with beta
           do j=nexa+1,i-1
              ij=lin(i,j)
              jo=iconfb(j-nexa,1)
              jv=iconfb(j-nexa,2)
              jjv=jv-nob
              jwrk=(jo-1)*nvb+jjv
              ek=sdot(ncent,qk,1,qiab(1,jwrk),1) ! ek = (ia|jb)
              apb(ij)=2.0*ek
           enddo
           ij=lin(i,i)
           ek=sdot(ncent,qk,1,qiab(1,iwrk),1)
           ambsqr(ij)=sqrt(real(edb(i-nexa))-ek) ! diagonal element of (A-B)^0.5
           apb(ij)=real(edb(i-nexa))+ek ! diagonal element of A+B
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
        allocate(qj(ncent),stat=ierr)
        if(ierr.ne.0)stop 'allocation failed for qk vector'
c calculate A+B and A-B
c alpha-alpha block
        ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek,qj,ej)
!$omp do
        do i = 1,nexa
           io=iconfa(i,1)
           iv=iconfa(i,2)
           iiv=iv-noa
           iwrk=(io-1)*nva + iiv
           qk(1:ncent)=piaa(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconfa(j,1)
              jv=iconfa(j,2)
              jjv=jv-noa
              jwrk=(jo-1)*nva+jjv
              ek=sdot(ncent,qk,1,qiaa(1,jwrk),1) ! ek = (ia|jb)
              qj(1:ncent)=pija(1:ncent,lin(io,jo))
              ej=sdot(ncent,qj,1,qaba(1,lin(iiv,jjv)),1) !  ej = (ij|ab)
              ambsqr(ij)=2.0*ek-ej
              jwrk=(io-1)*nva+jjv
              qj(1:ncent)=piaa(1:ncent,jwrk)
              jwrk=(jo-1)*nva+iiv
              ek=sdot(ncent,qj,1,qiaa(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ! A+B part
              ambsqr(ij)=ambsqr(ij)-ax*ek 
              ! we first use apb as A-B 
              apb(ij)=ax*ek-ej 
           enddo
           ij=lin(i,j)
           ek=sdot(ncent,qk,1,qiaa(1,iwrk),1)
           apb(ij)=real(eda(i))-ek+ax*ek! diagonal element of A-B
           ambsqr(ij)=real(eda(i))-ax*ek+ek ! diagonal element of A+B
        enddo
!$omp end do
!$omp end parallel

! beta...
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek,qj,ej)
!$omp do
        do i = nexa+1,nex
           io=iconfb(i-nexa,1)
           iv=iconfb(i-nexa,2)
           iiv=iv-nob
           iwrk=(io-1)*nvb + iiv
           qk(1:ncent)=piab(1:ncent,iwrk)
! ...alpha block
           do j = 1,nexa
              ij=lin(i,j)
              jo=iconfa(j,1)
              jv=iconfa(j,2)
              jjv=jv-noa
              jwrk=(jo-1)*nva+jjv
              ek=sdot(ncent,qk,1,qiaa(1,jwrk),1) ! ek = (ia|jb)
              ambsqr(ij)=2.0*ek
           enddo
! ...beta block
           do j = nexa+1,i-1
              ij=ij+1
              jo=iconfb(j-nexa,1)
              jv=iconfb(j-nexa,2)
              jjv=jv-nob
              jwrk=(jo-1)*nvb+jjv
              ek=sdot(ncent,qk,1,qiab(1,jwrk),1) ! ek = (ia|jb)
              qj(1:ncent)=pijb(1:ncent,lin(io,jo))
              ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1) !  ej = (ij|ab)
              ambsqr(ij)=2.0*ek-ej
              jwrk=(io-1)*nvb+jjv
              qj(1:ncent)=piab(1:ncent,jwrk)
              jwrk=(jo-1)*nvb+iiv
              ek=sdot(ncent,qj,1,qiab(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ambsqr(ij)=ambsqr(ij)-ax*ek ! scaled by ax
           enddo
           ij=lin(i,j)
           ek=sdot(ncent,qk,1,qiab(1,iwrk),1)
           ambsqr(ij)=real(edb(i-nexa))+ek-ax*ek ! diagonal element of A+0.5*B  (beta part)
        enddo
!$omp end do
!$omp end parallel

        write(52)ambsqr
        
c this the time determining step, since A-B needs to be diagonalized in full nci space
        write(*,*) ' calculating (A-B)^0.5 ...'
        write(*,'('' estimated time (min) '',f8.2)')
     . (float(nexb)**2*float(nexb)+float(nexa)**2*float(nexa))/4.d+8/60.

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

c beta-beta block
        ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek,qj,ej)
!$omp do
        do i = 1,nexb
           io=iconfb(i,1)
           iv=iconfb(i,2)
           iiv=iv-nob
           iwrk=(io-1)*nvb + iiv
           qk(1:ncent)=piab(1:ncent,iwrk)
           do j = 1,i-1
              ij=lin(i,j)
              jo=iconfb(j,1)
              jv=iconfb(j,2)
              jjv=jv-nob
              jwrk=(jo-1)*nvb + jjv
              qj(1:ncent)=pijb(1:ncent,lin(io,jo))
              ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1) !  ej = (ij|ab)
              jwrk=(io-1)*nvb+jjv
              qj(1:ncent)=piab(1:ncent,jwrk)
              jwrk=(jo-1)*nvb+iiv
              ek=sdot(ncent,qj,1,qiab(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ! we first use apb as A-B 
              apb(ij)=ax*ek-ej
           enddo
           ij=lin(i,i)
           ek=sdot(ncent,qk,1,qiab(1,iwrk),1)
           apb(ij)=real(edb(i))-ek+ax*ek ! diagonal element of A-B
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
        deallocate(apb,qj)
      endif ! GGA/hybrid case 
      close(52)
      deallocate(qk)

      return

      end subroutine urpamat
***********************************************************************

***********************************************************************
* set up A matrix in UKS case 
***********************************************************************
       subroutine utdamat(nci,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa,
     .                    mxcnfb,iconfa,iconfb,dax,eda,edb,piaa,qiaa,
     .                    piab,qiab,pija,qaba,pijb,qabb,hci)
       use omp_lib
      implicit none
      integer, intent(in) :: nci,nexa,nexb,ncent,noa,nva,nob,nvb,mxcnfa
      integer, intent(in) :: mxcnfb,iconfa(mxcnfa,2),iconfb(mxcnfb,2)
      real*4, intent(in)  :: qiaa(ncent,mxcnfa),piaa(ncent,mxcnfa)
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qaba(ncent,nva*(nva+1)/2)
      real*4, intent(in)  :: qiab(ncent,mxcnfb),piab(ncent,mxcnfb)
      real*4, intent(in)  :: pijb(ncent,nob*(nob+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*4, intent(out)  :: hci(nci,nci)
      real*8, intent(in)  :: dax,eda(mxcnfa),edb(mxcnfb)
      real*4, allocatable :: qj(:),qk(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,k,iiv,jjv,iwrk,jwrk
      real*4 ek,ej,sdot,ax,de
      allocate(qk(ncent),qj(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj crashed'
! calculate CIS matrix A
      hci=0.0e0
      ax=real(dax)
!$omp parallel private(i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
c alpha-alpha block
      do i = 1,nexa
         io=iconfa(i,1)
         iv=iconfa(i,2)
         iiv=iv-noa
         iwrk=(io-1)*nva + iiv
         qk(1:ncent)=piaa(1:ncent,iwrk)
         do j=1,i-1
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
            qj(1:ncent)=pija(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qaba(1,lin(iiv,jjv)),1)
            hci(j,i)=ek-ej
            hci(i,j)=hci(j,i)
         enddo
         hci(i,i)=real(eda(i))                         
      enddo
!$omp end do
!$omp end parallel 

c beta blocks
!$omp parallel private(i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
      do i = nexa+1,nci
         io=iconfb(i-nexa,1)
         iv=iconfb(i-nexa,2)
         iiv=iv-nob
         iwrk=(io-1)*nvb + iiv
         qk(1:ncent)=piab(1:ncent,iwrk)
! alpha-beta
         do j = 1,nexa
            jo=iconfa(j,1)
            jv=iconfa(j,2)
            jjv=jv-noa
            jwrk=(jo-1)*nva + jjv
            ek=sdot(ncent,qk,1,qiaa(1,jwrk),1)
            hci(j,i)=ek
            hci(i,j)=ek
         enddo
! beta-beta
         do j = nexa+1,i-1
            jo=iconfb(j-nexa,1)
            jv=iconfb(j-nexa,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
            ek=sdot(ncent,qk,1,qiab(1,jwrk),1)
            qj(1:ncent)=pijb(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1)
            hci(j,i)=ek-ej
            hci(i,j)=hci(j,i)
         enddo
         hci(i,i)=real(edb(i-nexa))
      enddo
!$omp end do
!$omp end parallel

      deallocate(qk,qj)

      return

      end subroutine utdamat
***********************************************************************

      subroutine somo_shift(nsomo,io,iv,de,ek)
      use kshiftcommon
      implicit none
      real*8, intent(inout) :: de ! A(ia,ia)
      integer,intent (in) :: nsomo,io,iv    
      real*8, intent (in) :: ek ! K(ia,ia)
      real*8  shft,wau,sau,ekia,dmp
      logical l1,l2

      l1=io.gt.nsomo
      l2=iv.le.nsomo
     
      if(l1.or.l2)then
      ! shiftmax and width are given eV, convert to a.u.       
      sau=shftmax_somo*0.036749320d0
      wau=shftwidth*0.036749320d0*5.0d0
      shft = 1.0d0 + (ek/wau)**shftsteep
      de = de + sau/shft ! shift A(ia,ia)
      endif

      end
