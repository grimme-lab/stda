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

! written by Marc de Wegifosse 2018-2019

      SUBROUTINE sfstda(ncent,nmo,nao,xyz,c,eps,occ,csp,iaoat,thr,thrp,
     .            ax,alphak,betaj,fthr,nvec)
      use commonlogicals
      use omp_lib
      IMPLICIT NONE
c input:
      integer ncent,nmo,nao
      integer iaoat(*)
      real*8  thr,thrp,ax,othr,vthr
      real*8  c(*),eps(*),occ(*),xyz(4,*)
      integer csp(*)
      logical ggavec,ex

c local variables:

      integer :: STATUS

c l=dipole lengths, v=dipole velocity (d/dxyz), m=angular momentum
      real*8, allocatable ::xla(:),yla(:),zla(:)
!       real*8, allocatable ::xva(:),yva(:),zva(:)
!       real*8, allocatable ::xma(:),yma(:),zma(:)
      real*8, allocatable ::xlb(:),ylb(:),zlb(:)
!       real*8, allocatable ::xvb(:),yvb(:),zvb(:)
!       real*8, allocatable ::xmb(:),ymb(:),zmb(:)
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
!      real*4, allocatable ::gamk(:,:)
!      real*4, allocatable ::qiaa(:,:),qiab(:,:),piaa(:,:),piab(:,:)
      real*4, allocatable ::pija(:,:),qija(:,:)
      real*4, allocatable ::qabb(:,:)
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
      integer maxconfb
      real*8, allocatable ::edb(:)
      real*8, allocatable ::edptb(:)
      integer, allocatable :: iconfb(:,:)
      integer, allocatable :: kconfb(:,:)

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
!      real*8 alp_real(6),sumf
      real*4 :: start_time, end_time, stda_time
! variables for vector printout
      integer, allocatable :: vecchka(:),vecchkb(:)
      integer nvec,jhomo,jhomoa,jhomob
c atomic Hubbard parameters
      real*8 eta(94)
c atomic masses
      common /amass/ ams(107)
      real*8 ams
      character*79 dummy
c For the computation of <S**2>
      real*8, allocatable :: i_alpha_j_beta(:,:)
      real*8, allocatable :: i_alpha_a_beta(:,:)
      real*8 :: sum_ij,S2_UHF
      real*8, allocatable :: Xia(:),Xiaib(:),Xiaka(:)
      real*8, allocatable :: S2(:)
      integer kk

c just a printout
      call header('s T D A',0)
      call cpu_time(start_time)
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

      allocate(xla(mocia*(mocia+1)/2),yla(mocia*(mocia+1)/2),
     .         zla(mocia*(mocia+1)/2),
     .         xlb(mocib*(mocib+1)/2),ylb(mocib*(mocib+1)/2),
     .         zlb(mocib*(mocib+1)/2),help(nao*(nao+1)/2),
     .         clowa(nao*mocia),clowb(nao*mocib),
     .         scr(nao*nao),dum(nao*nao),x(nao,nao),
     .         gamj(ncent,ncent),
     .         cca(nao*mocia),epsia(mocia),
     .         ccb(nao*mocib),epsib(mocib)
     .        )

      write(*,*) 'Active space MOs in TDA : ',moci,mocia,mocib

      mocia=0
      mocib=0
      moci =0

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

      nv  = moci - no
      nva = mocia - noa
      nvb = mocib - nob

      write(*,*) 'Occupied active MOs in TDA: ', no,noa,nob
      if(noa<=nob)then
      write(*,*)'warning number of alpha < beta'
      write(*,*)'only with alpha > beta'
      call exit(status)
      endif
      write(*,*) 'Virtual active MOs in TDA: ', nv,nva,nvb
      if((noa.eq.0.or.nva.eq.0).and.(nob.eq.0.or.nvb.eq.0)) then
        stop 'no CSF, increase energy threshold (-e option)'
      endif

      ! SPIN-FLIP i_alpha -> b_beta
      maxconfb=noa*nvb
      allocate(edb(maxconfb),edptb(maxconfb), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for ed and edpt'
      allocate(iconfb(maxconfb,2),kconfb(maxconfb,2), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for iconf and kconf'
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
!
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! c
! c            magnetic dipole
! c
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      open(unit=34,file='xmint',form='unformatted',status='old')
!       read(34) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,xma)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,xmb)
      close(34,status='delete')
!
      open(unit=35,file='ymint',form='unformatted',status='old')
!       read(35) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,yma)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,ymb)
      close(35,status='delete')
!
      open(unit=36,file='zmint',form='unformatted',status='old')
!       read(36) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,zma)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,zmb)
      close(36,status='delete')
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! c
! c            velocity dipole
! c
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      open(unit=37,file='xvint',form='unformatted',status='old')
!       read(37) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,xva)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,xvb)
      close(37,status='delete')
!
      open(unit=38,file='yvint',form='unformatted',status='old')
!       read(38) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,yva)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,yvb)
      close(38,status='delete')
!
      open(unit=39,file='zvint',form='unformatted',status='old')
!       read(39) help
!       call onetri(-1,help,dum,scr,cca,nao,mocia)
!       call shrink(mocia,dum,zva)
!       call onetri(-1,help,dum,scr,ccb,nao,mocib)
!       call shrink(mocib,dum,zvb)
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
      if(sf_s2)then
      ! for the computation of <S**2>
      allocate(i_alpha_j_beta(noa,nob))
      i_alpha_j_beta=0.0
      sum_ij=0.0
      do i=1,noa
        do j=1,nob
          do k=1,nao
          i_alpha_j_beta(i,j)=i_alpha_j_beta(i,j)+
     .         clowa(k+(i-1)*nao)*clowb(k+(j-1)*nao)
          enddo
          sum_ij=sum_ij+i_alpha_j_beta(i,j)**2.0
        enddo
      enddo
      S2_UHF=(noa-nob)/2.0*((noa-nob)/2.0+1.0)+nob-sum_ij
      write(*,*) 'Unrestricted ground state <S**2>=',S2_UHF
      allocate(i_alpha_a_beta(noa,nvb))
      i_alpha_a_beta=0.0
      do i=1,noa
        do j=1,nvb
          do k=1,nao
          i_alpha_a_beta(i,j)=i_alpha_a_beta(i,j)+
     .         clowa(k+(i-1)*nao)*clowb(k+(j+nob-1)*nao)
          enddo
        enddo
      enddo
      endif
      deallocate(scr,dum,help,x)


      call cpu_time(time)
      allocate(q1a(ncent),q1b(ncent),q2(ncent))
      q1a=0.0
      q1b=0.0
      q2=0.0
      allocate(qija(ncent,noa),
     .         qabb(ncent,nvb),stat=ierr)
      if(ierr.ne.0) stop 'error in diag. J charges allocation'
      qija=0.0
      qabb=0.0

      do i=1,noa
         call lo12pop(i,i,ncent,nao,iaoat,clowa,q2)
         q1a(1:ncent)=q1a(1:ncent)+q2(1:ncent)
         qija(1:ncent,i)=q2(1:ncent)
      enddo

!       write(*,'(/'' SCF atom population (using active MOs):'')')
!       write(*,'(10F7.3)') q1a(1:ncent)
!       write(*,'(10F7.3)') q1b(1:ncent)
!       write(*,'(10F7.3)') (q1a(1:ncent)+q1b(1:ncent))
!       write(*,*)
!       write(*,'('' # of alpha and beta electrons in TDA:'',2F8.3)')
!      .sum(q1a(1:ncent)),sum(q1b(1:ncent))
!       write(*,*)

      do i=nob+1,mocib
         j=i-nob
         call lo12pop(i,i,ncent,nao,iaoat,clowb,q2)
         qabb(1:ncent,j)=q2(1:ncent)
      enddo

c  UCIS scaling factor for K(ia|jb):
      ak = 1.0d0

c the global parameters of the method:
      beta1 = 0.3
      beta2 = 1.0


      if(betaj.lt.-99.0d0) then ! if no beta parameter was read in
      betaj=beta1+beta2*ax
      endif


      write(*,*)
      write(*,*) 'ax(DF)  : ',ax
      write(*,*) 's^K      : ',ak
      write(*,*) 'beta  (J): ',betaj
      write(*,*)

c set the hardness table
      call setrep(eta)
      write(*,*) 'hardness table read.'

      write(*,*) 'setting up gammas ...'
c compute gamma(j)
      xmolw=0
      do i = 1,ncent
         ii = idint(xyz(4,i))
c ams is the atomic mass (mol weight for output file)
         xmolw = xmolw+ams(ii)
         do j=1,i
            jj = idint(xyz(4,j))
            xj  = 0.50d0*(eta(ii)+eta(jj)) * (ax+ax*0.4) !!!!!!!!!!!!!!!!!!!!!!   SF this seems resonable
            rabx = dsqrt((xyz(1,i)-xyz(1,j))**2
     .                 +(xyz(2,i)-xyz(2,j))**2
     .                 +(xyz(3,i)-xyz(3,j))**2)
            gamj(j,i)=real(1./(rabx**betaj+1./xj**betaj)**(1.0d0/betaj))
            gamj(i,j) = gamj(j,i)
         enddo
      enddo


! ! compute intermediates q's refer to charges, p's to gam*q (i.e., contracted)
!       imem1=noa
!       imem2=nva
!       imem3=(noa*nva+nob*nvb)*2
!       imem1=noa*(noa+1)/2+nob*(nob+1)/2
!       imem2=nva*(nva+1)/2+nvb*(nvb+1)/2
!       imem3=imem3+imem1+imem2
!       hilf=dble(imem3)/1024.0**2
!       hilf=dble(4*ncent)*hilf
!       imem1=idint(hilf)
!       write(*,*)'memory needed for qa and qb data (Mb) ',imem1
!       write(*,*)'computing q(ij,n) ...'

! first alpha part

!     Coulomb type terms: calc qij*gam^J
      allocate(pija(ncent,noa), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii| intermediate'
      pija=0.0
      call ssymm('l','l',ncent,noa,1.0,gamj,ncent,qija,ncent,0.0,pija
     .           ,ncent)
      allocate(uci(nvb,noa), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for (ii|aa) matrix'
      uci=0.0
      ! now calc (ii|aa)^J SPIN-FLIP (alpha|beta)
      call sgemm('t','n',nvb,noa,ncent,1.0,qabb,ncent,pija,ncent,0.0,
     .           uci,nvb)
      deallocate(pija)

      nex=noa*nvb


c determine singles which are lower than thr using A(ia,ia)  !! SPIN-FLIP

      k=0
      j=0
      l=0
      nsomoa=noa+nuhf
      do io=1,noa
         do iv=nob+1,mocib
c compute A(ia,ia) = eps(a)-eps(i) - (ii|aa)    (alpha,beta)
            de=epsib(iv)-epsia(io)
            l=iv-nob
            ej=0.0d0
            ej=dble(uci(l,io))
            de=de-ej

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
      nex = nexb
      nexpt = nexptb

      write(*,*)
      write(*,*) '# CSF included by energy:',nex
      write(*,*)
      write(*,*) '# CSF considered in PT2:',nexpt

c errors and warning
      if(nex.lt.1) stop 'No CSF, increase energy threshold (-e option)'
!      if(nexa.eq.maxconf1.or.nexb.eq.maxconf1)
!     .   stop 'Primary CSF space exceeded. use -e option!'
!      if(nexpta.eq.maxconf2.or.nexptb.eq.maxconf2)
!     .   write(*,*)'CSF PT2 space exceeded. try -p option!'

c sort for E diag in each spin manifold

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
      write(*,*)'        eV'
      j=max(1,noa-10)
      do i=j,noa
         write(*,'(i4,F10.3,F8.1)') i,epsia(i)*27.21139
      enddo
      write(*,*)'Unoccupied:'
      j=min(mocia,noa+11)
      do i=noa+1,j
         write(*,'(i4,F10.3,F8.1)') i,epsia(i)*27.21139
      enddo

      write(*,*) ' '
      write(*,*)'Ordered frontier beta orbitals:'
      write(*,*)'        eV'
      j=max(1,nob-10)
      do i=j,nob
         write(*,'(i4,F10.3,F8.1)') i,epsib(i)*27.21139
      enddo
      write(*,*)'Unoccupied:'
      j=min(mocib,nob+11)
      do i=nob+1,j
         write(*,'(i4,F10.3,F8.1)') i,epsib(i)*27.21139
      enddo

      write(*,*)
      write(*,*)'            Lowest CSF states:'
      write(*,*)'         Alpha to Beta Excitations:'
      write(*,*)'      eV     nm      excitation i->a               eV'
      do i=1,min(nexb,25)
         io=iconfb(i,1)
         iv=iconfb(i,2)
         q1a(1:ncent)=qija(1:ncent,io)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1a,1,0.0,q2,1)
         jii=sdot(ncent,q1a,1,q2,1)
         k=iv-nob
         q1b(1:ncent)=qabb(1:ncent,k)
         ej=sdot(ncent,q1b,1,q2,1)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1b,1,0.0,q2,1)
         jaa=sdot(ncent,q1b,1,q2,1)
! de is now the Kia shift
         de=0
         loc=ej/sqrt(jii*jaa) ! locality
          write(*,16) i,27.211*edb(i),
     .   1.d+7/(edb(i)*2.19474625d+5),iconfb(i,1:2),
     .   27.211*(epsib(iv)-epsia(io)),27.211*ej,27.211*de,loc
      enddo

      deallocate(qija,qabb)
!!! now set up pij -alpha first
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

      deallocate(qija,gamj)
      if(XsTD.eqv..false.)then
      deallocate(clowa)
      endif
!!! now beta part of qab !!!

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
      if(XsTD.eqv..false.)then
      deallocate(clowb)
      endif
!!!!!PTselect, for spin-flip!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*)
      write(*,*)'selecting CSF ...'
c determine singles which interact most strongly with all basis
c CSF which have been selected by the  diagonal element A(ia,ia)
      deallocate(q1a,q1b,q2)
      call sf_ptselect_uks(nexb,ncent,noa,nva,nob,nvb,nexptb,
     .                  maxconfb,iconfb,kconfb,
     .                  ak,ax,edb,edptb,
     .                  pija,qabb,thrp,newb)

      deallocate(edptb,kconfb)

      new = newb


c nroot at this point is the number of primary CSF. The
c number of roots in the range 0-thr (as desired) is not
c known but will be determined by the diag routine.
       nroot = nex
      write(*,*) 'CSF included by PT:',new
      nexb = nexb + newb
      nex= nexb
!!!!!PTselect, for spin-flip!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nci = nex
      write(*,*) 'CSF in total:',nex

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
      if(XsTD)then
      call SF_Xstda_mat(nci,nexb,ncent,noa,nva,nob,nvb,maxconfb,iconfb,
     .ax,edb,hci,betaj,xyz,nao,mocia,mocib,clowa,clowb,epsia,epsib)
      else
      call sfstdamat(nci,nexb,ncent,noa,nob,nvb,
     .             maxconfb,iconfb,ax,edb,
     .             pija,qabb,hci)
      endif
********************************************************************************


c big arrays not needed anymore
      deallocate(pija,qabb,edb)

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

      vl=-10.0
      vu=thr
      call ssyevr('V','V','U',nci,hci,nci,vl,vu,il,iu,1.e-6,
     .            nfound,eci,uci,nci,isuppz,
     .            work,lwork,iwork,liwork,info)
      nroot=nfound
      if(nfound.lt.1) stop 'internal error in diag'
      if(sf_s2)then
      ! Compute <S**2> for each SF states
      allocate(Xia(nroot),Xiaib(nroot),Xiaka(nroot),S2(nroot))
      Xia=0.0
      Xiaib=0.0
      Xiaka=0.0
      Do i=1,nroot
        !X_iAaB <iA aB>
        Do j=1,nci
        Xia(i)=Xia(i)+dble(uci(j,i))*
     .         i_alpha_a_beta(iconfb(j,1),iconfb(j,2)-nob)
        enddo
        Xia(i)=Xia(i)**2.0
        !XjAaB XjAbB <iA aB><iA bB>
        Do j=1,nci
          Do kk=1,nci
            if(iconfb(kk,1)==iconfb(j,1))then
            Do k=1,noa
          Xiaib(i)=Xiaib(i)+dble(uci(j,i))*dble(uci(kk,i))*
     .             i_alpha_a_beta(k,iconfb(j,2)-nob)*
     .             i_alpha_a_beta(k,iconfb(kk,2)-nob)
            enddo
            endif
          enddo
        enddo
        !XiAaB XkAbB <iA jB><jA kB>
        Do j=1,nci
          Do kk=1,nci
            if(iconfb(kk,2)==iconfb(j,2))then
            Do k=1,nob
          Xiaka(i)=Xiaka(i)+dble(uci(j,i))*dble(uci(kk,i))*
     .             i_alpha_j_beta(iconfb(j,1),k)*
     .             i_alpha_j_beta(iconfb(kk,1),k)
            enddo
            endif
          enddo
        enddo
        S2(i)=nob+1-sum_ij-Xiaib(i)+Xiaka(i)+Xia(i)+
     .            0.25*(noa-nob-2)*(noa-nob)
      enddo
      deallocate(Xia,Xiaib,Xiaka,i_alpha_j_beta,i_alpha_a_beta)
      endif

      write(*,'(i5,'' SF-roots found, lowest/highest eigenvalue : '',
     .2F8.3,i4)') nroot,eci(1)*27.21139,eci(nroot)*27.21139,info
      if(info.gt.0) stop 'diag error (ssyevr)'

      allocate(umerk(14,nroot))

      ! largest configurations

      Do i=1,nroot
         umax=-1
         kmem=1
      ! first largest
      do j=1,nci
            io=iconfb(j,1)
            iv=iconfb(j,2)
            idum1=max(io,iv)
            idum2=min(io,iv)
            ij=idum2+idum1*(idum1-1)/2
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if(dabs(uu*pp).gt.umax)then
               umax=uu*pp
               imax=io
               jmax=iv
               kmem=j
            endif
      enddo
      umerk(1,i)=imax
      umerk(2,i)=jmax
      umerk(3,i)=dble(uci(kmem,i))
      ! second largest
         umax=-1
         lmem=1
         do j=1,nci
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if((uu*pp).gt.umax.and.j.ne.kmem)then
               umax=uu*pp
               imax=iconfb(j,1)
               jmax=iconfb(j,2)
               lmem=j
            endif
         enddo
      umerk(4,i)=imax
      umerk(5,i)=jmax
      umerk(6,i)=dble(uci(lmem,i))
      ! third largest
         umax=-1
         jmem=-1
         do j=1,nci
            uu=dble(uci(j,i))
            pp=dble(uci(j,i))
            if((uu*pp).gt.umax.and.j.ne.kmem.and.j.ne.lmem)then
               umax=uu*pp
               imax=iconfb(j,1)
               jmax=iconfb(j,2)
               jmem=j
            endif
         enddo
      umerk(7,i)=imax
      umerk(8,i)=jmax
      umerk(9,i)=dble(uci(jmem,i))
      enddo

 20   format(i5,2f9.5,3(' ',i4,'(a) ->',i4,'(b) ',F6.2),f9.5)
 21   format(i5,2f9.5,3(' ',i4,'(a) ->',i4,'(b) ',F6.2))

      write(*,*)'SF states:'
      if(sf_s2)then
      Do i=1, nroot
      write(*,20)i,eci(i),eci(i)*27.21139,
     .int(umerk(1:2,i)),umerk(3,i),
     .int(umerk(4:5,i)),umerk(6,i),
     .int(umerk(7:8,i)),umerk(9,i),S2(i)
      enddo
      else
      Do i=1, nroot
      write(*,21)i,eci(i),eci(i)*27.21139,
     .int(umerk(1:2,i)),umerk(3,i),
     .int(umerk(4:5,i)),umerk(6,i),
     .int(umerk(7:8,i)),umerk(9,i)
      enddo
      endif

      call sf_lresp_ESA(nci,iconfb,maxconfb,xla,yla,zla,mocia,
     .                                      xlb,ylb,zlb,mocib,
     .                noa,nva,nob,nvb,eci,uci,nroot,xmolw,thr)


      if(nto)then


      call SFprint_nto(uci,cca,mocia,nci,nroot,nao,
     .                   noa,nva,ccb,mocib,iconfb,
     .                   maxconfb,nob,nvb)

      endif


      deallocate(hci)



 11   format(i5,f9.3,f8.1,2f11.4)
 12   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3,
     .       3x,'Kshft:',f8.3)
 13   format(                   41x,f8.3,13x,f8.3)
 14   format(3x,F6.2,'(',i4,'a','->',i4,'a',')')
 15   format(3x,F6.2,'(',i4,'b','->',i4,'b',')')
 16   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J:',2f8.3,
     .       3x,'Kshft:',f8.3,2x,'locality:',f6.3)
      call cpu_time(end_time)
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      write(*,*) 'SF-sTDA done.'
      call EXIT(STATUS)

      return
      end SUBROUTINE sfstda



      subroutine sfstdamat(nci,nexb,ncent,noa,nob,nvb,
     .                    mxcnfb,iconfb,dax,edb,
     .                    pija,qabb,hci)
       use omp_lib
      implicit none
      integer, intent(in) :: nci,nexb,ncent,noa,nob,nvb
      integer, intent(in) :: mxcnfb,iconfb(mxcnfb,2)
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*4, intent(out)  :: hci(nci,nci)
      real*8, intent(in)  :: dax,edb(mxcnfb)
      real*4, allocatable :: qj(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,k,iiv,jjv,iwrk,jwrk
      real*4 ek,ej,sdot,ax,de
      allocate(qj(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qk crashed'
! calculate CIS matrix A
      hci=0.0e0
      ax=real(dax)

c beta blocks
!$omp parallel private(i,j,io,iv,jo,jv,iiv,jjv,jwrk,qj,ej)
!$omp do
      do i = 1,nci
         io=iconfb(i,1)
         iv=iconfb(i,2)
         iiv=iv-nob
         do j = 1,i-1
            jo=iconfb(j,1)
            jv=iconfb(j,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
            qj(1:ncent)=pija(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1)
            hci(j,i)=-ej
            hci(i,j)=hci(j,i)
         enddo
         hci(i,i)=real(edb(i))
      enddo
!$omp end do
!$omp end parallel

      deallocate(qj)

      return

      end subroutine sfstdamat


      subroutine sf_ptselect_uks(nexb,ncent,noa,nva,nob,nvb,
     .                  nexptb,mxcnfb,iconfb,
     .                  kconfb,dak,dax,edb,edptb,
     .                  pija,qabb,thrp,newb)
      use omp_lib
      implicit none
      integer, intent(in) :: nexb,ncent,nva,nvb,nexptb
      integer, intent(in) :: noa,nob,mxcnfb
      integer, intent(in) :: kconfb(mxcnfb,2)
      integer, intent (inout) :: iconfb(mxcnfb,2)
      integer, intent (out) :: newb
      real*4, intent(in)  :: pija(ncent,noa*(noa+1)/2)
      real*4, intent(in)  :: qabb(ncent,nvb*(nvb+1)/2)
      real*8, intent(in)  :: dak,dax,thrp
      real*8, intent(in)  :: edptb(mxcnfb)
      real*8, intent(inout) :: edb(mxcnfb)
      real*4, allocatable :: qj(:)
      real*8, allocatable :: ptb(:),pt2b(:)
      integer i,j,k,l,io,iv,jo,jv,ierr,lin,iwrk,jwrk,iiv,jjv,nex
      real*4 ej,sdot,tmpi,tmpj
      real*8 de,pert,amat
      logical, allocatable :: incl_conf(:)

      allocate(ptb(nexb),pt2b(nexb),
     .         qj(ncent),incl_conf(nexptb),stat=ierr)
      if(ierr.ne.0)stop 'allocation for PT intermediates crashed'
      incl_conf=.false.
      qj=0.0
      nex=nexb
      newb = 0
      ptb = 0.0d0
      pt2b = 0.0d0

!$omp parallel private(i,k,io,iv,iiv,iwrk,de,j,jo,jv,jjv,jwrk,ej,
!$omp&                l,qj,ptb,amat,pert) reduction (+:pt2b)
!$omp do
c outer loop over beta S-CSF
      do k=1,nexptb
         io=kconfb(k,1)
         iv=kconfb(k,2)
         iiv=iv-nob
         iwrk=(io-1)*nvb + iiv
         de=edptb(k)
c loop over beta P-CSF
         do j=1,nexb
            jo=iconfb(j,1)
            jv=iconfb(j,2)
            jjv=jv-nob
            jwrk=(jo-1)*nvb + jjv
c J coupling
            qj(1:ncent)=pija(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qabb(1,lin(iiv,jjv)),1)
            ptb(j)=0.0d0
            amat=-ej
c PT2 - e(2) = -((-J)(ia,jb)**2)/(-A(ia,ia))
            ptb(j)=amat**2/(de-edb(j)+1.d-10)
         enddo
         pert= sum(ptb(1:nexb))
c if sum > threshold include the beta S-CSF in the beta P-CSF
         if(pert.gt.thrp)then
            incl_conf(k)=.true.
         else
c else accumulate in the E(PT2) contribution to the alpha and beta P-CSF
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
      do i=1,nexb
         pert=pert+pt2b(i)
         edb(i)=edb(i)-pt2b(i)
      enddo

      write(*,'('' average/max(b) PT2 energy lowering (eV):'',
     .3F10.3)')
     .             27.21139*pert/float(nex),
     .maxval(pt2b)*27.21139


      deallocate(qj,ptb,pt2b)
      return

      end subroutine sf_ptselect_uks
