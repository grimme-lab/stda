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

      SUBROUTINE stda(ncent,nmo,nao,xyz,c,eps,occ,iaoat,thr,thrp,
     .            ax,alphak,betaj,fthr,nvec,nprims)
      use commonlogicals
      use commonresp
      use omp_lib
      IMPLICIT NONE
c input:
      integer ncent,nmo,nao,nprims
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

      integer :: nocc
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
      real*4, allocatable ::q1(:),q2(:)
      real*4 sdot

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

      If(Xcore)then
        do i=1,nmo
           if(occ(i).gt.1.990d0.and.i.le.Ecore2.and.i.ge.Ecore)then
           moci=moci+1
           endif
           if(occ(i).lt.0.010d0.and.eps(i).lt.vthr)moci=moci+1
        enddo

      else
      do i=1,nmo
         if(occ(i).gt.1.990d0.and.eps(i).gt.othr)moci=moci+1
         if(occ(i).lt.0.010d0.and.eps(i).lt.vthr)moci=moci+1
      enddo
      endif

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

      if(Xcore)then
        if(eigvec.or.nto) then ! we want eigenvectors to be printed out
        allocate(vecchk(nmo), stat=ierr)
        if(ierr.ne.0)stop 'allocation failed for vecchk'
          vecchk=0
          moci=0
          do i=1,nmo
             if(occ(i).gt.1.990d0.and.i.le.Ecore2.and.i.ge.Ecore)then
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
             if(occ(i).gt.1.990d0.and.i.le.Ecore2.and.i.ge.Ecore)then
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
      else
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

      nocc=0
      do i=1,nmo
      nocc=nocc+idint(occ(i))
      enddo
      nocc=nocc/2
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
      close(33,status='delete')
      write(*,*)

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
      print '("Transform to MO space = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
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
      call cpu_time(start_time)
      read(40) help
      call cpu_time(end_time)
      print '("Reading sint = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      write(*,*) 'ints done.'
      close(40,status='delete')


      write(*,*) 'S^1/2 ...'
      call cpu_time(start_time)
      call makel(nao,help,x)
      call cpu_time(end_time)
      print '("S^1/2 = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0

      call cpu_time(start_time)
      call dgemm('n','n',nao,moci,nao,1.d0,X,nao,CA,nao,0.d0,SCR,nao)
      write(*,*) 'S^1/2 orthogonalized MO coefficients done.'
      call cpu_time(end_time)
      print '("orthogonalized = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
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
      allocate(q1(ncent),q2(ncent))
      q1=0.0
      q2=0.0
      allocate(qij(ncent,ihomo),qab(ncent,nv),stat=ierr)
      if(ierr.ne.0) stop 'error in diag. J charges allocation'
      qij=0.0
      qab=0.0

      do i=1,ihomo
         call lo12pop(i,i,ncent,nao,iaoat,clow,q2)
         q1(1:ncent)=q1(1:ncent)+q2(1:ncent)
         qij(1:ncent,i)=q2(1:ncent)
      enddo
      write(*,'(/'' SCF atom population (using active MOs):'')')
      write(*,'(10F7.3)')q1(1:ncent)*2.0
      write(*,*)
      write(*,'('' # electrons in TDA:'',F8.3)') 2.0*sum(q1(1:ncent))
      write(*,*)
      do i=no+1,moci
         j=i-no
         call lo12pop(i,i,ncent,nao,iaoat,clow,q2)
         qab(1:ncent,j)=q2(1:ncent)
      enddo

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

c     call prmat4(6,gamj,ncent,ncent,'gamj')
c     call prmat4(6,gamk,ncent,ncent,'gamk')
c     stop

      imem1=no
      imem2=nv
      imem3=imem1*imem2*2
      imem1=imem1*(imem1+1)/2
      imem2=imem2*(imem2+1)/2
      imem3=imem3+imem1+imem2
      hilf=dble(imem3)/1024.0**2
      hilf=dble(4*ncent)*hilf
      imem1=idint(hilf)
      write(*,*)'memory needed for q data (Mb) ',imem1
      write(*,*)'computing q(ij,n) ...'



! compute intermediates q's refer to charges, p's to gam*q (i.e., contracted)

!     Coulomb type terms: calc qij*gam^J
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


      nex=no*nv
      allocate(pia(ncent,nex),qia(ncent,nex), stat=ierr)
      if(ierr.ne.0)stop 'allocation failed for intermediates'
      qia=0.0
      pia=0.0
!       K-type terms
      k=0
!$omp parallel private(i,j,k,q2)
!$omp do
      do i=1,no
         do j=no+1,moci
!            k=k+1
            k=(i-1)*nv+j-no
            call lo12pop(i,j,ncent,nao,iaoat,clow,q2)
            qia(1:ncent,k)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      pia=0.0
      call ssymm('l','l',ncent,nex,1.0,gamk,ncent,qia,ncent,0.0
     .             ,pia,ncent)
      deallocate(gamk)

c determine singles which are lower than thr
      k=0
      j=0
      l=0
      i=0
      do io=1,no ! occ loop
         do iv=no+1,moci ! virt loop
            de=epsi(iv)-epsi(io)
            l=iv-no
            ej=0.0d0
            ej=dble(uci(l,io))
            de=de-ej
            ek=0.0d0
            i=nv*(io-1)+l
            q1(1:ncent)=pia(1:ncent,i)
            ek=sdot(ncent,q1,1,qia(1,i),1)
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
      deallocate(uci)
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
      write(*,*)'        eV     # centers'
      j=max(1,no-10)
      do i=j,no
         xc=0
         do k=1,ncent
            xc=xc+qij(k,i)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsi(i)*27.21139,1./(xc+1.d-8)
      enddo
      write(*,*)
      j=min(moci,no+11)
      do i=no+1,j
         xc=0
         do k=1,ncent
            xc=xc+qab(k,i-no)**2
         enddo
         write(*,'(i4,F10.3,F8.1)') i,epsi(i)*27.21139,1./(xc+1.d-8)
      enddo


      write(*,*)
      write(*,*)'            lowest CSF states'
      write(*,*)'      eV     nm      excitation i->a               eV'
      do i=1,min(nex,25)
         io=iconf(i,1)
         iv=iconf(i,2)
         q1(1:ncent)=qij(1:ncent,io)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1,1,0.0,q2,1)
         jii=sdot(ncent,q1,1,q2,1)
         k=iv-no
         q1(1:ncent)=qab(1:ncent,k)
         ej=sdot(ncent,q1,1,q2,1)
         call ssymv('l',ncent,1.0e0,gamj,ncent,q1,1,0.0,q2,1)
         jaa=sdot(ncent,q1,1,q2,1)
         l=nv*(io-1)+k
         q1(1:ncent)=pia(1:ncent,l)
         q2(1:ncent)=qia(1:ncent,l)
         ek=sdot(ncent,q1,1,q2,1)
! de is now the Kia shift
         de=0
         loc=ej/sqrt(jii*jaa) ! locality
         if(dokshift) call kshift_to_ediag(de,ek)
         write(*,14) i,27.211*ed(i),
     .   1.d+7/(ed(i)*2.19474625d+5),iconf(i,1:2),
     .   27.211*(epsi(iv)-epsi(io)),27.211*ej,27.211*ek,27.211*de,loc
      enddo
 14   format(i5,f6.2,f8.1,       5x,i4,' ->',i4,5x,'gap,J,K:',3f8.3,
     .       3x,'Kshft:',f8.3,2x,'locality:',f6.3,E12.5)


      deallocate(qij,qab)
      ihilf=no*(no+1)/2
      allocate(qij(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qij allocation'
      qij=0.0
      ij=0
!$omp parallel private(i,j,ij,q2)
!$omp do
      do i=1,no
         do j=1,i
            ij=lin(i,j)
            call lo12pop(i,j,ncent,nao,iaoat,clow,q2)
            qij(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      allocate(pij(ncent,ihilf), stat=ierr)
        if(ierr.ne.0)stop 'allocation failed for (ij| intermediate'
        pij=0.0
        call ssymm('l','l',ncent,ihilf,1.0,gamj,ncent,qij,ncent,0.0
     .             ,pij,ncent)


      q2=0.0
      ihilf=nv*(nv+1)/2
      allocate(qab(ncent,ihilf),stat=ierr)
      if(ierr.ne.0) stop 'error in qab allocation'
      ij=0
!$omp parallel private(i,j,l,k,ij,q2)
!$omp do
      do i=no+1,moci
         do j=no+1,i
           k=i-no
           l=j-no
           ij=lin(k,l)
           call lo12pop(i,j,ncent,nao,iaoat,clow,q2)
           qab(1:ncent,ij)=q2(1:ncent)
         enddo
      enddo
!$omp end do
!$omp end parallel
      if(XsTD.eqv..false.)then
      deallocate(clow) ! orthogonalized MO coefficients not needed any more
      endif
!      call cpu_time(hilf)
!      write(*,*) 'time elapsed:',hilf-time

      write(*,*)
      write(*,*)'selecting CSF ...'

      deallocate(q1,q2)
      call ptselect(nex,ncent,no,nv,nexpt,maxconf,iconf,kconf,
     .              ak,ax,ed,edpt,pia,qia,pij,qab,thrp,new)
      deallocate(edpt,kconf)

c nroot at this point is the number of primary CSF. The
c number of roots in the range 0-thr (as desired) is not
c known but will be determined by the diag routine.
      nroot=nex
      write(*,*)new,'CSF included by PT.'
      nci=nex+new
      write(*,*)nci,'CSF in total.'

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
      if(XsTD)then
      if(RSH_flag)then
      if(RSH_sub)then
      call xstd_rpamat_RSH2(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      else
      call xstd_rpamat_RSH(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      endif
      else
      call xstd_rpamat(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,apb,ambsqr,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      endif
      if(FULL2PA.eqv..false.)then
      deallocate(clow)
      endif
      else
      if(full)then
      write(*,*)'(  |  ) integrals are computed exactly'
      call rrpamat_full(nci,ncent,no,nv,maxconf,iconf,
     .             ak,ax,ed,apb,ambsqr,ca,nao,moci,nprims,epsi)
      elseif(direct_full)then
      call rrpamat_full_direct(nci,ncent,no,nv,maxconf,iconf,
     .             ak,ax,ed,apb,ambsqr,ca,nao,moci,nprims,epsi)
      else
      call rrpamat(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed,pia
     .             ,qia,pij,qab,apb,ambsqr)
      endif
      endif

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
      if(velo_OR.eqv..false.)call optrot(nci,apb,amb,iconf,maxconf,
     .xl,yl,zl,moci,no,nv,xm,ym,zm,xmolw)
      if(velo_OR.eqv..true.)call optrot_velo(nci,apb,amb,iconf,maxconf,
     .xv,yv,zv,moci,no,nv,xm,ym,zm,xmolw)
      call cpu_time(end_time)
      print '("Opt. Rot.   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      if(nto)then
      nnto=num_freq
      allocate(uci(nci,nnto))
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in print NROs
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      Do i=1,3
      if(i==1)write(*,*)'X axis'
      if(i==2)write(*,*)'Y axis'
      if(i==3)write(*,*)'Z axis'
      Do ii=1,nnto
      write(dummy,'(i1,a,i0)')i,'-',ii
      open(unit=14,file=dummy)
      read(14,*)uci(:,ii)
      close(14,status='delete')
      enddo
      call print_nto_resp_new(uci,ca,moci,nci,nnto,nao,iconf,maxconf,
     .no,nv,i,xm,ym,zm,.true.)
      enddo
      endif
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
      print '("sTD-DFT Time = ",f12.2," minutes.")'
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

      call lresp(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
     .                     no,nv)

      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      if(nto)then
      nnto=num_freq+1
      allocate(uci(nci,nnto))
      write(*,*)
      write(*,*)'====================================================
     .=================='
      write(*,*)'               Welcome in print NROs
     . program'
      write(*,*)'====================================================
     .=================='
      write(*,*)
      Do i=1,3
      if(i==1)write(*,*)'X axis'
      if(i==2)write(*,*)'Y axis'
      if(i==3)write(*,*)'Z axis'
      Do ii=1,nnto
      write(dummy,'(i1,a,i0)')i,'-',ii
      open(unit=14,file=dummy)
      read(14,*)uci(:,ii)
      close(14,status='delete')
      enddo
      call print_nto_resp(uci,ca,moci,nci,nnto,nao,iconf,maxconf,no,nv,i
     .,xl,yl,zl,.false.)
      enddo
      endif
!       call lresp_noinv(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv)
!       call lresp_noinv1(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
!      .                     no,nv)
!
!       call cpu_time(end_time)
!       print '("Lresp   Time = ",f12.2," minutes.")'
!      .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      CALL EXIT(STATUS)
      endif



c big arrays not needed anymore
      deallocate(pij,qab,pia,ed)
      if(.not.printexciton)deallocate(qia)
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
      if(FULL2PA)then
      if(triplet) stop 'not available'
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53,status='delete')
      call lresp_2PA_full(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
     .                     no,nv,eci,uci,hci,nroot,ncent,ax,nao,clow)
      deallocate(clow)
      call cpu_time(end_time)
      print '("Lresp   Time = ",f12.2," minutes.")'
     .      ,(end_time-start_time)/60.0
      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      write(*,*)
      CALL EXIT(0)
      else
      if(TPA)then
      if(triplet) stop 'not available'
      call cpu_time(start_time)
      allocate( amb(nci*(nci+1)/2), stat=ierr )
      if(ierr.ne.0)stop 'allocation failed for A-B'
      open(unit=53,file='amb',form='unformatted',status='old')
      read(53) amb
      close(53,status='delete')
      call lresp_2PA_SP(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,
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

      print '("sTD-DFT Time = ",f12.2," minutes.")'
     .      ,(end_time-stda_time)/60.0
      write(*,*)
      CALL EXIT(0) !important because vectors are not normalized the same way
      endif

      endif
      deallocate(apb,ambsqr)
!********************************************************************************
      open(unit=53,file='amb',form='unformatted',status='old')
      close(53,status='delete')
      else
      deallocate(gamj,qij)
!********************************************************************************
cccccccccccccccccccccccccc
c standard TDA procedure c
cccccccccccccccccccccccccc

!********************************************************************************
! construct ( 0.5 * B ) for X trafo (velocity correction) and print to file *
!********************************************************************************
      if(velcorr) then
        if(XsTD)then
        if(RSH_flag)then
        call RSH_Xsrtdacorr(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed
     .                ,clow,nao,moci)
        else
        call Xsrtdacorr(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed
     .                ,clow,nao,moci)
        endif
        else
        call rtdacorr(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed
     .                ,pia,qia,pij,qab)
        endif
      endif
!********************************************************************************

      allocate( hci(nci,nci), stat=ierr  )
      if(ierr.ne.0)stop 'allocation failed for TDA matrix'
      write(*,*)'calculating TDA matrix ...'
      if(XsTD)then
      if(RSH_flag)then
      if(RSH_sub)then
      call Xstda_mat_RSH2(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      else
      call Xstda_mat_RSH(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      endif
      else
      call Xstda_mat(nci,ncent,no,nv,maxconf,iconf,
     .ak,ax,ed,hci,alphak,betaj,xyz,nao,moci,clow,alphak,betaj,
     .epsi)
      endif
      deallocate(clow)
      else
      call rtdamat(nci,ncent,no,nv,maxconf,iconf,ak,ax,ed,pia
     .             ,qia,pij,qab,hci)
      endif
!      call prmat4(6,hci,nci,nci,'A-Matrix')
!********************************************************************************

c big arrays not needed anymore
      deallocate(pij,qab,pia,ed)
      if(.not.printexciton)deallocate(qia)

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
      if(SOS_2PA)then
      call tpa_sos(nroot,nci,eci,uci,hci,xl,yl,zl,moci,
     .maxconf,iconf,no,nv) !idem, it gives also information about how lresp relaxes
      endif

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
      write(*,*) 'sTD-DFT done.'
      else
      write(*,*) 'sTDA done.'
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
*                                                                     *
*  make the loewdin orthogonalization matrix x = u' s1/2 u            *
*                                                                     *
***********************************************************************

      subroutine makel(nao, s, x)
      use omp_lib
      use commonlogicals
      implicit none
      real*8, intent (in) :: s(nao*(nao+1)/2)
      real*8, intent (out) :: x(nao,nao)
      integer, intent (in) :: nao
      real*8, allocatable ::aux(:),vecs(:,:),e(:),cc(:,:)
      integer lwork,k,i,j,lin,m,info

      lwork  = 1 + 6*nao + 2*nao**2
      allocate (vecs(nao,nao),e(nao),aux(lwork),cc(nao,nao))

      k=0
      do i=1,nao
         do j=1,i
            k=k+1
            x(j,i)=s(k) ! the dsyev routine stores only upper triangular matrix
            x(i,j)=s(k)
         enddo
      enddo

      call dsyev ('V','U',nao,x,nao,e,aux,lwork,info)

c     call dHQRII(s,nao,nao,e,vecs)

      do i=1,nao

         if(bfw)then
         if(e(i).lt.0)then
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) 'normaly must stop in S^1/2!'
         write(*,*) 'because e(',i,')=',e(i)
         write(*,*) 'assuming that e(i)=0'
         write(*,*) 'YOU ASK FOR IT'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         e(i)=0.0
         endif
         else
         if(e(i).lt.0)then
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) 'must stop in S^1/2!'
         write(*,*) 'because e(',i,')=',e(i)
         write(*,*) 'If e(i) is near 0 try -BFW option'
         write(*,*) 'and check your output !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         write(*,*) '!!!!! BIG FAT WARNING !!!!!'
         stop
         endif
         endif
         e(i)=dsqrt(e(i))
      enddo
!$omp parallel private (i,m)
!$omp do
      do m=1,nao
         do i=1,nao
         vecs(i,m)=     x(i,m)
         cc(i,m)=e(m)*x(i,m)
         enddo
      enddo
!$omp end do
!$omp end parallel

      call dgemm('N','T',nao,nao,nao,1.0d0,vecs,
     .                   nao,cc,nao,0.0d0,x,nao)

      deallocate(e,aux,cc,vecs)

      return
      end

***********************************************************************

      subroutine shrink(n,x,y)
      integer n,i,j,k,lin
      real*8 x(n,n),y(n*(n+1)/2)
      k=0
      do i=1,n
         do j=1,i
            k=k+1
            y(k)=x(j,i)
         enddo
      enddo
      return
      end

***********************************************************************
c for single precision
      subroutine shrink4(n,x,y)
      integer n,i,j,k,lin
      real*4 x(n,n),y(n*(n+1)/2)
      k=0
      do i=1,n
         do j=1,i
            k=k+1
            y(k)=x(j,i)
         enddo
      enddo
      return
      end

***********************************************************************
* atomic Löwdin population for a (transition) density with MOs ij
***********************************************************************

      subroutine lo12pop(i,j,ncent,nao,iaoat,ca,q)
      IMPLICIT NONE
      integer i,j,nao,ncent,iaoat(*)
      real*8 ca(*)
      real*4 q(*)

      integer k,ii,jj,iat

      q(1:ncent) = 0.0

      ii=(i-1)*nao
      jj=(j-1)*nao
      do k=1,nao
         iat=iaoat(k)
         q(iat)=q(iat)+ca(k+ii)*ca(k+jj)
      enddo

      end
***********************************************************************

***********************************************************************
* shift diagonal elements by K(ia,ia) -> xTB option
***********************************************************************
      subroutine kshift_to_ediag(de,ek)
      use kshiftcommon
      implicit none
      real*8, intent(inout) :: de ! A(ia,ia)
      real*8, intent (in) :: ek ! K(ia,ia)
      real*8  shft,wau,sau,ekia,dmp
      ! we apply a rationally damped, K(ia,ia) dependent shift to A(ia,ia)
      !
      ! f [ K(ia,ia) ]  =   shiftmax /( 1 + (K(ia,ia)/width)**steepness
      !
      ! A(ia,ia)=A(ia,ia) + f [ K(ia,ia) ]

      ! shiftmax and width are given eV, convert to a.u.
      sau=shftmax*0.036749320d0
      wau=shftwidth*0.036749320d0
      shft = 1.0d0 + (ek/wau)**shftsteep
      de = de + sau/shft ! shift A(ia,ia)

      end

***********************************************************************
* select important CSFs beyond energy threshold
***********************************************************************
      subroutine ptselect(nex,ncent,no,nv,nexpt,mxcnf,iconf,kconf,
     .                    dak,dax,ed,edpt,pia,qia,pij,qab,thrp,new)
      use omp_lib
      implicit none
      integer, intent(in) :: nex,ncent,nexpt,mxcnf
      integer, intent(in) :: kconf(mxcnf,2),no,nv
      integer, intent (inout) :: iconf(mxcnf,2)
      integer, intent (out) :: new
      real*4, intent(in)  ::qia(ncent,mxcnf),pia(ncent,mxcnf)
      real*4, intent(in) ::qab(ncent,nv*(nv+1)/2),pij(ncent,no*(no+1)/2)
      real*8, intent(in)  :: dak,dax,thrp,edpt(mxcnf)
      real*8, intent(inout) :: ed(mxcnf)
      real*4, allocatable :: qj(:),qk(:)
      real*8, allocatable :: pt(:),pt2(:)
      integer i,j,k,io,iv,jo,jv,ierr,lin,iwrk,jwrk,iiv,jjv
      real*4 ek,ej,sdot,tmpi,tmpj,ak,ax
      real*8 de,pert,amat
      logical, allocatable :: incl_conf(:)
      ak=real(dak)
      ax=real(dax)
      allocate(qj(ncent),qk(ncent),pt(nex),pt2(nex),incl_conf(nexpt),
     .         stat=ierr)
      if(ierr.ne.0)stop 'allocation for PT intermediates crashed'
      incl_conf=.false.
      new=0
      pt2=0.0d0
      qk=0.0
      qj=0.0

! loop over secondary/neglected CSF, this is done omp-parallel
!$omp parallel private(i,k,io,iv,iiv,iwrk,qk,de,j,jo,jv,jjv,jwrk,ek,ej,
!$omp&                 qj,pt,amat,pert) reduction (+:pt2)
!$omp do
      do k=1,nexpt
         io=kconf(k,1)
         iv=kconf(k,2)
         iiv=iv-no
         iwrk=(io-1)*nv + iiv
         qk(1:ncent)=pia(1:ncent,iwrk)
         de=edpt(k)
c loop over primary CSF
         do j=1,nex
            jo=iconf(j,1)
            jv=iconf(j,2)
            jjv=jv-no
            jwrk=(jo-1)*nv + jjv
            ek=sdot(ncent,qk,1,qia(1,jwrk),1)
c coupling
            qj(1:ncent)=pij(1:ncent,lin(io,jo))
            ej=sdot(ncent,qj,1,qab(1,lin(iiv,jjv)),1)
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


      deallocate(qj,qk,pt,pt2,incl_conf)
      return

      end subroutine ptselect
***********************************************************************


***********************************************************************
* set up A+B and A-B (packed form) in RKS case
***********************************************************************
      subroutine rrpamat(nci,ncent,no,nv,mxcnf,iconf,dak,dax,ed,
     .                  pia,qia,pij,qab,apb,ambsqr)
      use commonlogicals
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*4, intent(in)  :: qia(ncent,mxcnf),pia(ncent,mxcnf)
      real*4, intent(in)  :: pij(ncent,no*(no+1)/2)
      real*4, intent(in)  :: qab(ncent,nv*(nv+1)/2)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: apb(nci*(nci+1)/2),ambsqr(nci*(nci+1)/2)
      real*4, allocatable :: qk(:),qj(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk
      real*4 ek,ej,sdot,ak,ax,de
      allocate(qk(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for q1 crashed'
      ak=real(dak)
      ax=real(dax)
! calculate A+B and A-B
      apb=0.0e0
      ambsqr=0.0e0

! if ax=0, A-B is diagonal and its off-diagonal elements do not need to be calculated
      if(abs(dax).lt.1.0d-6) then
        ij=0
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,ek,de)
!$omp do
        do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           qk(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv+jjv
              ek=sdot(ncent,qk,1,qia(1,jwrk),1) ! ek = (ia|jb)
              apb(ij)=2.0*ak*ek
              ambsqr(ij)=0.0
           enddo ! j
           de=real(ed(i))
           ij=lin(i,i)
           ek=sdot(ncent,qk,1,qia(1,iwrk),1)
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
        allocate(qj(ncent), stat=ierr)
        ij=0
        ! for now ambsqr=A+B and apb=A-B, since we need to take the sqrt of A-B (but want to save memory)
!$omp parallel private(ij,i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
        do i=1,nci
           io=iconf(i,1)
           iv=iconf(i,2)
           iiv=iv-no
           iwrk=(io-1)*nv + iiv
           qk(1:ncent)=pia(1:ncent,iwrk)
           do j=1,i-1
              ij=lin(i,j)
              jo=iconf(j,1)
              jv=iconf(j,2)
              jjv=jv-no
              jwrk=(jo-1)*nv+jjv
              ek=sdot(ncent,qk,1,qia(1,jwrk),1) ! ek = (ia|jb)
              qj(1:ncent)=pij(1:ncent,lin(io,jo))
              ej=sdot(ncent,qj,1,qab(1,lin(iiv,jjv)),1) !  ej = (ij|ab)
              ambsqr(ij)=2.0*ak*ek
              jwrk=(io-1)*nv+jjv
              qj(1:ncent)=pia(1:ncent,jwrk)
              jwrk=(jo-1)*nv+iiv
              ek=sdot(ncent,qj,1,qia(1,jwrk),1) ! now ek = (ib|aj), results from Fock-exchange, thus we scale by ax
              ambsqr(ij)=ambsqr(ij)-ax*ek-ej
              apb(ij)=ax*ek-ej
           enddo ! j
           ij=lin(i,i)
           ek=sdot(ncent,qk,1,qia(1,iwrk),1)
           apb(ij)=real(ed(i))-ak*ek+ax*ek    ! diagonal element of A-B
           ambsqr(ij)=real(ed(i))-ax*ek+ak*ek ! diagonal element of A+B
        enddo ! i
!$omp end do
!$omp end parallel




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
        deallocate(qj)
      endif ! GGA/hybrid case


      deallocate(qk)
      return

      end subroutine rrpamat
***********************************************************************

***********************************************************************
* set up A matrix in RKS case
***********************************************************************
      subroutine rtdamat(nci,ncent,no,nv,mxcnf,iconf,dak,dax
     .                   ,ed,pia,qia,pij,qab,hci)
      use omp_lib
      implicit none
      integer, intent(in) :: nci,ncent,no,nv,mxcnf,iconf(mxcnf,2)
      real*4, intent(in)  :: qia(ncent,mxcnf),pia(ncent,mxcnf)
      real*4, intent(in)  :: pij(ncent,no*(no+1)/2)
      real*4, intent(in)  :: qab(ncent,nv*(nv+1)/2)
      real*8, intent(in)  :: dak,dax,ed(mxcnf)
      real*4, intent(out) :: hci(nci,nci)
      real*4, allocatable :: qk(:),qj(:)
      integer i,j,ij,io,iv,jo,jv,ierr,lin,iiv,jjv,iwrk,jwrk
      real*4 ek,ej,sdot,ak,ax,de
      allocate(qk(ncent),qj(ncent), stat=ierr)
      if(ierr.ne.0)stop 'allocation for qkj crashed'
      ak=real(dak)
      ax=real(dax)
! calculate CIS matrix A
      hci=0.0e0
!$omp parallel private(i,j,io,iv,jo,jv,iiv,iwrk,jjv,jwrk,qk,qj,ek,ej)
!$omp do
      do i=1,nci
         io=iconf(i,1)
         iv=iconf(i,2)
         iiv=iv-no
         iwrk=(io-1)*nv + iiv
         qk(1:ncent)=pia(1:ncent,iwrk)
         do j=1,i-1
            jo=iconf(j,1)
            jv=iconf(j,2)
            jjv=jv-no
            jwrk=(jo-1)*nv + jjv
            ek=0.0
            ek=sdot(ncent,qk,1,qia(1,jwrk),1)
            qj(1:ncent)=pij(1:ncent,lin(io,jo))
            ej=0.0
            ej=sdot(ncent,qj,1,qab(1,lin(iiv,jjv)),1)
            hci(j,i)=ak*ek-ej
            hci(i,j)=hci(j,i)
         enddo
         hci(i,i)=real(ed(i))
      enddo
!$omp end do
!$omp end parallel
      deallocate(qk,qj)
      return
      end subroutine rtdamat
***********************************************************************


      subroutine setrep(rep)
      real*8 rep(94)

cSemiempirical Evaluation of the GlobalHardness of the Atoms of 103
cElementsof the Periodic Table Using the MostProbable Radii as their Size
cDescriptors
cDULAL C. GHOSH, NAZMUL ISLAM
c2009 in Wiley InterScience (www.interscience.wiley.com).DOI 10.1002/qua.22202

c values in the paper multiplied by two because (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
c definition they use is 1/2d^2 E/dN^2 (in Eh)
      rep( 1)=  0.472592880d0
      rep( 2)=  0.922033910d0
      rep( 3)=  0.174528880d0
      rep( 4)=  0.257007330d0
      rep( 5)=  0.339490860d0
      rep( 6)=  0.421954120d0
      rep( 7)=  0.504381930d0
      rep( 8)=  0.586918630d0
      rep( 9)=  0.669313510d0
      rep(10)=  0.751916070d0
      rep(11)=  0.179641050d0
      rep(12)=  0.221572760d0
      rep(13)=  0.263485780d0
      rep(14)=  0.305396450d0
      rep(15)=  0.347340140d0
      rep(16)=  0.389247250d0
      rep(17)=  0.431156700d0
      rep(18)=  0.473082690d0
      rep(19)=  0.171054690d0
      rep(20)=  0.202762440d0
      rep(21)=  0.210073220d0
      rep(22)=  0.217396470d0
      rep(23)=  0.224710390d0
      rep(24)=  0.232015010d0
      rep(25)=  0.239339690d0
      rep(26)=  0.246656380d0
      rep(27)=  0.253982550d0
      rep(28)=  0.261288630d0
      rep(29)=  0.268594760d0
      rep(30)=  0.275925650d0
      rep(31)=  0.307629990d0
      rep(32)=  0.339315800d0
      rep(33)=  0.372359850d0
      rep(34)=  0.402735490d0
      rep(35)=  0.434457760d0
      rep(36)=  0.466117080d0
      rep(37)=  0.155850790d0
      rep(38)=  0.186493240d0
      rep(39)=  0.193562100d0
      rep(40)=  0.200633110d0
      rep(41)=  0.207705220d0
      rep(42)=  0.214772540d0
      rep(43)=  0.221846140d0
      rep(44)=  0.228918720d0
      rep(45)=  0.235986210d0
      rep(46)=  0.243056120d0
      rep(47)=  0.250130180d0
      rep(48)=  0.257199370d0
      rep(49)=  0.287847800d0
      rep(50)=  0.318486730d0
      rep(51)=  0.349124310d0
      rep(52)=  0.379765930d0
      rep(53)=  0.410408080d0
      rep(54)=  0.441057770d0
      rep(55)=  0.050193320d0
      rep(56)=  0.067625700d0
      rep(57)=  0.085044450d0
      rep(58)=  0.102477360d0
      rep(59)=  0.119911050d0
      rep(60)=  0.137327720d0
      rep(61)=  0.154762970d0
      rep(62)=  0.172182650d0
      rep(63)=  0.189612880d0
      rep(64)=  0.207047600d0
      rep(65)=  0.224467520d0
      rep(66)=  0.241896450d0
      rep(67)=  0.259325030d0
      rep(68)=  0.276760940d0
      rep(69)=  0.294182310d0
      rep(70)=  0.311595870d0
      rep(71)=  0.329022740d0
      rep(72)=  0.345922980d0
      rep(73)=  0.363880480d0
      rep(74)=  0.381305860d0
      rep(75)=  0.398774760d0
      rep(76)=  0.416142980d0
      rep(77)=  0.433645100d0
      rep(78)=  0.451040140d0
      rep(79)=  0.468489860d0
      rep(80)=  0.485845500d0
      rep(81)=  0.125267300d0
      rep(82)=  0.142686770d0
      rep(83)=  0.160116150d0
      rep(84)=  0.177558890d0
      rep(85)=  0.194975570d0
      rep(86)=  0.212407780d0
      rep(87)=  0.072635250d0
      rep(88)=  0.094221580d0
      rep(89)=  0.099202950d0
      rep(90)=  0.104186210d0
      rep(91)=  0.142356330d0
      rep(92)=  0.163942940d0
      rep(93)=  0.185519410d0
      rep(94)=  0.223701390d0

      return
      end


****************************************************************************
* Given orbital indeces i1 and i2, lin() returns index in the linear array *
****************************************************************************

      integer function lin(i1,i2)
      integer i1,i2
      integer idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      lin=idum2+idum1*(idum1-1)/2
      return
      end

      integer*8 function lin8(i1,i2)
      integer i1,i2
      integer*8 idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      lin8=idum2+idum1*(idum1-1)/2
      return
      end

      subroutine cofc(nat,xyz,coc)
!      calculates center of nuclear charge
      implicit none
      integer, intent (in) :: nat
      real*8, intent (out) :: coc(3)
      real*8, intent (in) :: xyz(4,nat)
      integer iat,iatyp
      real*8  ctot
!      real*8 atwt(111)
!      data  atwt / 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107,
!     .        14.0067, 15.9994,  18.9984032, 20.1797, 22.98976928,
!     .        24.305, 26.9815386, 28.0855, 30.973762, 32.065,
!     .        35.453, 39.948, 39.0983, 40.078, 44.955912, 47.867,
!     .        50.9415, 51.9961, 54.938045, 55.845, 58.933195,
!     .        58.6934, 63.546, 65.38, 69.723, 72.64, 74.9216,
!     .        78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
!     .        91.224, 92.90638, 95.96, 98.0, 101.07, 102.90550,
!     .        106.42, 107.8682, 112.411, 114.818, 118.71, 121.76,
!     .        127.6, 126.90447, 131.293, 132.9054519, 137.327,
!     .        138.90547, 140.116, 140.90765, 144.242, 145.0, 150.36,
!     .        151.964, 157.25, 158.92535, 162.5, 164.93032,
!     .        167.259, 168.93421, 173.054, 174.9668, 178.49,
!     .        180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,
!     .        196.966569, 200.59, 204.3833, 207.2, 208.9804, 209.0,
!     .        210.0, 222.0, 223.0, 226.0, 227.0, 232.03806,
!     .        231.03588, 238.02891, 237.0, 244.0, 243.0, 247.0,
!     .        247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0,
!     .        267.0, 268.0, 271.0, 272.0, 270.0, 276.0, 281.0, 280.0 /

       iatyp=0
       ctot=0.0d0
       coc(1:3)=0.0d0
       ctot=xyz(4,1)
       coc(1:3)=ctot*xyz(1:3,1)

       do iat=2,nat
          ctot=ctot+xyz(4,iat)
          coc=coc+xyz(4,iat)*xyz(1:3,iat)
       end do
       coc=coc/ctot

      end subroutine cofc

      subroutine print_tdadat(nroot,xmolw,energy,trstr,thr,fname)
      implicit none
      integer, intent( in )      :: nroot
      real*4, intent( in )       :: energy(*)
      real*8, intent( in )       :: trstr(4,nroot),thr,xmolw
      character*80, intent( in ) :: fname
      integer i

      open(unit=26,file=trim(fname))
      write(26,*)'NM'
      write(26,*)'UV'
      write(26,*)'MMASS'
      write(26,*)xmolw
      write(26,*)'LFAKTOR'
      write(26,*)' 0.5'
      write(26,*)'RFAKTOR'
      write(26,*)' 1.0'
      write(26,*)'WIDTH'
      write(26,*)' 0.20'
      write(26,*)'SHIFT'
      write(26,*)' 0.00'
      write(26,*)'DATXY'
      do i=1,nroot
         if(energy(i).gt.thr) cycle
         write(26,'(i4,F10.4,4f13.6)')
     . i,energy(i)*27.21139,trstr(1,i),trstr(2,i),trstr(3,i),
     . trstr(4,i)
      enddo
      close (26)

      end subroutine print_tdadat
