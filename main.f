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
      program acis_prog
      use stdacommon ! mostly input and primitive data
      use kshiftcommon ! kshiftvariables
      use commonlogicals
      use commonresp
      implicit real*8 (a-h,o-z)

      real*8, allocatable ::cc(:)

      integer, allocatable :: ccspin(:)
      real*8, allocatable ::xyz(:,:)
      real*8  xx(10),alpha,beta,ptlim
      character*79 dummy
      character*79 fname
      character*8 method
      integer imethod,inpchk,mform,nvec
      logical molden,da,chkinp,xtbinp
      integer, dimension(8) :: datetimevals
      !!! libcint
      double precision, allocatable :: overlap_AO(:),mu(:,:),mag(:,:)
      double precision, allocatable :: velo(:,:),quadrupole(:,:)
      real*8,allocatable :: help(:),help1(:),help2(:),help3(:)

      call date_and_time(VALUES=datetimevals)
      print '(I0,"-",I0,"-",I0,1X,I0,":",I0,":",I0,".",I3)',
     .      datetimevals(1:3), datetimevals(5:8)

      write(*,'(//
     .          17x,''*********************************************'')')
      write(*,'(17x,''*                                           *'')')
      write(*,'(17x,''*               s  t  d  2                  *'')')
      write(*,'(17x,''*                                           *'')')
      write(*,'(17x,''*                S. Grimme                  *'')')
      write(*,'(17x,''* Mulliken Center for Theoretical Chemistry *'')')
      write(*,'(17x,''*             Universitaet Bonn             *'')')
      write(*,'(17x,''*                                           *'')')
      write(*,'(17x,''*             M. de Wergifosse              *'')')
      write(*,'(17x,''*                MOST/IMCN                  *'')')
      write(*,'(17x,''*     Universite Catholique de Louvain      *'')')
      write(*,'(17x,''*                                           *'')')
      write(*,'(17x,''*               Version 2.0                 *'')')
      write(*,'(17x,''*     Tue 07 Jan 2025 09:09:41 AM CET       *'')')
      write(*,'(17x,''*********************************************'')')
      write(*,*)
      write(*,'('' Please cite:'')')
      write(*,*)
      write(*,'('' For sTDA and sTD-DFT,'')')
      write(*,'('' S. Grimme, J. Chem. Phys. 138 (2013) 244104'')')
      write(*,'('' C. Bannwarth, S. Grimme Comput. Theor. Chem.
     . 1040–1041 (2014) 45-53'')')
      write(*,'('' M. de Wergifosse, S. Grimme, J. Phys. Chem A
     . 125 (2021) 18 3841-3851'')')
      write(*,*)
      write(*,'('' For XsTDA and XsTD-DFT,'')')
      write(*,'('' M. de Wergifosse, S. Grimme, J. Chem. Phys.
     . 160 (2024) 204110'')')
      write(*,'('' M. de Wergifosse, J. Phys. Chem. Lett. 15
     . (2024) 51 12628–12635'')')
      write(*,*)
      write(*,'('' With contributions from:'')')
      write(*,'('' C. Bannwarth, P. Shushkov, P. Beaujean'')')
      write(*,*)
      write(*,'(a,a)')'===============================================',
     .                 '======================='
      write(*,*)
c defaults

      ptlim=1.7976931348623157d308 ! energy range that will be scanned by PT (we use just a large number)
      thre=7.0                     ! energy range for primary CSF space
      alpha=-100.0d0               ! alpha & beta are large negative numbers and can be changed by user input
      beta=-100.0d0                ! otherwise global hybrid defaults will be used

c the following value yields converged UV spectra for several members of
c of the DYE12 set
      thrp=1.d-4

      mform=1 ! mform is the "style" specifier for Molden input, by default try TM input: ORCA/XTB = 0, TM=1,Molpro=2, Terachem/Gaussian=3

      rpachk=.false.  ! sTD-DFT wanted?
      triplet=.false. ! triplet excitations wanted?
      eigvec=.false. ! eigenvector printout wanted?
      nvec=0 ! if so, how many vecs?
      printexciton=.false. ! print information for exciton coupling program
      velcorr=.true. ! by default: use  Rv(corr) in TDA
      aniso=.false. ! print  anisotropic f/R

      chkinp=.false. ! perform input check?
      fname='molden.input'
      xtbinp=.false. !use xtbinput?
      screen=.false. ! prescreen configurations (by Schwarz-type screening)

      ! Kia shifting defaults
      dokshift=.false.
      shftmax=0.50d0 ! maximum Kia shift in eV
      shftwidth=0.10d0 ! damping threshold in eV
      shftsteep=4.0d0 ! steepness

c read the tm2xx file, otherwise (-f option) the tm2molden file
      molden=.true.
      ax=-1
      imethod=0
      inpchk=0
      resp=.false.
      TPA=.false.
      aresp=.false.
      ESA=.false.
      smp2=.false.
      bfw=.false.
      spinflip=.false.
      sf_s2=.false.
      rw=.false.
      rw_dual=.false.
      pt_off=.false.
      optrota=.false.
      XsTD=.false.
      cint=.true.
      RSH_flag=.false.
      SOS_2PA=.false.
      Xcore=.false.
      FULL2PA=.false.
      multipole=.false.

! check for input file
      inquire(file='.STDA',exist=da)
      if(da)then
        call readinp(ax,thre,alpha,beta)
      endif

      do i=1,command_argument_count()
      call getarg(i,dummy)
      if(index(dummy,'-fo').ne.0)then
         call getarg(i+1,fname)
         molden=.false.
         inpchk=inpchk+1
      endif
      if(index(dummy,'-f').ne.0.and.index(dummy,'-fo').eq.0)then
         call getarg(i+1,fname)
         molden=.true.
         inpchk=inpchk+1
      endif

      if(index(dummy,'-ax').ne.0)then     ! get amout of Fock-exchange
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) ax=xx(1)
      endif
      if(index(dummy,'-e').ne.0.and.index(dummy,'-exc').eq.0)then      ! get energy threshold
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) thre=xx(1)
      endif
      if(index(dummy,'-p').ne.0)then      ! get PT threshold
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) thrp=10.0**(-xx(1))
      endif
      if(index(dummy,'-sty').ne.0)then      ! inputform of Molden input
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) mform=dnint(xx(1))
      endif
      if(index(dummy,'-al').ne.0)then      ! get alpha -> parameter for K term
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) alpha=xx(1)
      endif

      if(index(dummy,'-be').ne.0)then      ! get beta -> parameter for J term
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) beta=xx(1)
      endif

      if(index(dummy,'-lpt').ne.0)then ! PT limit -> negelct CSFs beyond this completely
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) ptlim=xx(1)
      endif

      if(index(dummy,'-xtb').ne.0) then ! xTB input, set defaults
                                        ! two other dirty paramters (not
                                        ! really fitted in function
                                        ! kshift_to_ediag
         ! set sTDA parameters
         ax=0.50d0
         alpha=2.0d0  ! the La/Lb splitting becomes bad for values > 2.5
         beta= 4.0d0
         ! set Kshift parameters
         dokshift=.true.
         shftmax=0.500d0 ! maximum Kia shift in eV
         shftmax_somo=1.000d0 ! maximum additional Kia shift in eV for *->SOMOs and SOMOs->* excitations
         shftwidth=0.10d0 ! damping threshold in eV
         shftsteep=4.0d0 ! steepness
         ! input settings
         mform=0
         inpchk=1
         molden=.false.
         xtbinp=.true.
         chkinp=.false.

      endif
      if(index(dummy,'-kshift').ne.0) dokshift=.true.
      if(index(dummy,'-t').ne.0) triplet=.true. ! triplet excitations wanted
      if(index(dummy,'-rpa').ne.0) rpachk=.true. ! do stddft

      if(index(dummy,'-resp').ne.0)then ! Do response function
         resp=.true.
         rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         num_freq=int(xx(1))
      endif

      if(index(dummy,'-aresp').ne.0) then
         aresp=.true.
         rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         num_freq=int(xx(1))
      endif

      if(index(dummy,'-2PA').ne.0)then ! Do response function
         TPA=.true.
         rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         num_trans=int(xx(1))
      endif

      if(index(dummy,'-SOS2PA').ne.0)then ! Do response function
         SOS_2PA=.true.
         rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         num_trans=int(xx(1))
      endif

      if(index(dummy,'-hxc2PA').ne.0)then ! Do response function
         FULL2PA=.true.
         rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         num_trans=int(xx(1))
         XsTD=.true.
         write(*,*) 'hxc2PA, only with XsTD-DFT'
         cint=.true.
      endif

      if(index(dummy,'-s2s').ne.0)then ! Do response function
         ESA=.true.
         !rpachk=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         state2opt=int(xx(1))
      endif

!       if(index(dummy,'-MP2').ne.0)then ! Do mp2
!          smp2=.true.
!       endif

      if(index(dummy,'-rw').ne.0)then
         rw=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(xx(1)==1) pt_off=.true.
      endif

      if(index(dummy,'-dual').ne.0)then
         rw_dual=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(xx(1)==1) pt_off=.true.
      endif

      if(index(dummy,'-BFW').ne.0)then
         bfw=.true.
      endif

      if(index(dummy,'-sf').ne.0)then ! Do spinflip
         spinflip=.true.
         if(xtbinp) then
         beta= 3.0d0
         ax=0.36d0
         endif
      endif
      if(index(dummy,'-spin').ne.0)then
      sf_s2=.true.
      endif

      if(index(dummy,'-oprot').ne.0)then
      rpachk=.true.
      optrota=.true.
      velo_OR=.false.
      call getarg(i+1,dummy)
      call readl(79,dummy,xx,nn)
      if(xx(1)==1)velo_OR=.true.
      endif

      if(index(dummy,'-nto').ne.0)then ! Do nto
         nto=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         Nnto=int(xx(1))
      endif

      if(index(dummy,'-XsTD').ne.0)then
      XsTD=.true.
      if(xtbinp.eqv..true.)stop 'xTB with XsTD is not implemented'
      cint=.true.
      write(*,*)'**********************'
      write(*,*)'You choose to use XsTD'
      write(*,*)'**********************'
      !dokshift=.false.
         if(rpachk.eqv..false.)then
      write(*,*)'Velocity correction deactivated with XsTDA by default'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-libcintOFF').ne.0)then
      cint=.false.
      endif

      if(index(dummy,'-CORE').ne.0)then
      Xcore=.true.
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         Ecore=int(xx(1))
         call getarg(i+2,dummy)
         call readl(79,dummy,xx,nn)
         Ecore2=int(xx(1))
      write(*,*)'Doing core to valence excitations'
      write(*,*)'using an occupied space from'
      write(*,*)Ecore, ' to ',Ecore2
      endif
!      if(index(dummy,'-RSH').ne.0)then ! Do range-separated hybrid only with XsTD
!         RSH_flag=.true.
!         XsTD=.true.
!         cint=.true.
!         dokshift=.false.
!         call getarg(i+1,dummy)
!         call readl(79,dummy,xx,nn)
!         if(nn.gt.0) RSH=xx(1)
!         call getarg(i+2,dummy)
!         call readl(79,dummy,xx,nn)
!         if(nn.gt.0) RSH_ax=xx(1)
!         call getarg(i+3,dummy)
!         call readl(79,dummy,xx,nn)
!         if(nn.gt.0) RSH_beta=xx(1)
!         call getarg(i+4,dummy)
!         call readl(79,dummy,xx,nn)
!         if(nn.gt.0) RSH_sub=xx(1)
!         write(*,*)RSH_sub
!         RSH_sub=.true.
!         write(*,*)'Using a long-range corrected functional with XsTD'
!         write(*,*)'Range-separated hybrid parameter is ', RSH,RSH_ax,
!     .   RSH_beta,RSH_sub
!      endif
 12   format(a,3f9.4,L)
      if(index(dummy,'-CAMB3LYP').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.33d0
         RSH2=RSH
         RSH_ax=0.19d0
         RSH_beta=0.46d0
         RSH_sub=.true.
         ax=0.38d0
         alpha=0.9d0
         beta=1.86d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using CAM-B3LYP with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'Velocity correction deactivated with XsTDA by default'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-wB97XD2').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.2d0
         RSH2=RSH
         RSH_ax=0.222d0
         RSH_beta=0.888d0
         RSH_sub=.true.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using wB97X-D with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'Velocity correction deactivated with XsTDA by default'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-wB97XD3').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.25d0
         RSH2=RSH
         RSH_ax=0.1957d0
         RSH_beta=0.8043d0
         RSH_sub=.true.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using wB97X-D3 with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'Velocity correction deactivated with XsTDA by default'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-wB97MV').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.3d0
         RSH2=RSH
         RSH_ax=0.15d0
         RSH_beta=0.85d0
         RSH_sub=.true.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using wB97MV with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'Velocity correction deactivated with XsTDA by default'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-Bvel').ne.0)then
      velcorr=.true.
      write(*,*)'Velocity correction asked for XsTDA'
      endif

      if(index(dummy,'-oldwB97XD3').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.25d0
         RSH2=RSH
         RSH_ax=0.1957d0
         RSH_beta=1.0d0
         RSH_sub=.false.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using oldwB97X-D3 with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'No velocity correction with XsTDA/RSH for the moment'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-oldwB97MV').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.3d0
         RSH2=RSH
         RSH_ax=0.15d0
         RSH_beta=1.0d0
         RSH_sub=.false.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using oldwB97MV with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'No velocity correction with XsTDA/RSH for the moment'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-oldwB97XD2').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH=0.2d0
         RSH2=RSH
         RSH_ax=0.222d0
         RSH_beta=1.0d0
         RSH_sub=.false.
         ax=0.51d0
         alpha=4.51d0
         beta=8.0d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using oldwB97X-D with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'No velocity correction with XsTDA/RSH for the moment'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-SRC2R1').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH2=0.69d0 !short range mu
         RSH=1.02d0  !long  range mu
         RSH_ax=0.55d0
         RSH_beta=0.08d0
         RSH_sub=.false.
         !B3LYP parameters for the CSF screening
         ax=0.2d0
         alpha=1.516d0
         beta=0.566d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using SRC2-R1 with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'No velocity correction with XsTDA/RSH for the moment'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-SRC2R2').ne.0)then ! Do range-separated hybrid only with XsTD
         RSH_flag=.true.
         XsTD=.true.
         cint=.true.
         dokshift=.false.
         RSH2=2.20d0 !short range mu
         RSH=1.80d0  !long  range mu
         RSH_ax=0.91d0
         RSH_beta=0.28d0
         RSH_sub=.false.
         !B3LYP parameters for the CSF screening
         ax=0.2d0
         alpha=1.516d0
         beta=0.566d0
         write(*,*)'**********************'
         write(*,*)'You choose to use XsTD'
         write(*,*)'**********************'
         write(*,*)'Using SRC2-R2 with XsTD'
         write(*,12)'Range-separated hybrid parameters are ', RSH,
     .   RSH_ax,RSH_beta,RSH_sub
         if(rpachk.eqv..false.)then
      write(*,*)'No velocity correction with XsTDA/RSH for the moment'
         velcorr=.false.
         endif
      endif

      if(index(dummy,'-intfull').ne.0)then
      if(XsTD.eqv..true.)stop 'you cannot ask for full
     .integrals and XsTD at the same time'
      full=.true.
      cint=.true.
      rpachk=.true.
      dokshift=.false.
      endif

      if(index(dummy,'-directintfull').ne.0)then
      if(XsTD.eqv..true.)stop 'you cannot ask for full
     .integrals and XsTD at the same time'
      direct_full=.true.
      cint=.true.
      rpachk=.true.
      dokshift=.false.
      endif

      if(index(dummy,'-chk').ne.0) chkinp=.true. ! do input check
      if(index(dummy,'-vectm').ne.0)then
         eigvec=.true. ! print eigenvectors
         call getarg(i+1,dummy)
         call readl(79,dummy,xx,nn)
         if(nn.gt.0) nvec=dnint(xx(1))
      endif
      ! print transition dipole moments
      if(index(dummy,'-excprint').ne.0) printexciton=.true.
      if(index(dummy,'-oldtda').ne.0) velcorr=.false.
      if(index(dummy,'-aniso').ne.0) aniso=.true.
      enddo

!      if(ptlim.gt.1.0d308) ptlim = 3.0d0*thre ! set ptlimit to a multiple of thre

      if(chkinp) mform=0 ! if input check is done, start with 0

ccccccccccccccccccccccccccccccc
c check the input
ccccccccccccccccccccccccccccccc
      if(inpchk.gt.1) stop 'please specify only one input file!'
      if(inpchk.lt.1) stop 'no input file specified!'
      if(molden) then
       write(*,*) 'reading a molden input...'
      else if (xtbinp) then
       write(*,*) 'reading an xTB output...'
      else
       write(*,*) 'reading a tm2stda file...'
       chkinp=.false.
      endif
      if(ax.lt.0.and..not.chkinp) then
        stop 'specify Fock exchange via -ax <real>'
      endif

      if(rpachk) velcorr=.false.

! set imethod to 1 if old reading version is used = only RKS possible
      if(.not.molden.and.imethod.lt.1) imethod=1

  888 continue
ccccccccccccccccccccccccccccccc
c first call to get dimensions
ccccccccccccccccccccccccccccccc
      if(molden)then
        inpchk=0 ! use this integer now to determine UKS/RKS
        call readmold0(ncent,nmo,nbf,nprims,fname,inpchk)
        if(imethod.eq.0)imethod=inpchk ! if UKS/RKS has not been specified
c compare input and inpchk
        if(inpchk.ne.imethod)stop'U/R-KS input doesnot match moldenfile'

      else if(xtbinp) then ! read parameters from xTB input
       call readxtb0(imethod,ncent,nmo,nbf,nprims)

      else
        call readbas0a(0,ncent,nmo,nbf,nprims,fname)
      endif

      if(nprims.eq.0.or.ncent.eq.0.or.nmo.eq.0)
     .stop 'read error'

ccccccccccccccccccccccccccccccc
c allocate mo vectors
ccccccccccccccccccccccccccccccc
       icdim = nmo*nbf
      if(imethod.eq.2.and..not.molden) then
       icdim=2*nmo*nbf
       nmo = 2*nmo
      endif

      allocate(cc(icdim),stat=ierr)
      if(ierr.ne.0) stop 'allocation failed in main for cc'

      allocate(ccspin(nmo),stat=ierr)
      if(ierr.ne.0) stop 'allocation failed in main for ccspin'

*****************************
* allocate common variables *
*****************************
      allocate(co(ncent,4),exip(nprims),cxip(nprims),occ(nmo),eps(nmo))
      allocate(ipat(nprims),ipty(nprims),ipao(nprims),iaoat(nbf))
      allocate(atnam(ncent),eta(nprims,25))

!      if(ncent.gt.maxat) then
!         write(*,*)ncent,maxat
!         stop 'too many centers'
!      endif
!      if(nmo   .gt.ndi22) then
!         write(*,*) nmo,ndi22
!         stop 'too many mos    '
!      endif
!      if(nprims.gt.ndi22) then
!         write(*,*) nprims,ndi22
!         stop 'too many prims  '
!      endif

ccccccccccccccccccccccccccccccccc
c read vectors and basis and ..
ccccccccccccccccccccccccccccccccc
      if(molden)then
        call readmold(mform,imethod,ncent,nmo,nbf,nprims,cc,
     .  ccspin,icdim,fname)
      else if(xtbinp) then
        call readxtb(imethod,ncent,nmo,nbf,nprims,cc)
        if(imethod.eq.2) then
          do i=1,nmo/2
            ccspin(i)=1
          enddo
          do i=nmo/2+1,nmo
            ccspin(i)=2
          enddo
        endif
      else
       if(imethod.eq.1) call readbasa(1,imethod,ncent,nmo,nbf,nprims,cc,
     .icdim,fname,iaobas)
       if(imethod.eq.2) call readbasb(1,imethod,ncent,nmo,nbf,nprims,cc,
     .ccspin,icdim,fname,iaobas)
      endif
      if(imethod.eq.1) deallocate( ccspin )

ccccccccccccccccccccccccccccccccc
c precalculate primitive data
ccccccccccccccccccccccccccccccccc
      if(cint)then
      write(*,*)'*****************************'
      write(*,*)'*using libcint integral deck*'
      write(*,*)'*****************************'
      write(*,*) 'SIZEOF(int)', SIZEOF(nbf)
      allocate(overlap_AO(nbf*(nbf+1)/2))
      call overlap(ncent,nprims,nbf,overlap_AO)
      allocate(mu(nbf*(nbf+1)/2,1:3))
      call dipole_moment(ncent,nprims,nbf,mu)
! unfortunately libcint does not give the same result for the magnetic moment.
! thus, we still use the old integral deck for the magnetic moment
!       allocate(mag(nbf*(nbf+1)/2,1:3))
!       call magnetic_moment(ncent,nprims,nbf,mag)
      call intslvm2(ncent,nmo,nbf,nprims)
      nao=nbf
      if(nao.eq.0)nao=nprims
      allocate(velo(nbf*(nbf+1)/2,1:3))
      call velo_moment(ncent,nprims,nbf,velo)
      open(unit=40,file='sint', form='unformatted',status='replace')
      open(unit=31,file='xlint',form='unformatted',status='replace')
      open(unit=32,file='ylint',form='unformatted',status='replace')
      open(unit=33,file='zlint',form='unformatted',status='replace')
!       open(unit=34,file='xmint',form='unformatted',status='replace')
!       open(unit=35,file='ymint',form='unformatted',status='replace')
!       open(unit=36,file='zmint',form='unformatted',status='replace')
      open(unit=37,file='xvint',form='unformatted',status='replace')
      open(unit=38,file='yvint',form='unformatted',status='replace')
      open(unit=39,file='zvint',form='unformatted',status='replace')
      write(40)overlap_AO
      close(40)
      write(31)mu(:,1)
      write(32)mu(:,2)
      write(33)mu(:,3)
      close(31)
      close(32)
      close(33)
!       write(34)mag(:,1)
!       write(35)mag(:,2)
!       write(36)mag(:,3)
!       close(34)
!       close(35)
!       close(36)
      write(37)velo(:,1)
      write(38)velo(:,2)
      write(39)velo(:,3)
      close(37)
      close(38)
      close(39)
      deallocate(overlap_AO,mu,velo)
      else
      call intslvm(ncent,nmo,nbf,nprims)
      nao=nbf
      if(nao.eq.0)nao=nprims
      endif



      if(nto)then
      !
      ! Have everything necessary to compute nice NTOs
      !
      open(unit=11,file='NTOao')
      Do i=1, nprims
      write(11,21)ipat(i),ipty(i),ipao(i),exip(i),cxip(i)
      enddo
      close(11)
      open(unit=11,file='NTOat')
      Do i=1,ncent
      write(11,*)atnam(i),co(i,1:4)
      enddo
      close(11)
      open(unit=11,file='NTOvar')
      write(11,*)ncent,nprims,nmo,icdim,nbf,imethod
      close(11)
      if(imethod==2)then
      open(unit=11,file='NTOspin')
      Do i=1,nmo
      write(11,*)ccspin(i)
      enddo
      close(11)
      endif
      !call molden_file(ncent,nprims,nmo,icdim,nbf,imethod,cc,ccspin)
      if(xtbinp.or.mform==1)then
      open(unit=12,file='fnorm')
      Do i=1,nprims
      write(12,*)1.0
      enddo
      close(12)
      endif
      else
      open(unit=12,file='fnorm')
      close(12,status='delete')
      endif

 21   format(3i10,3x,2f24.9)

 11   format(3i5,2f9.4)

      deallocate(ipty,exip,cxip,atnam,eta)


cccccccccccccccccccccccccccccccccccccc
c  if input check (-chk) is used     c
c use calculated ovlp from intslvm   c
c to check input (Mulliken-pop)      c
cccccccccccccccccccccccccccccccccccccc
      if(chkinp) then
        call mulpopcheck(nbf,nmo,cc,occ,idum)
        if(idum.ne.0) then ! restart from input read if norm is wrong
          mform=mform+1
          close(36,status='delete')
          deallocate( cc )
          if(imethod.eq.2) deallocate( ccspin )
          if(mform.gt.3) stop 'unreadable format'
          write(*,*)'Restarting and trying different input style...'
          deallocate(iaoat,occ,eps,co)
          goto 888
        else
          call inputcheck_printout(mform,.true.,imethod,nbf,nmo)
          call exit(0) ! leave program
        endif
      else
         if(.not.xtbinp) then
           call inputcheck_printout(mform,.false.,imethod,nbf,nmo)
         endif
      endif
cccccccccccccccccccccccccccccccccccccccccc
c make an additional printout whether    c
c restricted or unrestricted calculation c
c will be done                           c
cccccccccccccccccccccccccccccccccccccccccc
       write(*,*)' '
      if (imethod.eq.1) then
       write(*,'(A)',advance='no')'Restricted MOs found...'
       if(rpachk) then
        write(*,'(A)',advance='yes')' RKS-sTD-DFT will be performed'
       else
        write(*,'(A)',advance='yes')' RKS-sTDA will be performed'
       endif
      else
       write(*,'(A)',advance='no')'Unrestricted MOs found...'
       if(rpachk) then
        write(*,'(A)',advance='yes')' UKS-sTD-DFT will be performed'
       else
        write(*,'(A)',advance='yes')' UKS-sTDA will be performed'
       endif
      endif
       write(*,*) ' '
ccccccccccccccccccccccccccccccccc
c stda
ccccccccccccccccccccccccccccccccc
      allocate(xyz(4,ncent),stat=ierr)
      if(ierr.ne.0) stop 'allocation failed in main for xyz'

      do i=1,ncent
      xyz(1:4,i)=co(i,1:4)
      enddo
      deallocate(co)

!       if(smp2)then
!       call sMP(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,
!      .        ax,alpha,beta)
!       endif

      if(spinflip)then
      call sfstda(ncent,nmo,nao,xyz,cc,eps,occ,ccspin,iaoat,thre,
     .           thrp,ax,alpha,beta,ptlim,nvec)
      endif

      if (imethod.eq.1) then
      if(rw)then
      call stda_rw(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,thre,
     .        thrp,ax,alpha,beta,ptlim,nvec,nprims)
      else
      if(rw_dual)then
      call stda_rw_dual(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,thre,
     .        thrp,ax,alpha,beta,ptlim,nvec,ipat,ipao,nprims)
      else
      call stda(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,thre,
     .        thrp,ax,alpha,beta,ptlim,nvec,nprims)
      endif
      endif
      else
      call sutda(ncent,nmo,nao,xyz,cc,eps,occ,ccspin,iaoat,thre,
     .           thrp,ax,alpha,beta,ptlim,nvec)
      endif
      deallocate(ipat,ipao)

      call date_and_time(VALUES=datetimevals)
      print '(I0,"-",I0,"-",I0,1X,I0,":",I0,":",I0,".",I3)',
     .       datetimevals(1:3), datetimevals(5:8)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine byteout(s,mem)
      integer*8 mem
      character*(*)s
      if(mem/1024**3.gt.0) then
      write(*,'('' memory needed for '',a,'':'',i4,'' Gb'')')
     .s,mem/1024**3
      return
      endif
      if(mem/1024**2.gt.0) then
      write(*,'('' memory needed for '',a,'':'',i4,'' Mb'')')
     .s,mem/1024**2
      return
      endif
      if(mem/1024.gt.0) then
      write(*,'('' memory needed for '',a,'':'',i4,'' Kb'')')s,mem/1024
      return
      endif
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine performs a mulliken population analysis to check if the
! input data is correct
      subroutine mulpopcheck(nbf,nmo,c,occ,ireturn)
      implicit none
      real*8, intent( in ) :: c(nbf*nmo),occ(nmo)
      real*8 summe1,summe2
      integer i,j,k
      integer, intent( in ) :: nmo,nbf
      integer, intent( out ) :: ireturn
      real*8, allocatable ::pmul(:,:),ovvec(:),ovmat(:,:),dmat(:,:)
      allocate( dmat(nbf,nbf), ovmat(nbf,nbf), ovvec(nbf*(nbf+1)/2),
     .         pmul(nbf,nbf) )
      ireturn=0
      write(*,*) ' '
      call header('I N P U T   C H E C K',0)
! for debugging
!      k=0
!      do i=1,nmo
!        write(*,*) 'MO', i
!        do j=1,nbf
!        k=k+1
!        write(*,*)j,c(k)
!        enddo
!      enddo

! check whether input data is reasonable
      open(unit=40,file='sint',form='unformatted',status='old')
      read(40)ovvec ! read overlap matrix elements
      close(40)

      call blow(nbf,ovvec,ovmat)  ! blow up from vector to matrix
!      call prmat(6,ovvec,nbf,0,'ovlp') ! optional: printout
      deallocate(ovvec) ! free memory
      ! construct density matrix
      do i=1,nbf
        do j=1,i
          dmat(i,j)=0.0d0
          dmat(j,i)=0.0d0
          do k=1,nmo
            dmat(i,j)=dmat(i,j)+occ(k)*c(nbf*(k-1)+i)*c(nbf*(k-1)+j)
            dmat(j,i)=dmat(i,j)
          enddo
        enddo
      enddo
!      call prmat(6,dmat,nbf,nbf,'density matrix')
      summe2=0.0d0
      ! get number of electrons in the system (reference)
      do i=1,nmo
        summe2=summe2+occ(i)
      enddo
      call dsymm('l','l',nbf,nbf,1.d0,dmat,nbf,ovmat,nbf,0.d0,pmul,nbf) ! P * S , for Mulliken population
      ! now compute # of electrons from  tr( P * S ) and compare to reference
      summe1=0.0d0
      do i=1,nbf
        summe1=summe1+pmul(i,i)
      enddo
!      call prmat(6,pmul,nbf,nbf,'P * S')
      deallocate(dmat,ovmat,pmul) ! free memory
      if(dabs(summe1-summe2).gt.nmo*1.d-5)then ! since Gaussian comes with accuracy up to 5 decimal points, take a system-dependent accuracy
         write(*,*) 'Number of electrons from Mulliken population'
         write(*,'(a,x,f18.6,x,f18.6)') 'does not match:',summe1,summe2
         write(*,*) 'WARNING: Input style is incorrect!'
         write(*,*) ' '
         ireturn=1
         return
      endif
! print out which input format was detected
      if(ireturn.eq.0)then
          write(*,*)'Mulliken population correct'
          write(*,'(f18.6,x,a,x,f18.6)')summe1,'vs.',summe2
          write(*,*) ' '
      endif
      end subroutine mulpopcheck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! input check printouts !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inputcheck_printout (ires,typ,imethod,nbf,nmo)
      implicit none
      integer, intent( in ) :: ires,imethod,nbf,nmo
      logical typ
!     check run
      if (typ) then
        write(*,*) ' '
        write(*,*) ' --- S U C C E S S --- '
        write(*,*) ' '
        select case(ires)
         case(0)
          write(*,*) ' GTO/MO data matches ORCA style'
          write(*,*) ' use "-sty 0" flag in the future'
         case(1)
          if(imethod*nbf.gt.nmo) then
           write(*,*) ' GTO/MO data matches TURBOMOLE style'
           write(*,*) ' use "-sty 1" flag in the future'
          else
           write(*,*) ' GTO/MO data matches Cartesian basis style'
           write(*,*) '  (compatible with TURBOMOLE/MOLPRO/G09)'
           write(*,*) '  use the following flag in the future'
           write(*,*) '   TURBOMOLE: "-sty 1" '
           write(*,*) '   GAUSSIAN 09: "-sty 1" '
           write(*,*) '   MOLPRO: "-sty 2" '
          endif
         case(2)
          write(*,*) ' GTO/MO data matches MOLPRO style'
          write(*,*) ' use "-sty 2" flag in the future'
         case(3)
          write(*,*) 'GTO/MO data matches TERACHEM/G09 style'
          write(*,*) ' use "-sty 3" flag in the future'
        end select
!   erase AO-files  before exiting program
        open(unit=40,file='sint',form='unformatted',status='old')
        close(40,status='delete')
        open(unit=31,file='xlint',form='unformatted',status='old')
        close(31,status='delete')
        open(unit=32,file='ylint',form='unformatted',status='old')
        close(32,status='delete')
        open(unit=33,file='zlint',form='unformatted',status='old')
        close(33,status='delete')
        open(unit=34,file='xmint',form='unformatted',status='old')
        close(34,status='delete')
        open(unit=35,file='ymint',form='unformatted',status='old')
        close(35,status='delete')
        open(unit=36,file='zmint',form='unformatted',status='old')
        close(36,status='delete')
        open(unit=37,file='xvint',form='unformatted',status='old')
        close(37,status='delete')
        open(unit=38,file='yvint',form='unformatted',status='old')
        close(38,status='delete')
        open(unit=39,file='zvint',form='unformatted',status='old')
        close(39,status='delete')
        return
      else
        ! standard run: no check
        write(*,*) ' '
        write(*,*) 'Skipping input check...'
        select case(ires)
           case(0)
           write(*,*) 'Assuming ORCA style w/ s,p functions (-sty 0)'
           case(1)
           if(imethod*nbf.gt.nmo) then
            write(*,*) 'Assuming TURBOMOLE style (-sty 1)'
           else
            write(*,*) 'Assuming TM/MOLPRO/G09 style (-sty 1)'
           endif
           case(2)
            write(*,*) 'Assuming MOLPRO style (-sty 2)'
           case(3)
            write(*,*) 'Assuming TERACHEM style (-sty 3)'
           case(4)
            write(*,*) 'Assuming spherical MOs (-sty 4)'
            write(*,*) ' ...this is format is not tested yet!'
        end select
        return
      endif
      end subroutine inputcheck_printout
