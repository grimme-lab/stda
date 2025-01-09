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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MOLDEN readout routine                !
! by C. Bannwarth, University of Bonn   !
! Thu Oct 23 10:44:18 CEST 2014         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some general notes:
!   1) The Molden input is not unique, i.e. the QC programs write it differently.
!   We therefore make a sanity check by means of a Mulliken population analysis (see main.f).
!   2) Molden uses a different arrangement of f-functions than stda, therefore when the MO
!   coefficients are read out, they need to be reordered for every set of f-functions.
!   That is why a marking vector "ifstart" is introduced which is declared when the AO basis is set up and
!   used as a checkpoint when to do the reordering in the MO readout part (which is done after reading the MO
!   coefficients for one set of f-functions)
!   3) Related to 1): We work with Cartesian basis functions, i.e. the dx**2,dy**2,dz**2 are not linearly independent
!   (similar for f-functions, except fxyz). This means that the GTO contraction coefficients need to be scaled by certain prefactors.
!   But these prefactors are somehow included into the MO coefficents (different in certain programs) and therefore, these need to be post-processed.

c reads the molden input file just to get the dimensions (actual readout follows later on)
      subroutine readmold0(ncent,nmo,nbf,nprims,wfn,idum)
      use strings
      implicit none

      character*(*), intent( in ) :: wfn
      integer, intent( out ) :: ncent,nmo,nbf,nprims

      character*79 line
      character*25 args(10)
      character*1 aang
      character*2 aatom

      integer iwfn,iatom,idum,jatom,mbasf,kbasf,lang,ibas,iostatus,nargs
      integer i,j,maxlen
      logical ex

      real*8 xdum,ydum,zdum

      write(*,*)
      write(*,*)'reading: ',wfn
      call header('M O / A O    I N P U T ',0)

      inquire(file=wfn,exist=ex)
      if(.not.ex)then
          write(*,*)'file:',wfn,' not found'
          stop 'input file not found'
      endif
      iwfn=29
      open(unit=iwfn,file=wfn,status='old')

c initialize some variables
      iatom=0
      idum=0
cccccccccccccccccccccccccccccccccccccccccccccc
c determine ncent from molden.input          c
cccccccccccccccccccccccccccccccccccccccccccccc
      ncent=0
      jatom=0
       call findstr('[Atoms]',7,iwfn,1,2)
       do
        read(iwfn,*,IOSTAT=iostatus)aatom,iatom,idum,xdum,ydum,zdum
        if(iostatus.ne.0) exit
        if(iatom.le.jatom) exit
        jatom=iatom
        ncent=ncent+1
       enddo
      if(jatom.ne.ncent) stop 'error in ncent read - atoms not in range'
      rewind(iwfn)

cccccccccccccccccccccccccccccccccccccccccc
c make sure, there is no spherical basis c
cccccccccccccccccccccccccccccccccccccccccc
        do
           read(iwfn,'(A)',IOSTAT=iostatus) line
           if(iostatus.lt.0) exit
           call parse(line,' ',args,nargs)
           if(nargs.eq.1) then
            if(args(1).eq.'[5D]'.or.args(1).eq.'[5D10F]'.or.args(1).eq.
     .         '[5D7F]') then
             stop 'program requires Cartesian basis'
            endif
           endif
        enddo
        rewind(iwfn)

cccccccccccccccccccccccccccccccccccccccccccccc
c determine nbf and nprims from molden.input c
cccccccccccccccccccccccccccccccccccccccccccccc
      nbf=0
      nprims=0
      nmo=0
! search for [GTO] statement to read basis
      call findstr('[GTO]',5,iwfn,1,1)
! go through atoms and read number of aos,caos and prims
      jatom=0
      do i=1,ncent
  555   read(iwfn,*,iostat=iostatus)iatom,idum ! 1st, read atom identifier
        if(iostatus.ne.0) goto 555 ! if it does not fit, read next line
        if(ncent.lt.iatom) stop 'Error: no atoms does not fit to basis'
        if(jatom.lt.iatom)jatom=iatom
! read line, check whether s,p,d,f orbital is present (aang) and read no of prims (ibas)
        do
           read(iwfn,'(A)',IOSTAT=iostatus) line
           if(iostatus.ne.0) exit
! use string module
           call parse(line,' ',args,nargs)
           if(nargs.lt.2)exit
! special case of 'sp' declaration
           if(trim(args(1)).eq.'sp'.or.trim(args(1)).eq.'SP') then

             call value(args(2),ibas,iostatus)
             if(iostatus.ne.0) stop 'ibas: error in arg to integer conv'
             kbasf=4 ! a set of s and p orbitals, i.e., 1 + 3
             mbasf=4
           else
! general case
             aang=args(1)
             call value(args(2),ibas,iostatus)
             if(iostatus.ne.0) stop 'ibas: error in arg to integer conv'
             call aangchk(0,aang,mbasf,kbasf,lang)
             if(lang.lt.0)cycle

           endif

           do j=1,ibas
             read(iwfn,*)
           enddo
         nprims=nprims+kbasf*ibas
         nbf=nbf+kbasf
         nmo=nmo+mbasf ! preliminary MO number
        enddo
        backspace(iwfn) ! to be sure, go back one line
      enddo
! sanity check
      if(jatom.ne.ncent) stop 'int. error in ncent read!'

ccccccccccccccccccccccccccccc
! determine number of MOs   c
ccccccccccccccccccccccccccccc
      rewind(iwfn)
      mbasf=0 ! mbasf is now the counter for mos
      call findstr('[MO]',4,iwfn,1,1)

      do
       read(iwfn,'(A)',IOSTAT=iostatus) line
       if(iostatus.ne.0) exit
       call parse(line,' ',args,nargs)
        if(nargs.ne.2) cycle
        if(args(1).eq.'Ene=') then ! if energy is found increase counter
        mbasf=mbasf+1
        endif
      enddo

      nmo=mbasf

ccccccccccccccccccccccccccccccc
! check for closed/open shell c
ccccccccccccccccccccccccccccccc
      idum=1
      if(nmo.gt.nbf) then
       rewind(iwfn)
       call findstr('[MO]',4,iwfn,1,1)
      do
       read(iwfn,'(A)',IOSTAT=iostatus) line
       if(iostatus.lt.0) exit
       line=lowercase(line)
       call removesp(line)
       call parse(line,'=',args,nargs)
        if(nargs.ne.2) cycle
        if(index(line,'spin').ne.0) then
         if(index(line,'beta').ne.0) then
          idum=2
          exit
         else
          cycle
         endif
        endif
      enddo
      endif
! close molden input file
      rewind(iwfn)
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

cccccccccccccccccccccccccccccc
! actual read out routine    c
cccccccccccccccccccccccccccccc
      subroutine readmold(mform,imethod,ncent,nmo,nbf,nprims,cc,ccspin
     .,icdim,wfn)
      use strings
      use stdacommon
      implicit double precision (a-h,o-z)

      dimension cc(icdim)
      integer ccspin(nmo),ifstart(nbf)
      integer, intent( in ) :: mform,imethod
      integer, intent( in ) :: ncent,nmo,nbf,nprims

      character*(*), intent( in ) :: wfn
      character*25 args(10)
      character*1 aang
      character*79 line,dummy,dummy2,prntfrmt
      integer iostatus,maxlen

! data statements to read GTO basis
      data pi32 /5.56832799683170d0/
      data pt187 /1.875d+00/
      data pt5,pt75,pt656 /0.5d0,0.75d0,6.5625d0/


! determine length of ncent integer (for printout to prevent ***)
      maxlen=0
      call lenint(ncent,maxlen)
      prntfrmt=' '
      write(prntfrmt,'(a,i0,a)')'(2x,a2,x,i',maxlen,
     .               ',2x,3f14.8,3x,f10.2)'
      iwfn=29
      open(unit=iwfn,file=wfn,status='old')
!      if(imethod.eq.2) nmo = 2*nmo ! not necessary - done automatically by read out

ccccccccccccccccccccccccccccccccccccccccccccccccc
! read in atoms & coordinates                   c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      call findstr('[Atoms]',7,iwfn,1,2)
      backspace(iwfn)
! either coordinates are given in a.u.
      read(iwfn,'(A)',IOSTAT=iostatus) line
      call parse(line,' ',args,nargs)
      if(args(2).eq.'AU') then
        do i = 1,ncent
         read (iwfn,*) atnam(i),iatom,idum,co(i,1),co(i,2),co(i,3)
         co(i,4)=dble(idum)
         if(co(i,4).lt.1.0d0) atnam(i)='xx'
         write(*,prntfrmt) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
!         write(*,203) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
        enddo
! or they are given in Angström
      else if(args(2).eq.'Angs') then
        do i = 1,ncent
         read (iwfn,*) atnam(i),iatom,idum,co(i,1),co(i,2),co(i,3)
         co(i,4)=dble(idum)
         do j=1,3
          co(i,j)=co(i,j)/0.52917721092d0
         enddo
         if(co(i,4).lt.1.0d0) atnam(i)='xx'
         write(*,prntfrmt) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
!         write(*,203) atnam(i),i,co(i,1),co(i,2),co(i,3),co(i,4)
        enddo
      endif
! reading geometry data finished! rewind input file
      rewind(iwfn)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Now read basis set data                                    c
!                                                            c
! go through atoms and assign orbitals (contracted & prims)  c
! search for [GTO] statement to read basis set data          c
! for each set of atomic orbitals                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      icount=0 !counter for primitives
      jcount=0 !counter for contractions
      call findstr('[GTO]',5,iwfn,1,1)
cccccccccccccccccccc
c go through atoms c
cccccccccccccccccccc
      do i=1,ncent
        do
  444    read(iwfn,*,iostat=iostatus)iatom,idum
         if(iostatus.ne.0) goto 444 ! if it does not fit, read next line
!  check whether there is a bf declaration in next line
         read(iwfn,'(A)',IOSTAT=iostatus) dummy
         if(iostatus.ne.0) stop 'unknown format'
! parse line
         call parse(dummy,' ',args,nargs)
         if(nargs.ne.3) cycle
           aang=args(1)
           if(aang.eq.'s'.or.aang.eq.'p'.or.aang.eq.'d'
     .       .or.aang.eq.'f'.or.trim(args(1)).eq.'sp') then
             backspace(iwfn) ! if in correct line, exit and start reading basis on atom
             exit
           endif
        enddo

      if(ncent.lt.iatom) stop 'Atom number does not fit to basis'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  check whether s,p,d,f orbital is present (aang) and read no of prims (ibas) c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do
           read(iwfn,'(A)',IOSTAT=iostatus) dummy2
           if(iostatus.ne.0) exit
! use string module
           call parse(dummy2,' ',args,nargs)
           if(nargs.lt.2)exit
           aang=args(1)
           call value(args(2),ibas,iostatus)
           if(iostatus.ne.0) stop 'ibas:error in arg to integer conv'

! special 'sp' case
           if(trim(args(1)).eq.'sp') then

! set mult. factor for higher ang. momentum functions
! we assume Cartesian coordinates - spherical ones are not implemented
! read and assign 1st component (s,px,dxx,fxxx)
             jcount=jcount+1
             ifstart(jcount)=0
             do j=1,ibas ! read data for each primitive
              icount=icount+1
              read(iwfn,*) exip(icount),cxip(icount),cxip(icount+ibas) !     exip - exponents of primitives, first s part, then p part
              exip(icount+ibas)=exip(icount)           ! s and p have same exponent in 'sp'

! s part
              ee=2.0d0*exip(icount)
              facs = pi32 / (ee * dsqrt(ee))
              cxip(icount)=cxip(icount)/dsqrt(facs)
! p part
              fac = pt5   * facs / ee
              cxip(icount+ibas)=cxip(icount+ibas)/dsqrt(fac)

              ipat(icount)=iatom
              ipat(icount+ibas)=iatom
!     ipat - primitive to atom
              ipty(icount)=1
              ipty(icount+ibas)=2
!     ipty - angular momemtum type of primitive
              ipao(icount)=jcount
              ipao(icount+ibas)=jcount+1
!     ipao - primitive to contracted
             enddo
             kbasf=3
             iang=1
             aang='p'
             mbasf=2
             icount=icount+ibas ! increase # prims
             jcount=jcount+1 ! increase # AOs
             ifstart(jcount)=0

           else
! general case
! set mult. factor for higher ang. momentum functions
! we assume Cartesian coordinates - spherical ones are not implemented
! aang:  character of function type, i.e., s,p,d,f
! mbasf: interger of ang momentum increased by one; i.e., s=1,p=2,d=3,...
! kbasf: # of Cartesian functions of this type: s=1,p=3,d=6,f=10
! lang: interger of ang momentum i.e., s=0,p=1,d=2,...
             call aangchk(1,aang,mbasf,kbasf,lang)
! read and assign 1st component (s,px,dxx,fxxx)
             jcount=jcount+1
             ifstart(jcount)=0
             do j=1,ibas ! read data for each primitive
              icount=icount+1
              read(iwfn,*) exip(icount),cxip(icount) !     exip - exponents of primitives
              ipat(icount)=iatom  !     ipat - primitive to atom
              ipty(icount)=mbasf  !     ipty - angular momemtum type of primitive
              ipao(icount)=jcount !     ipao - primitive to contracted
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! for TM & Molpro: contract exponent of primitive into contraction coefficient             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              if(mform.gt.0) then
                ee=2.0d0*exip(icount)
                facs = pi32 / (ee * dsqrt(ee))
                if (lang.eq.0) fac = facs
                if (lang.eq.1) fac = pt5   * facs / ee
                if(imethod*nbf.gt.nmo) then ! spherical MO basis
                  if (lang.eq.2) fac = pt75  * facs / (ee*ee*3.0d0)
                  if (lang.eq.3) fac = pt187 * facs / (ee*ee*ee*15.0d0)
                else if(imethod*nbf.eq.nmo) then !Cartesian MO basis
                  if (lang.eq.2) fac = pt75  * facs / (ee*ee)
                  if (lang.eq.3) fac = pt187 * facs / (ee*ee*ee)
                endif
                cxip(icount)=cxip(icount)/dsqrt(fac)
              endif
! end: contract exp. into coeff.
             enddo

          endif ! endif belongs to check if 'sp' ist present or not
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! if l > 0, assign remaining components (py,pz,dxy,...) accordingly  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if(imethod*nbf.gt.nmo) then ! for spherical mos
            do j=1,kbasf-1
              jcount=jcount+1
              mbasf=mbasf+1
              ifstart(jcount)=0
              if(mbasf.eq.19) ifstart(jcount)=1 ! set starter to rearrange f functions in MO coefficient readout part
              do k=1,ibas
               icount=icount+1
               exip(icount)=exip(icount-ibas)
               cxip(icount)=cxip(icount-ibas)
               ipat(icount)=iatom
               ipty(icount)=mbasf
               ipao(icount)=jcount
              enddo ! go through primitives
            enddo ! set other components for l>0
          else if(imethod*nbf.eq.nmo) then ! for cartesian mos
            select case (lang)
            case(2) !
              do j=1,kbasf-1
                jcount=jcount+1
                mbasf=mbasf+1
                ifstart(jcount)=0
                if(mbasf.eq.19) ifstart(jcount)=1 ! set starter to rearrange f functions in MO coefficient readout part
                do k=1,ibas
                 icount=icount+1
                 exip(icount)=exip(icount-ibas)
                 if(j.eq.3) then
                   cxip(icount)=cxip(icount-ibas)*dsqrt(3.0d0)
                 else
                   cxip(icount)=cxip(icount-ibas)
                 endif
                 ipat(icount)=iatom
                 ipty(icount)=mbasf
                 ipao(icount)=jcount
                enddo ! go through primitives
              enddo ! set other components for l>0

            case(3)
              do j=1,kbasf-1
                jcount=jcount+1
                mbasf=mbasf+1
                ifstart(jcount)=0
                if(mbasf.eq.19) ifstart(jcount)=1 ! set starter to rearrange f functions in MO coefficient readout part
                do k=1,ibas
                 icount=icount+1
                 exip(icount)=exip(icount-ibas)
                 if(j.eq.3) then
                   cxip(icount)=cxip(icount-ibas)*dsqrt(5.0d0)
                 else if(j.eq.9) then
                   cxip(icount)=cxip(icount-ibas)*dsqrt(3.0d0)
                 else
                   cxip(icount)=cxip(icount-ibas)
                 endif
                 ipat(icount)=iatom
                 ipty(icount)=mbasf
                 ipao(icount)=jcount
                enddo ! go through primitives
              enddo ! set other components for l>0

            case default
              do j=1,kbasf-1
                jcount=jcount+1
                mbasf=mbasf+1
                ifstart(jcount)=0
                if(mbasf.eq.19) ifstart(jcount)=1 ! set starter to rearrange f functions in MO coefficient readout part
                do k=1,ibas
                 icount=icount+1
                 exip(icount)=exip(icount-ibas)
                 cxip(icount)=cxip(icount-ibas)
                 ipat(icount)=iatom
                 ipty(icount)=mbasf
                 ipao(icount)=jcount
                enddo ! go through primitives
              enddo ! set other components for l>0

            end select

          endif

         enddo ! read AO on atom
         backspace(iwfn) ! be safe not to go to far
      enddo ! read/go through atom

      if(jcount.ne.nbf) then
        write(0,*) jcount,'vs',nbf
        write(*,*) jcount,'vs',nbf
        stop 'int. error in nbf read!'
      endif

      if(icount.ne.nprims) stop 'int. error in nprims read!'


      write(*,*) ' '
      select case(mform)
      case(0)
      write(*,*) 'interpreted GTO in ORCA/xTB style'
      case(1)
      if(imethod*nbf.gt.nmo) then
        write(*,*) 'interpreted GTO in TURBOMOLE style'
      else
        write(*,*) 'interpreted GTO in TURBOMOLE/MOLPRO/GAUSSIAN style'
      endif
      case(2)
        write(*,*) 'interpreted GTO in MOLPRO style'
      case(3)
        write(*,*) 'interpreted GTO in TERACHEM/GAUSSIAN style'
        write(*,'(a)',advance='no') '  ...normalize GTO contractions...'
        call normalize(.true.,nprims,ipao,ipty,exip,cxip)
        write(*,'(a)',advance='yes') 'done!'
      end select


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! read occupation number of orbitals and orbital energies c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      rewind(iwfn)
      call findstr('[MO]',4,iwfn,1,1)
! icount is now number of mos (by counting energies)
! jcount is now number of mos (by counting occupations)
      icount=0
      jcount=0
! kcount is the coefficient counter
      kcount=0
      do
       read(iwfn,'(A)',IOSTAT=iostatus) line
       if(iostatus.ne.0) exit
       call parse(line,' ',args,nargs)
        if(nargs.ne.2) cycle
        if(args(1).eq.'Ene=') then ! if energy is found, read it and search for occupation
        icount=icount+1
        call value(args(2),eps(icount),iostatus)
        if(iostatus.ne.0) stop 'eps: error in arg to real*8 conv'

! optional in uhf case: read spin
         if(imethod.eq.2) then
         read(iwfn,'(A)',IOSTAT=iostatus) dummy2
           if(iostatus.ne.0) exit
           dummy2=lowercase(dummy2)
           call removesp(dummy2)
           call parse(dummy2,'=',args,nargs)
           if(nargs.ne.2) cycle
           if(index(dummy2,'spin').ne.0) then
           if(index(dummy2,'alpha').ne.0) then
            ccspin(icount) = 1
           else
            ccspin(icount) = 2
           endif
           endif
         endif
! end: spin read

        do
         read(iwfn,'(A)',IOSTAT=iostatus) dummy
         if(iostatus.ne.0) exit
         call parse(dummy,' ',args,nargs)
         if(nargs.ne.2) cycle
         if(args(1).eq.'Occup=') then
         jcount=jcount+1
         call value(args(2),occ(jcount),iostatus)
         if(iostatus.ne.0) stop 'occ: error in arg to real*8 conv'


! read LCAO-MO coefficients
          jdum=0
          do
           read(iwfn,'(A)',IOSTAT=iostatus) dummy2
           if(iostatus.ne.0) exit
           call parse(dummy2,' ',args,nargs)
           if(nargs.ne.2) exit
           call value(args(1),idum,iostatus)
           if(iostatus.ne.0) exit
           if(idum.le.0.or.idum.gt.nbf)exit
           kdum=idum-jdum

           do i=1,kdum
             kcount=kcount+1
             cc(kcount)=0.0d0
           if(i.eq.kdum) call value(args(2),cc(kcount),iostatus)
           if(iostatus.ne.0) stop 'MO-coef: error in arg to real*8 conv'

! rearrange MO coefficients of f functions that are different
            if(ifstart(jdum+i).eq.1) then
             fyyz=cc(kcount)
             fyzz=cc(kcount-1)
             fxzz=cc(kcount-2)
             fxxz=cc(kcount-3)
             fxxy=cc(kcount-4)
             fxyy=cc(kcount-5)
             cc(kcount)=fyzz
             cc(kcount-1)=fxzz
             cc(kcount-2)=fyyz
             cc(kcount-3)=fxyy
             cc(kcount-4)=fxxz
             cc(kcount-5)=fxxy
           endif
! end of rearranging

           enddo
           jdum=idum
          enddo ! this goes back to occ loop
          backspace(iwfn)

! e.g. in Terachem, only nonzero MO coefficients are printed - fill up remaining coefficients
          if(jdum.lt.nbf) then
          kdum=nbf-jdum
          do i=1,kdum
           kcount=kcount+1
           cc(kcount)=0.0d0
           if(ifstart(jdum+i).eq.1) then ! rearrange MO coefficients of f functions that are different
            fyyz=cc(kcount)
            fyzz=cc(kcount-1)
            fxzz=cc(kcount-2)
            fxxz=cc(kcount-3)
            fxxy=cc(kcount-4)
            fxyy=cc(kcount-5)
            cc(kcount)=fyzz
            cc(kcount-1)=fxzz
            cc(kcount-2)=fyyz
            cc(kcount-3)=fxyy
            cc(kcount-4)=fxxz
            cc(kcount-5)=fxxy
            endif
! end of rearranging

          enddo
          endif
! end fill-up

          backspace(iwfn) ! backspace,because some molden inputs give no symmetry (then energy is overridden)
          exit ! go back to epsilon loop if coeffs have been read
         endif
        enddo

        endif
      enddo

cccccccccccccccccccccccccccc
! postprocess coefficients c (factor introduced for d and f in contraction coefficients needs to be multiplied/divided here)
cccccccccccccccccccccccccccc
! Turbomole style
       if(imethod*nbf.gt.nmo) then ! only if MOs generated from spherical basis
         if(mform.le.1) then
          kcount=0
          do j=1,nmo
          do i=1,nprims+1
           iprimao=ipao(i)
           iprimtyp=ipty(i)
           if(i.gt.1)then
            if(i.eq.nprims+1) iprimao=0
            if(iprimao.ne.jprimao) then
              kcount=kcount+1
              if(jprimtyp.eq.8.or.jprimtyp.eq.9.or.jprimtyp.eq.10) then
               cc(kcount)=cc(kcount)*dsqrt(3.0d0)
              endif
              if(jprimtyp.gt.13) then
               cc(kcount)=cc(kcount)*dsqrt(5.0d0)
              endif
              if(jprimtyp.eq.20) then
               cc(kcount)=cc(kcount)*dsqrt(3.0d0)
              endif
            endif
           endif
           jprimao=iprimao
           jprimtyp=iprimtyp
          enddo
          enddo
! Molpro style
         else if(mform.eq.2) then
          kcount=0
          do j=1,nmo
          do i=1,nprims+1
           iprimao=ipao(i)
           iprimtyp=ipty(i)
           if(i.gt.1)then
            if(i.eq.nprims+1) iprimao=0
            if(iprimao.ne.jprimao) then
              kcount=kcount+1
              if(jprimtyp.eq.5.or.jprimtyp.eq.6.or.jprimtyp.eq.7) then
               cc(kcount)=cc(kcount)/dsqrt(3.0d0)
              endif
              if(jprimtyp.gt.10.and.jprimtyp.lt.14) then
               cc(kcount)=cc(kcount)/dsqrt(5.0d0)
              endif
              if(jprimtyp.gt.10.and.jprimtyp.lt.20) then
               cc(kcount)=cc(kcount)/dsqrt(3.0d0)
              endif
            endif
           endif
           jprimao=iprimao
           jprimtyp=iprimtyp
          enddo
          enddo
         endif
      endif

      rewind(iwfn)
      close(iwfn)

!      call header('Orbitals',0)
!      write(*,'(/,A,/)') '  Occupancy,  Energy (eV), Orbital Spin'
!      do i = 1, nmo
!       write(*,'(F8.2,F12.4,I4)') occ(i),eps(i)*27.21139
!,ccspin(i)
!      enddo

!      endif

      write(*,95) ncent,nmo,nprims,nbf
 95   format (/,1x,'# atoms          =',i5,/,
     .          1x,'# mos            =',i5,/,
     .          1x,'# primitive  aos =',i5,/,
     .          1x,'# contracted aos =',i5,/)

      if(imethod*nbf.gt.nmo)then
         write(*,*) 'spherical AO basis'
         if(mform==1.or.mform==2)then
         spherical=.true.
         else
         spherical=.false.
         endif
      else
         write(*,*) 'cartesian AO basis'
         spherical=.false.
      endif

      call etafill(nprims)

!203   format(2x,a2,i3,2x,3f14.8,3x,f10.2)

      end

! search for string in input
      subroutine findstr(str,lenstr,ifile,i,n)
      use strings
      implicit none
      integer, intent( in ) :: lenstr,i,n,ifile
      character(len=lenstr), intent( in ) :: str
      character*25 arg(10)
      integer narg,ios
      character*79 line
      if(n.le.0) then
       do
        read(ifile,'(A)',IOSTAT=ios) line
        if(ios.lt.0) write(*,*) 'No ',str,' flag found'
        if(ios.lt.0) stop 'end of file reached'
        call parse(line,' ',arg,narg)
        if(arg(i).eq.str) exit
       enddo

      else
       do
        read(ifile,'(A)',IOSTAT=ios) line
        if(ios.lt.0) write(*,*) 'No ',str,' flag found'
        if(ios.lt.0) stop 'end of file reached'
        call parse(line,' ',arg,narg)
        if(narg.ne.n) cycle
        if(arg(i).eq.str) exit
       enddo
      endif
      end



      subroutine aangchk(modus,chr,mbasf,kbasf,iang)
      implicit none
      character*(*), intent( in ) :: chr
      integer, intent( in ) :: modus
      integer, intent( out ) :: mbasf,kbasf,iang
      iang=0
      if(modus.eq.0) then
         select case(chr)
          case('s')
           kbasf=1
           mbasf=1
          case('p')
           kbasf=3
           mbasf=3
          case('d')
           kbasf=6
           mbasf=5
          case('f')
           kbasf=10
           mbasf=7
          case('g')
           stop'ang. momentum > f not implemented! Exiting!'
          case('h')
           stop'ang. momentum > f not implemented! Exiting!'
          case('i')
           stop'ang. momentum > f not implemented! Exiting!'
          case default
           iang=-1
         end select
      else if(modus.eq.1) then
        select case(chr)
         case('s')
          kbasf=1 ! # of cartesian bfs
          mbasf=1 ! starter for ipty
          iang=0
         case('p')
          kbasf=3
          mbasf=2
          iang=1
         case('d')
          kbasf=6
          mbasf=5
          iang=2
         case('f')
          kbasf=10
          mbasf=11
          iang=3
         case default
          continue
        end select
      endif
      end

      subroutine lenint(iin,iout)
      implicit none
      integer, intent( in ) :: iin
      integer, intent( out) :: iout
      integer mdim,jdim
cccccccccccccccccccccccccccccccccccccccccccccc
c
c specifies the length of an integer
c
cccccccccccccccccccccccccccccccccccccccccccccc
      jdim=0
      iout=0
      do
      iout=iout+1
      jdim=10*jdim+9
      mdim=iin-jdim
      if(mdim.le.0) exit
      enddo
      return
      end


c fill up eta array
      subroutine etafill(nprims)
      use stdacommon
      implicit double precision (a-h,o-z)

      common /carte  / lmn(0:3,0:3,0:3)


c=======================================================================
c cartesian gaussian functions (6d,10f...)
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c=======================================================================

      lmn(0,0,0)=1
      lmn(1,0,0)=2
      lmn(0,1,0)=3
      lmn(0,0,1)=4
      lmn(2,0,0)=5
      lmn(0,2,0)=6
      lmn(0,0,2)=7
      lmn(1,1,0)=8
      lmn(1,0,1)=9
      lmn(0,1,1)=10
      lmn(3,0,0)=11
      lmn(0,3,0)=12
      lmn(0,0,3)=13
      lmn(2,1,0)=14
      lmn(2,0,1)=15
      lmn(1,2,0)=16
      lmn(0,2,1)=17
      lmn(1,0,2)=18
      lmn(0,1,2)=19
      lmn(1,1,1)=20

      do i=1,nprims
         iat=ipat(i)
         eta(i,1)=co(iat,1)
         eta(i,2)=co(iat,2)
         eta(i,3)=co(iat,3)
         eta(i,4)=exip(i)
         eta(i,5)=float(ipty(i))
         if(ipty(i).gt.20) then
          write(*,*) ' '
          write(*,*)'WARNING WARNING WARNING WARNING WARNING WARNING'
          write(*,*)'     Functions > f present in basis but not'
          write(*,*)'  implemented! The program will stop here. '
          write(*,*)'      Redo the SCF without functions > f!'
          write(*,*)'WARNING WARNING WARNING WARNING WARNING WARNING'
          write(0,*)'WARNING: Functions > f in basis. Program stopped!'
          stop
         endif
        do m=1,20
            eta(i,m+5)=0.0d0
         enddo
         eta(i,5+ipty(i))=1.00d0
      enddo

      return
      end

      subroutine readinp(ax,thre,alp,bet)
      use strings
      implicit none
      real*8 ax,thre,alp,bet
      integer ios,narg
      character*79 line
      character*25 arg(10)
      open(unit=89,file='.STDA')
  777 read(89,'(A)',iostat=ios) line
      if(ios.lt.0) stop '.STDA file empty'
      call parse(line,' ',arg,narg)
      if(arg(1).eq.'#') goto 777
      if(narg.ge.1) call value(arg(1),ax,ios)
      if(ios.ne.0) stop 'ax from file unreadable'
      if(narg.ge.2) call value(arg(2),thre,ios)
      if(ios.ne.0) stop 'Ethr from file unreadable'
      if(narg.ge.3) call value(arg(3),alp,ios)
      if(ios.ne.0) stop 'alpha from file unreadable'
      if(narg.ge.4) call value(arg(4),bet,ios)
      if(ios.ne.0) stop 'beta from file unreadable'
      close(89)
      end
