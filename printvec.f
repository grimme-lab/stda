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
      subroutine printvectda(rpachk,nci,nroot,uci,eci)
      ! This routine prints out eigenvectors in TM format
      integer, intent ( in ) :: nci,nroot
      real*4, intent ( in ) :: uci(nci,nci),eci(nci)
      logical, intent ( in ) :: rpachk
      integer i,j,ij,k,jhilf,ii,ihilf,jhomo,nvec,itype,nmo,na
      integer jhomoa,jhomob,nexa,nexb,khilf,lhilf
      integer, allocatable :: iconf(:,:),vecchk(:),icsf(:)
      integer, allocatable :: iconfb(:,:),vecchkb(:),icsfb(:) ! beta variables
      character*79 fname

      open(unit=39,file='TmPvEcInFo',status='old')
      ! read from temporary file
      rewind(39)
      read(39,*)itype
      if(itype.ne.2) then
! restricted case
        read(39,*) nvec,nmo,jhomo
        allocate(vecchk(nmo), iconf(nci,2),icsf(jhomo*(nmo-jhomo)),
     .         stat=i)
        if(i.ne.0)stop 'allocation failed in printz'
        do i=1,nmo
          read(39,*) vecchk(i)
        enddo
         do i=1,nci
          read(39,*) iconf(i,1),iconf(i,2)
        enddo
      else
! unrestricted case
        read(39,*) nvec,nmo,jhomo
        read(39,*) nexa,nexb,jhomoa,jhomob
        na=nmo/2 ! number of alpha (and beta) mos
        ij=jhomoa*(na-jhomoa)+jhomob*(na-jhomob)
        allocate(vecchk(na), iconf(nexa,2),iconfb(nexb,2),
     .   vecchkb(na),icsf(jhomoa*(na-jhomoa)),icsfb(jhomob*(na-jhomob))
     .   , stat=i)
        if(i.ne.0)stop 'allocation failed in printz'
         do i=1,na
          read(39,*) vecchk(i)
         enddo
        do i=1,na
          read(39,*) vecchkb(i)
        enddo
        do i=1,nexa
          read(39,*) iconf(i,1),iconf(i,2)
        enddo
        do i=1,nexb
          read(39,*) iconfb(i,1),iconfb(i,2)
        enddo
      endif

      close(39,status='delete')

      if (nvec.lt.1.or.nvec.gt.nroot) nvec=nroot
      select case(itype)
       case(0) ! RKS singlet
        fname='ciss_a'
        if(rpachk) fname='sing_a' ! UKS RPA
       case(1) ! RKS triplet
        fname='cist_a'
        if(rpachk) fname='trip_a' ! UKS RPA
       case(2) ! UKS TDA
        fname='ucis_a'
        if(rpachk) fname='unrs_a' ! UKS RPA
      end select
      open(unit=29,file=trim(fname),status='replace')
      write(29,'(a)',advance='yes') '$title'
      write(29,'(a)',advance='yes') '$symmetry c1'
      write(29,'(a)',advance='no') '$tensor space dimension'
      if(itype.eq.2) then
        na=nmo/2
        ij=jhomoa*(na-jhomoa)+jhomob*(na-jhomob)
        write(29,'(x,i8)',advance='yes') ij
      else
        write(29,'(x,i8)',advance='yes') jhomo*(nmo-jhomo)
      endif
      select case(itype)
       case(0) ! RKS singlet
        if(rpachk) then
          write(29,'(a)',advance='yes') '$scfinstab rpas'
        else
          write(29,'(a)',advance='yes') '$scfinstab ciss'
        endif
       case(1) ! RKS triplet
        if(rpachk) then
          write(29,'(a)',advance='yes') '$scfinstab rpat'
        else
          write(29,'(a)',advance='yes') '$scfinstab cist'
        endif
       case(2) ! UKS RPA
        if(rpachk) then
          write(29,'(a)',advance='yes') '$scfinstab urpa'
        else
          ! UKS TDA
          write(29,'(a)',advance='yes') '$scfinstab ucis'
        endif
      end select
      write(29,'(a)',advance='no') '$current subspace dimension'
      write(29,'(x,i8)',advance='yes') nvec
      write(29,'(a)',advance='yes') '$current iteration converged'
      write(29,'(a)',advance='no') '$eigenpairs'

! resort
! since TURBOMOLE uses order occ1->virt1, occ1->virt2,....occn->virtn
      icsf=0
      ! sort for restricted case
      if(itype.ne.2) then
        do i=1,nci ! go through CSFs
          ihilf=iconf(i,1)  ! get occ and virt from CSF (in sTDA sorting)
          jhilf=iconf(i,2)
          l=0
          do j=1,jhomo ! now go through all orbitals (in input sorting)
            if (vecchk(j).eq.ihilf) then ! see whether this orbital corresponds to the occupied one from CSF
              do k=jhomo+1,nmo
                l=l+1
                if(vecchk(k).eq.jhilf) then   ! see whether this orbital corresponds to the virtual one from CSF
                  icsf(l)=i ! now set marker
!                  l=l+nmo-k
                  exit
                endif
              enddo
!              l=l+(jhomo-j)*(nmo-jhomo)
              exit
            else
              l=l+(nmo-jhomo)
            endif
          enddo
        enddo
      else
      !sort for unrestricted case
        icsfb=0
        do i=1,nexa
          ihilf=iconf(i,1)
          jhilf=iconf(i,2)
          l=0
          ! alpha-alpha part
          do j=1,jhomoa
            if (vecchk(j).eq.ihilf) then
              do k=jhomoa+1,na
                l=l+1
                if(vecchk(k).eq.jhilf) then
                  icsf(l)=i
!                  l=l+na-k ! if found, skip remaining loops
                  exit
                endif
              enddo
!              l=l+(jhomoa-j)*(na-jhomoa)
              exit
            else
              l=l+na-jhomoa
            endif
          enddo
        enddo
          ! beta-beta part
        do i=1,nexb
          khilf=iconfb(i,1)
          lhilf=iconfb(i,2)
          l=0
          do j=1,jhomob
            if (vecchkb(j).eq.khilf) then
              do k=jhomob+1,na
                l=l+1
                if(vecchkb(k).eq.lhilf) then
                  icsfb(l)=i+nexa
!                  l=l+na-k ! if found, skip remaining loops
                  exit
                endif
              enddo
!              l=l+(jhomob-j)*(na-jhomob)
            else
              l=l+na-jhomob
            endif
          enddo
        enddo

      endif
!! done sorting

      ihilf=1
      do i = 1, nvec
         if (ihilf.le.4.and.ihilf.gt.0) then
           ihilf=0
           write(29,'(a)',advance='yes') ' '
         endif
         write(29,'(3x,i6,3x,a13,d23.16)') i, 'eigenvalue = ',
     &        eci(i)


! now print Z/TDA vecs

      if(itype.ne.2) then
      !restricted
        l=0
        do j=1,jhomo
         if(vecchk(j).eq.0) then

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop

         else

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              if(abs(uci(icsf(l),i)).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') dble(uci(icsf(l),i))
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! llop over virts
        endif
       end do ! loop over occs

      else
      ! unrestricted
        l=0
        do j=1,jhomoa
         if(vecchk(j).eq.0) then

          ! alpha
          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop-alpha

         else

          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              if(abs(uci(icsf(l),i)).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') dble(uci(icsf(l),i))
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! loop over virts-alpha
        endif
        end do ! loop over occs-alpha

       ! beta
       l=0
       do j=1,jhomob
        if(vecchkb(j).eq.0) then
         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         enddo ! virt loop-alpha

        else

         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsfb(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              if(abs(uci(icsf(l),i)).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') dble(uci(icsfb(l),i))
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         end do ! llop over virts-beta
        endif
       end do ! loop over occs-beta

      endif ! check UKS/RKS

      end do ! loop over roots/nvec



      if(ihilf.le.3.and.ihilf.gt.0) write(29,'(a)',advance='yes') ' '
      write(29,'(a)',advance='yes') '$end'
      close(29)

      deallocate(vecchk,iconf,icsf)
      if(itype.eq.2) deallocate(iconfb,vecchkb,icsfb)
! note the printout in stdout
      write(*,'(a)',advance='no') ' eigenvectors printed to '
      write(*,'(a)',advance='yes') trim(fname)

      end


      subroutine printvecrpa(nci,nroot,uci,hci,eci)
      ! This routine prints out eigenvectors in TM format
      integer, intent ( in ) :: nci,nroot
      real*4, intent ( in ) :: uci(nci,nci),hci(nci,nci),eci(nci)
      real*8 tmp
      integer i,j,ij,k,jhilf,ii,ihilf,jhomo,nvec,itype,nmo,na
      integer jhomoa,jhomob,nexa,nexb,khilf,lhilf
      integer, allocatable :: iconf(:,:),vecchk(:),icsf(:)
      integer, allocatable :: iconfb(:,:),vecchkb(:),icsfb(:) ! beta variables
      character*79 fname

      open(unit=39,file='TmPvEcInFo',status='old')
      ! read from temporary file
      rewind(39)
      read(39,*)itype
      if(itype.ne.2) then
! restricted case
        read(39,*) nvec,nmo,jhomo
        allocate(vecchk(nmo), iconf(nci,2),icsf(jhomo*(nmo-jhomo)),
     .         stat=i)
        if(i.ne.0)stop 'allocation failed in printz'
        do i=1,nmo
          read(39,*) vecchk(i)
        enddo
         do i=1,nci
          read(39,*) iconf(i,1),iconf(i,2)
        enddo
      else
! unrestricted case
        read(39,*) nvec,nmo,jhomo
        read(39,*) nexa,nexb,jhomoa,jhomob
        na=nmo/2
        ij=jhomoa*(na-jhomoa)+jhomob*(na-jhomob)
        allocate(vecchk(na), iconf(nexa,2),iconfb(nexb,2),
     .    vecchkb(na),icsf(jhomoa*(na-jhomoa)),icsfb(jhomob*(na-jhomob))
     .    , stat=i)
        if(i.ne.0)stop 'allocation failed in printz'
         do i=1,na
          read(39,*) vecchk(i)
         enddo
        do i=1,na
          read(39,*) vecchkb(i)
        enddo
        do i=1,nexa
          read(39,*) iconf(i,1),iconf(i,2)
        enddo
        do i=1,nexb
          read(39,*) iconfb(i,1),iconfb(i,2)
        enddo
      endif

      close(39,status='delete')

      if (nvec.lt.1.or.nvec.gt.nroot) nvec=nroot
      select case(itype)
       case(0) ! RKS singlet
        fname='sing_a' ! UKS RPA
       case(1) ! RKS triplet
        fname='trip_a' ! UKS RPA
       case(2) ! UKS TDA
        fname='unrs_a' ! UKS RPA
      end select
      open(unit=29,file=trim(fname),status='replace')
      write(29,'(a)',advance='yes') '$title'
      write(29,'(a)',advance='yes') '$symmetry c1'
      write(29,'(a)',advance='no') '$tensor space dimension'
      if(itype.eq.2) then
        na=nmo/2
        ij=jhomoa*(na-jhomoa)+jhomob*(na-jhomob)
        write(29,'(x,i8)',advance='yes') ij
      else
        write(29,'(x,i8)',advance='yes') jhomo*(nmo-jhomo)
      endif
      select case(itype)
       case(0) ! RKS singlet
          write(29,'(a)',advance='yes') '$scfinstab rpas'
       case(1) ! RKS triplet
          write(29,'(a)',advance='yes') '$scfinstab rpat'
       case(2) ! UKS RPA
          write(29,'(a)',advance='yes') '$scfinstab urpa'
      end select
      write(29,'(a)',advance='no') '$current subspace dimension'
      write(29,'(x,i8)',advance='yes') nvec
      write(29,'(a)',advance='yes') '$current iteration converged'
      write(29,'(a)',advance='no') '$eigenpairs'

! resort
! since TURBOMOLE uses order occ1->virt1, occ1->virt2,....occn->virtn
      icsf=0
      ! sort for restricted case
      if(itype.ne.2) then
        do i=1,nci
          ihilf=iconf(i,1)
          jhilf=iconf(i,2)
          l=0
          do j=1,jhomo
            if (vecchk(j).eq.ihilf) then
              do k=jhomo+1,nmo
                l=l+1
                if(vecchk(k).eq.jhilf) then
                  icsf(l)=i
!                  l=l+nmo-k
                  exit
                endif
              enddo
!              l=l+(jhomo-j)*(nmo-jhomo)
              exit
            else
              l=l+nmo-jhomo
            endif
          enddo
        enddo
      else
      !sort for unrestricted case
        icsfb=0
        do i=1,nexa
          ihilf=iconf(i,1)
          jhilf=iconf(i,2)
          l=0
          ! alpha-alpha part
          do j=1,jhomoa
            if (vecchk(j).eq.ihilf) then
              do k=jhomoa+1,na
                l=l+1
                if(vecchk(k).eq.jhilf) then
                  icsf(l)=i
!                  l=l+na-k ! if found, skip remaining loops
                  exit
                endif
              enddo
!              l=l+(jhomoa-j)*(na-jhomoa)
              exit
            else
              l=l+na-jhomoa
            endif
          enddo
        enddo
        ! beta-beta part
        do i=1,nexb
          khilf=iconfb(i,1)
          lhilf=iconfb(i,2)
          l=0
          do j=1,jhomob
            if (vecchkb(j).eq.khilf) then
              do k=jhomob+1,na
                l=l+1
                if(vecchkb(k).eq.lhilf) then
                  icsfb(l)=i+nexa
!                  l=l+na-k ! if found, skip remaining loops
                  exit
                endif
              enddo
!              l=l+(jhomob-j)*(na-jhomob)
            else
              l=l+na-jhomob
            endif
          enddo
        enddo

      endif
!! done sorting

      ihilf=1
      do i = 1, nvec
         if (ihilf.le.4.and.ihilf.gt.0) then
           ihilf=0
           write(29,'(a)',advance='yes') ' '
         endif
         write(29,'(3x,i6,3x,a13,d23.16)') i, 'eigenvalue = ',
     &        eci(i)
! print RPA vec

      if(itype.ne.2) then
      !restricted
      !X+Y
       l=0
       do j=1,jhomo
         if(vecchk(j).eq.0) then

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop

         else

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsf(l),i)+hci(icsf(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! llop over virts
         endif
       end do ! loop over occs

       !X-Y
       l=0
       do j=1,jhomo
         if(vecchk(j).eq.0) then

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop

         else

          do k=jhomo+1,nmo
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsf(l),i)-hci(icsf(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! llop over virts
         endif
       end do ! loop over occs
      else
      ! unrestricted
      ! X+Y
        l=0
        do j=1,jhomoa
         if(vecchk(j).eq.0) then

          ! alpha
          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop-alpha

         else

          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsf(l),i)+hci(icsf(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! llop over virts-alpha
        endif
        end do ! loop over occs-alpha

       ! beta
       l=0
       do j=1,jhomob
        if(vecchkb(j).eq.0) then
         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         enddo ! virt loop-alpha

        else

         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsfb(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsfb(l),i)+hci(icsfb(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         end do ! llop over virts-beta
        endif
       end do ! loop over occs-beta
       ! now X - Y , unrestricted
       l=0
       do j=1,jhomoa
         if(vecchk(j).eq.0) then

          ! alpha
          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          enddo ! virt loop-alpha

         else

          do k=jhomoa+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsf(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsf(l),i)-hci(icsf(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
          end do ! loop over virts-alpha
        endif
       end do ! loop over occs-alpha

       ! beta
       l=0
       do j=1,jhomob
        if(vecchkb(j).eq.0) then
         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            write(29,'(d20.14)',advance='no') 0.0d0
            if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         enddo ! virt loop-alpha

        else

         do k=jhomob+1,na
            l=l+1
            ihilf = ihilf+1
            if(icsfb(l).eq.0) then
              write(29,'(d20.14)',advance='no') 0.0d0
            else
              tmp=dble(uci(icsfb(l),i)-hci(icsfb(l),i))
              if(abs(tmp).lt.1.d-99) then
               write(29,'(d20.14)',advance='no') 0.0d0
              else
               write(29,'(d20.14)',advance='no') tmp
              endif
            endif
              if(ihilf.gt.3) then
              ihilf=0
              write(29,'(a)',advance='yes') ' '
            endif
         end do ! llop over virts-beta
        endif
       end do ! loop over occs-beta

      endif ! check UKS/RKS

      end do ! loop over roots/nvec



      if(ihilf.le.3.and.ihilf.gt.0) write(29,'(a)',advance='yes') ' '
      write(29,'(a)',advance='yes') '$end'
      close(29)

      deallocate(vecchk,iconf,icsf)
      if(itype.eq.2) deallocate(iconfb,icsfb,vecchkb)
! note the printout in stdout
      write(*,'(a)',advance='no') ' eigenvectors printed to '
      write(*,'(a)',advance='yes') trim(fname)

      end
