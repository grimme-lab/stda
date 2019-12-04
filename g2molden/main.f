! converts g09 output into a molden input that is readable by stda program 
      program g2molden  
      use strings
      implicit none
! xyz : cartesian coordinates
      real*8, allocatable :: xyz(:,:),coeff(:,:)
! attyp : atom type
      character*2, allocatable :: attyp(:)
! atnum : nuclear charge of atom
      integer, allocatable :: atnum(:)
! i,j,k,l,jatom some variables 
! ncent : # of atoms
! ifile : unit number of file
      integer i,j,k,l,ncent,jatom,ifile,ios,nargs,nbf,nmo
      integer nal,nbe, irun,jrem,iocc(5),ityp(5),isp
      real*8 tmp,eval(5)
      character*79 afile,dummy,frmt
      character*25 args(10),atmp(10)
      logical ex,uhf,nosym

      nmo=0
      uhf=.false.
      ncent=0
      nosym=.false.

! get file name
      call getarg(1,afile)
! file exists?
      inquire(file=trim(afile),exist=ex)
      if(.not.ex)then
          write(*,*)'file:',afile,' not found'
          stop 'input file not found'
      endif
! open file
      ifile=42
      open(unit=ifile,file=afile,status='old')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GEOMETRY DATA 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get number of atoms
      call findstr('NAtoms=',ifile,1,0)
      backspace(ifile)     
      read(ifile,'(A)',IOSTAT=ios) dummy
      if(ios.lt.0) stop 'error when reading NAtoms'
! use string module to read # atoms from line
      call parse(dummy,' ',args,nargs)
      if(nargs.lt.2) stop 'cannot find # of atoms'
      call value(args(2),ncent,ios)
      if(ios.lt.0) stop 'error - NAtoms not an integer'

! check for "NoSymm" flag
      rewind(ifile)     
      do
        read(ifile,'(A)',IOSTAT=ios) dummy
        dummy=lowercase(dummy) 
        if(index(dummy,'nosym').ne.0.or. 
     .     index(dummy,'symmetry=none').ne.0) then 
         nosym=.true.
         exit
        endif
        if(index(dummy,'natoms=').ne.0) exit
      enddo
      rewind(ifile) 

! now move to atom coordinates & read them
      if(nosym) then 
! nosym case
        call findstr('Input',ifile,1,2)
        backspace(ifile)
        read(ifile,'(A)',IOSTAT=ios) dummy
        if(ios.lt.0)stop 'reached end of file while searching for atoms'
        call parse(dummy,' ',args,nargs)
        if(nargs.ne.2) stop 'cannot find flag: "Standard orientation:"'
        if(trim(args(1)).ne.'Input') stop 'cannot find atoms'
        if(trim(args(2)).ne.'orientation:') stop 'cannot find atoms'
        read(ifile,*)
        read(ifile,'(A)',IOSTAT=ios) dummy
        call parse(dummy,' ',args,nargs)
        if(trim(args(1)).ne.'Center') stop 'cannot find atoms'
        read(ifile,'(A)',IOSTAT=ios) dummy
        call parse(dummy,' ',args,nargs)
        if(trim(args(1)).ne.'Number') stop 'cannot find atoms'
        read(ifile,*)
      else
! sym case
        call findstr('Standard',ifile,1,2)
        backspace(ifile)
        read(ifile,'(A)',IOSTAT=ios) dummy
        if(ios.lt.0)stop 'reached end of file while searching for atoms'
        call parse(dummy,' ',args,nargs)
        if(nargs.ne.2) stop 'cannot find flag: "Standard orientation:"'
        if(trim(args(1)).ne.'Standard') stop 'cannot find atoms'
        if(trim(args(2)).ne.'orientation:') stop 'cannot find atoms'
        read(ifile,*)
        read(ifile,'(A)',IOSTAT=ios) dummy
        call parse(dummy,' ',args,nargs)
        if(trim(args(1)).ne.'Center') stop 'cannot find atoms'
        read(ifile,'(A)',IOSTAT=ios) dummy
        call parse(dummy,' ',args,nargs)
        if(trim(args(1)).ne.'Number') stop 'cannot find atoms'
        read(ifile,*)
      endif

! allocate the atomic geometry data
      allocate( xyz(3,ncent), attyp(ncent), atnum(ncent), stat=ios)  
! now do the actual read-out
      do i=1,ncent
        read(ifile,*)j,atnum(i),l,(xyz(k,i),k=1,3)
        if (j.ne.i) stop 'error while reading atoms'
        call getelement(atnum(i),attyp(i))
      enddo

      write(*,'(A)',advance='yes') '[Molden Format]'
      write(*,'(A)',advance='yes') '[Atoms] Angs'
      do i=1,ncent
        write(*,100)attyp(i),i,atnum(i),(xyz(k,i),k=1,3)
      enddo
  100 format(a2,i5,i5,x,3f20.10) ! geometry printout
      deallocate(xyz,attyp,atnum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gaussian AO basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data is read in atom by atom and simply printed out as is
      call findstr('normalization):',ifile,12,12)
      write(*,'(A)',advance='yes') '[GTO]'
      do i=1,ncent 
       do 
         read(ifile,'(A)',IOSTAT=ios) dummy
         if(ios.lt.0) stop 'reached end of file during GTO readout'
         if(trim(dummy).eq.' ****') then 
           write(*,*) 
           exit  
         endif 
         call parse(lowercase(dummy),' ',args,nargs) 
         j=min(nargs,3)
         write(*,*) (' ',trim(args(k)),k=1,j)
       enddo 
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check whether UHF or RHF is present  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! assume restricted MOs and search for UHF flag
      do
       read(ifile,'(A)',IOSTAT=ios) dummy
       if(ios.lt.0) exit
       call parse(dummy,' ',args,nargs)
       if(nargs.ne.4) cycle
       if(trim(args(1)).ne.'UHF') cycle
       if(trim(args(2)).ne.'open') cycle
       if(trim(args(3)).ne.'shell') cycle
       if(trim(args(4)).ne.'SCF:') cycle
       uhf=.true.
      enddo      
      rewind(ifile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MO coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first find # of basis functions
      call findstr('functions,',ifile,3,10)
      backspace(ifile)
      read(ifile,'(A)',IOSTAT=ios) dummy
      if(ios.lt.0) stop 'End of file reached when reading # bf'
! use string module to read # of basis functions
      call parse(dummy,' ',args,nargs)
      if(nargs.lt.10) stop 'internal error # AO readout'
      call value(args(1),i,ios)
      if(ios.ne.0) stop '# AOs not an integer?'
      call value(args(7),j,ios)
      if(ios.ne.0) stop '# Cart. AOs not an integer?'
      if(i.ne.j) stop '# of AOs and MOs is different... exit here' 
! number of AOs = i
      nbf=i

! determine number of alpha and beta electrons
      read(ifile,'(A)',IOSTAT=ios) dummy
      if(ios.lt.0) stop 'End of file reached when reading # electrons'
      call parse(dummy,' ',args,nargs)
      call value(args(1),nal,ios)
      if(ios.ne.0) stop '# alpha electrons not an integer?'
      call value(args(4),nbe,ios)
      if(ios.ne.0) stop '# alpha electrons not an integer?'
      
! determine number of MOs (even if Cartesian, Gaussian sometimes removes MOs due to linear dependencies)
      call findstr('NBsUse=',ifile,1,0)
      backspace(ifile)
      read(ifile,'(A)',IOSTAT=ios) dummy
      if(ios.lt.0) stop 'End of file reached when reading # mos'
! use string module to read # of basis functions
      call parse(dummy,' ',args,nargs)
      if(nargs.lt.2) stop 'internal error # MO readout'
      call value(args(2),nmo,ios)
      if(ios.ne.0) nmo=nbf ! assume nmo=nbf


      allocate( coeff(nbf,5), stat=ios) 
      if(ios.ne.0) stop 'error allocating coeff'

      write(*,'(A)',advance='yes')'[MO]'

      if(.not.uhf) then 

        call findstr('Coefficients:',ifile,3,3) 
        backspace(ifile) 
        read(ifile,'(A)',IOSTAT=ios) dummy
        if(ios.lt.0) stop 'End of file reached searching for MOs'
        call parse(dummy,' ',args,nargs)
        if(trim(args(1)).ne.'Molecular') stop 'cannot find MOs'
        if(trim(args(2)).ne.'Orbital') stop 'cannot find MOs' 

! check how many blocks need to be read out
        tmp=dble(nmo)/5.0d0
        irun=int(tmp)
        jrem=nmo-5*irun

! now read a block of 5 columns each 
        do i=1,irun

! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,5)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,5
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=2 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,5 
              eval(j)=0.0d0
              call value(args(nargs-5+j),eval(j),ios)
              if(ios.ne.0) then  
                ! unformatted way of reading data does not work: switch to formatted readout
                backspace(ifile)
                read(ifile,101) (eval(k),k=1,5)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigval trafo'
                write(*,*) dummy
                write(*,*) args(nargs-5+k)
                stop 'error in character to eigval trafo'
              endif
           enddo

! now read coefficients
 101  format(21x,5f10.5)
           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,5
                coeff(j,k)=0.0d0
                call value(args(nargs-5+k),coeff(j,k),ios)
                if(ios.ne.0) then
                  ! fallback: formatted
                  backspace(ifile)
                  read(ifile,101) (coeff(j,l),l=1,5)
                  if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest
                  write(*,*) 'error in character to coeff1 trafo'
                  write(*,*) dummy
                  write(*,*) args(nargs-5+k)
                  stop 'error in character to coeff1 trafo'
                endif
              enddo
           enddo            

! now print out a set of 5 MOs

           do j=1,5
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Alpha'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,'(f12.8)')coeff(k,j)
            enddo

           enddo
 
        enddo  


       write(frmt,'(a,i1,a)') '(21x,',jrem,'f10.5)'
! same procedure for the last, incomplete block
       if(jrem.gt.0) then
         ! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,jrem)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,jrem
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=2 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,jrem
              eval(j)=0.0d0
              call value(args(nargs-jrem+j),eval(j),ios)
              if(ios.ne.0) then
                ! fallback: formatted 
                backspace(ifile)
                read(ifile,frmt) (eval(k),k=1,jrem)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigval trafo'
                write(*,*) dummy
                write(*,*) args(nargs-jrem+j)
                stop 'error in character to eigval trafo'
              endif
           enddo

! now read coefficients

           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,jrem
                coeff(j,k)=0.0d0
                call value(args(nargs-jrem+k),coeff(j,k),ios)
                if(ios.ne.0) then
                 ! fallback: formatted
                 backspace(ifile)
                 read(ifile,frmt) (coeff(j,l),l=1,jrem)
                 if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest
                 write(*,*) 'error in character to coeff1 trafo'
                 write(*,*) dummy
                 write(*,*) args(nargs-jrem+j)
                 stop 'error in character to coeff1 trafo'
                endif
              enddo
           enddo

! now print out the set of remaining MOs

           do j=1,jrem
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Alpha'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,*)coeff(k,j)
            enddo

           enddo
      
       endif



      else 

! UHF case
       do isp=1,2
        call findstr('Coefficients:',ifile,4,4)
        backspace(ifile)
        read(ifile,'(A)',IOSTAT=ios) dummy
        if(ios.lt.0) stop 'End of file reached searching for MOs'
        call parse(dummy,' ',args,nargs)
        if(trim(args(2)).ne.'Molecular') stop 'cannot find MOs'
        if(trim(args(3)).ne.'Orbital') stop 'cannot find MOs'
        if(trim(args(1)).eq.'Alpha') then

! check how many blocks need to be read out
        tmp=dble(nmo)/5.0d0
        irun=int(tmp)
        jrem=nmo-5*irun

! now read a block of 5 columns each 
        do i=1,irun

! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,5)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,5
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=1 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,5 
              eval(j)=0.0d0
              call value(args(nargs-5+j),eval(j),ios)
              if(ios.ne.0) then
                backspace(ifile)
                read(ifile,101) (eval(k),k=1,5)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigvala trafo'
                write(*,*) dummy
                write(*,*) args(nargs-5+k)
                stop 'error in character to eigvala trafo'
              endif
           enddo
           
! now read coefficients

           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,5
                coeff(j,k)=0.0d0
                call value(args(nargs-5+k),coeff(j,k),ios)
                if(ios.ne.0) then 
                  ! fallback: formatted
                  backspace(ifile)
                  read(ifile,101) (coeff(j,l),l=1,5)
                  if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest
                  write(*,*) 'error in character to coeffa trafo'
                  write(*,*) dummy
                  write(*,*) args(nargs-5+k)
                  stop 'error in character to coeffa trafo'
                endif
              enddo
           enddo            

! now print out a set of 5 MOs

           do j=1,5
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Alpha'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,*)coeff(k,j)
            enddo

           enddo
 
        enddo  

! same procedure for the last, incomplete block
       write(frmt,'(a,i1,a)') '(21x,',jrem,'f10.5)'
       if(jrem.gt.0) then
         ! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,jrem)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,jrem
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=1 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,jrem
              eval(j)=0.0d0
              call value(args(nargs-jrem+j),eval(j),ios)
              if(ios.ne.0) then
                ! fallback: formatted 
                backspace(ifile)
                read(ifile,frmt) (eval(k),k=1,jrem)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigvala trafo'
                write(*,*) dummy
                write(*,*) args(nargs-jrem+j)
                stop 'error in character to eigvala trafo'
              endif
           enddo

! now read coefficients

           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,jrem
                coeff(j,k)=0.0d0
                call value(args(nargs-jrem+k),coeff(j,k),ios)
                if(ios.ne.0) then
                 ! fallback: formatted
                 backspace(ifile)
                 read(ifile,frmt) (coeff(j,l),l=1,jrem)
                 if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest                 
                 write(*,*) 'error in character to coeffa trafo'
                 write(*,*) dummy
                 write(*,*) args(nargs-jrem+k)
                 stop 'error in character to coeffa trafo'
                endif
              enddo
           enddo

! now print out the set of remaining MOs

           do j=1,jrem
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Alpha'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,*)coeff(k,j)
            enddo

           enddo
      
        endif
         
        else if(trim(args(1)).eq.'Beta') then

! check how many blocks need to be read out
        tmp=dble(nmo)/5.0d0
        irun=int(tmp)
        jrem=nmo-5*irun

! now read a block of 5 columns each 
        do i=1,irun

! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,5)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,5
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=1 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,5 
              eval(j)=0.0d0
              call value(args(nargs-5+j),eval(j),ios)
              if(ios.ne.0) then
                ! fallback: formatted 
                backspace(ifile)
                read(ifile,101) (eval(k),k=1,5)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigvalb trafo'
                write(*,*) dummy
                write(*,*) args(nargs-5+j)
                stop 'error in character to eigvalb trafo'
              endif
           enddo
           
! now read coefficients

           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,5
                coeff(j,k)=0.0d0
                call value(args(nargs-5+k),coeff(j,k),ios)
                if(ios.ne.0) then
                  ! fallback: formatted
                  backspace(ifile)
                  read(ifile,101) (coeff(j,l),l=1,5)
                  if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest
                  write(*,*) 'error in character to coeffb trafo'
                  write(*,*) dummy
                  write(*,*) args(nargs-5+k)
                  stop 'error in character to coeffb trafo'
                endif
              enddo
           enddo            

! now print out a set of 5 MOs

           do j=1,5
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Beta'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,*)coeff(k,j)
            enddo

           enddo
 
        enddo  

! same procedure for the last, incomplete block
       write(frmt,'(a,i1,a)') '(21x,',jrem,'f10.5)'
       if(jrem.gt.0) then
         ! read "symmetry" (here we just use orbital numbers)
           read(ifile,*) (ityp(j),j=1,jrem)
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (MO read)'
           call parse(dummy,' ',args,nargs)
           do j=1,jrem
              iocc(j)=0
              call parse(args(j),'-',atmp,k)
              if(trim(atmp(k)).eq.'O'.or.trim(atmp(k)).eq.'o')iocc(j)=1 ! get occupation
           enddo

! read orbital energies
           read(ifile,'(A)',IOSTAT=ios) dummy
           if(ios.lt.0) stop 'end of file (eigval read)'
           call parse(dummy,' ',args,nargs)
           if(trim(args(1)).ne.'Eigenvalues')stop 'error in eigval read'
           do j=1,jrem
              eval(j)=0.0d0
              call value(args(nargs-jrem+j),eval(j),ios)
              if(ios.ne.0) then
                ! fallback: formatted 
                backspace(ifile)
                read(ifile,frmt) (eval(k),k=1,jrem)
                if (abs(eval(j)-eval(j)).lt.1.d-5) exit ! do the same for the rest
                write(*,*) 'error in character to eigvalb trafo'
                write(*,*) dummy
                write(*,*) args(nargs-jrem+j)
                stop 'error in character to eigvalb trafo'
              endif
           enddo

! now read coefficients

           do j=1,nbf
              read(ifile,'(A)',IOSTAT=ios) dummy
              if(ios.lt.0) stop 'end of file (MO coeff read)'
              call parse(dummy,' ',args,nargs)
              do k=1,jrem
                coeff(j,k)=0.0d0
                call value(args(nargs-jrem+k),coeff(j,k),ios)
                if(ios.ne.0) then
                 ! fallback: formatted
                 backspace(ifile)
                 read(ifile,frmt) (coeff(j,l),l=1,jrem)
                 if (abs(coeff(j,k)-coeff(j,k)).lt.1.d-5) exit ! do the same for the rest                 
                 write(*,*) 'error in character to coeffb trafo'
                 write(*,*) dummy
                 write(*,*) args(nargs-jrem+k)
                 stop 'error in character to coeffb trafo'
                endif
              enddo
           enddo

! now print out the set of remaining MOs

           do j=1,jrem
            write(*,*)'Sym= ',ityp(j)
            write(*,*)'Ene= ',eval(j)
            write(*,'(A)')' Spin= Beta'
            write(*,*)'Occup= ',real(iocc(j))
            do k=1,nbf
             write(*,'(i6)',advance='no')k
             write(*,*)coeff(k,j)
            enddo

           enddo
      
       endif
       endif

       enddo

      endif

      deallocate( coeff ) 
     
      end

! search for string in input
      subroutine findstr(str,ifile,i,n)
      use strings
! str : string that should be searched for
! ifile : file unit where to read
! i : element in line where str is located
! n : total number of elements in line - if not known, set to zero
      implicit none
      character*(*), intent( in ) :: str
      integer, intent( in ) :: i,n,ifile
      character(len=:),allocatable :: arg(:)
      integer narg,ios
      real*8 tmp
      character*79 line
      if(n.lt.i.and.n.ne.0) stop 'error in FINDSTR: i > n!'
      allocate(character(len(str)) :: arg(79),stat=ios)
      if(ios.ne.0) stop 'allocation in findstr failed'
      if(n.le.0) then
       do
        read(ifile,'(A)',IOSTAT=ios) line
        if(ios.lt.0) write(*,*) 'No ',str,' flag found'
        if(ios.lt.0) stop 'end of file reached'
        call parse(line,' ',arg,narg)
        if(narg.lt.i) cycle
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
      deallocate(arg,stat=ios)
      if(ios.ne.0) stop 'deallocation in findstr failed'
      end subroutine findstr

      subroutine getelement(iatm,atm)
      implicit none
      integer, intent ( in ) :: iatm
      character*2, intent ( out ) :: atm
      select case(iatm)
       case( 1 )
        atm='H '
       case( 2 )
        atm ='He'
       case( 3 )
        atm ='Li'
       case( 4 )
        atm ='Be'
       case( 5 )
        atm ='B '
       case( 6 )
        atm ='C '
       case( 7 )
        atm ='N '
       case( 8 )
        atm ='O '
       case( 9 )
        atm ='F '
       case( 10 )
        atm ='Ne'
       case( 11 )
        atm ='Na'
       case( 12 )
        atm ='Mg'
       case( 13 )
        atm ='Al'
       case( 14 )
        atm ='Si'
       case( 15 )
        atm ='P '
       case( 16 )
        atm ='S '
       case( 17 )
        atm ='Cl'
       case( 18 )
        atm ='Ar'
       case( 19 )
        atm ='K '
       case( 20 )
        atm ='Ca'
       case( 21 )
        atm ='Sc'
       case( 22 )
        atm ='Ti'
       case( 23 )
        atm ='V '
       case( 24 )
        atm ='Cr'
       case( 25 )
        atm ='Mn'
       case( 26 )
        atm ='Fe'
       case( 27 )
        atm ='Co'
       case( 28 )
        atm ='Ni'
       case( 29 )
        atm ='Cu'
       case( 30 )
        atm ='Zn'
       case( 31 )
        atm ='Ga'
       case( 32 )
        atm ='Ge'
       case( 33 )
        atm ='As'
       case( 34 )
        atm ='Se'
       case( 35 )
        atm ='Br'
       case( 36 )
        atm ='Kr'
       case( 37 )
        atm ='Rb'
       case( 38 )
        atm ='Sr'
       case( 39 )
        atm ='Y '
       case( 40 )
        atm ='Zr'
       case( 41 )
        atm ='Nb'
       case( 42 )
        atm ='Mo'
       case( 43 )
        atm ='Tc'
       case( 44 )
        atm ='Ru'
       case( 45 )
        atm ='Rh'
       case( 46 )
        atm ='Pd'
       case( 47 )
        atm ='Ag'
       case( 48 )
        atm ='Cd'
       case( 49 )
        atm ='In'
       case( 50 )
        atm ='Sn'
       case( 51 )
        atm ='Sb'
       case( 52 )
        atm ='Te'
       case( 53 )
        atm ='I '
       case( 54 )
        atm ='Xe'
       case( 55 )
        atm ='Cs'
       case( 56 )
        atm ='Ba'
       case( 57 )
        atm ='La'
       case( 58 )
        atm ='Ce'
       case( 59 )
        atm ='Pr'
       case( 60 )
        atm ='Nd'
       case( 61 )
        atm ='Pm'
       case( 62 )
        atm ='Sm'
       case( 63 )
        atm ='Eu'
       case( 64 )
        atm ='Gd'
       case( 65 )
        atm ='Tb'
       case( 66 )
        atm ='Dy'
       case( 67 )
        atm ='Ho'
       case( 68 )
        atm ='Er'
       case( 69 )
        atm ='Tm'
       case( 70 )
        atm ='Yb'
       case( 71 )
        atm ='Lu'
       case( 72 )
        atm ='Hf'
       case( 73 )
        atm ='Ta'
       case( 74 )
        atm ='W '
       case( 75 )
        atm ='Re'
       case( 76 )
        atm ='Os'
       case( 77 )
        atm ='Ir'
       case( 78 )
        atm ='Pt'
       case( 79 )
        atm ='Au'
       case( 80 )
        atm ='Hg'
       case( 81 )
        atm ='Tl'
       case( 82 )
        atm ='Pb'
       case( 83 )
        atm ='Bi'
       case( 84 )
        atm ='Po'
       case( 85 )
        atm ='At'
       case( 86 )
        atm ='Rn'
       case( 87 )
        atm ='Fr'
       case( 88 )
        atm ='Ra'
       case( 89 )
        atm ='Ac'
       case( 90 )
        atm ='Th'
       case( 91 )
        atm ='Pa'
       case( 92 )
        atm ='U '
       case( 93 )
        atm ='Np'
       case( 94 )
        atm ='Pu'
       case( 95 )
        atm ='Am'
       case( 96 )
        atm ='Cm'
       case( 97 )
        atm ='Bk'
       case( 98 )
        atm ='Cf'
       case( 99 )
        atm ='Es'
       case( 100 )
        atm ='Fm'
       case( 101 )
        atm ='Md'
       case( 102 )
        atm ='No'
       case( 103 )
        atm ='Lr'
       case( 104 )
        atm ='Rf'
       case( 105 )
        atm ='Db'
       case( 106 )
        atm ='Sg'
       case( 107 )
        atm ='Bh'
       case( 108 )
        atm ='Hs'
       case( 109 )
        atm ='Mt'
       case( 110 )
        atm ='Ds'
       case( 111 )
        atm ='Rg'
       case( 112 )
        atm ='Cn'
!       case( 113 )
!        atm ='Uut'
       case( 114 )
        atm ='Fl'
!       case( 115 )
!        atm ='Uup'
       case( 116 )
        atm ='Lv'
!       case( 117 )
!        atm ='Uus'
!       case( 118 )
!        atm ='Uuo'
       case default
        write (*,*) 'non-recognizable element with number:',iatm
        atm='XX'
      end select
      return
      end

     
