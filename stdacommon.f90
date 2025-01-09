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
! wavefunction and basis common block
      module stdacommon

      real*8, allocatable :: co(:,:),exip(:),occ(:),eps(:)
      real*8, allocatable :: cxip(:),eta(:,:)
      integer, allocatable :: ipat(:),ipty(:),ipao(:),iaoat(:)
      character*2, allocatable :: atnam(:)
      logical :: spherical
      real*8 :: rsh,rsh_ax,rsh_beta,rsh2
      end module stdacommon

! xtb common block
      module kshiftcommon

      real*8 shftmax,shftwidth,shftsteep,shftmax_somo
      ! shftsteep: steepness, i.e. power of decay
      ! shftmax: maximum shift
      ! shftwidth: shiftwidth for damping

      end module kshiftcommon

! some logicals in order not to blow-up subroutine calls
      module commonlogicals
      logical triplet,rpachk,eigvec,screen,dokshift,printexciton,velcorr
      logical aniso
      logical resp,TPA,aresp,ESA,smp2,bfw,spinflip,rw,pt_off,nto,sf_s2,multipole
      logical optrota,velo_OR,rw_dual,XsTD,cint,RSH_flag,RSH_sub,SOS_2PA,Xcore
      logical FULL2PA,full,direct_full
      end module commonlogicals

! some variables for the response functions
      module commonresp
      integer :: num_freq
      integer :: num_trans
      integer :: state2opt
      integer :: Nnto, Ecore, Ecore2
      end module commonresp

      module commonlibcint
      integer :: n_env, nbas,n_at
      integer, allocatable :: atm(:,:), bas(:,:)
      integer, allocatable :: di_all(:), sum_di(:)
      double precision, allocatable :: env(:)
      integer,allocatable :: bas_start(:),bas_end(:),prim_at(:),counter_at(:)
      integer, allocatable :: Q_atom(:)
      end module commonlibcint
