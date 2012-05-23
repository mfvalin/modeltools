!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
module rpn_comm
  use iso_c_binding
  save
  integer MAX_OPTN
  parameter(MAX_OPTN=10)
!
!	domain boundary flags, LOGICAL
!	.true. if a PE(TILE) is on a domain edge
!	all 4 values .true. if single tile
!	normally set by routine rpn_comm_init
!
	logical bnd_east,bnd_west,bnd_north,bnd_south
!
!	characteristics of the local PE
!	plus some PE grid topology information
!	normally set by routine rpn_comm_init
!
!	pe_me	PE number of this PE
!	pe_mex	 x coordinate of this PE in domain (origin 0)
!	pe_mey	 y coordinate of this PE in domain (origin 0)
!	pe_myrow communicator for PEs in same ROW(same pe_mey)
!	pe_mycol communicator for PEs in same COLUMN(same pe_mex)
!	pe_tot	 total number of PEs involved
!	pe_nx	 number of PEs along x axis
!	pe_ny	 number of PEs along y axis
!	pe_pe0	 PE number of first PE in grid
!	pe_extra flag =1 if pe in compute grid, =0 otherwise
!	pe_grid	 communicator for PEs in grid and out of grid
!       pe_ingrid  communicator for PEs in grid
!       pe_outgrid communicator for PEs out of grid
!	pe_bloc	 communicator for internal blocs in grid
!       pe_blocmaster communicator for bloc_corner PEs in grid
!	pe_id	 matrix of PE numbers in grid
!		 pe_me=pe_id(pe_mex,pe_mey)
!	pe_optn  table of options ( option name)
!	pe_opiv  values of integer options
!	pe_oprv  values of real options
!	pe_opcv  values of character options
!
!       WORLD_COMM_MPI  replaces MPI_COMM_WORLD
!
  integer WORLD_COMM_MPI
  integer diag_mode
  logical :: WORLD_COMM_MPI_INIT=.false.
  logical :: RPM_COMM_IS_INITIALIZED=.false.
  integer :: deltai=1   ! deltai and deltaj are used by RPN_COMM_petopo to distribute grid over PEs
  integer :: deltaj=1   ! default PE distribution is X increasing, then Y increasing
  integer pe_me,pe_mex,pe_mey,pe_myrow,pe_mycol
  integer pe_tot,pe_nx,pe_ny,pe_pe0,pe_extra
  integer pe_gr_extra,pe_gr_myrow,pe_gr_mycol
  integer pe_bloc, pe_blocmaster, pe_defgroup
  integer pe_gr_bloc, pe_gr_blocmaster, pe_defcomm
  integer pe_gr_indomm, pe_gr_outdomm, pe_indomm, pe_outdomm
  integer pe_gr_indomms, pe_indomms ! multigrid countepart to pe_gr_indomm and pe_indomm
  integer pe_wcomm, pe_gr_wcomm
  integer pe_wcomms, pe_gr_wcomms  ! multigrid countepart to pe_wcomm and pe_gr_wcomm
  integer pe_dommtot, pe_medomm
  integer pe_all_domains, pe_gr_all_domains      ! all the domains
  integer pe_me_all_domains, pe_tot_all_domains  ! all the domains
  integer pe_a_domain, pe_gr_a_domain            ! all multigrids in a domain
  integer pe_me_a_domain, pe_tot_a_domain        ! all multigrids in a domain
  integer pe_multi_grid, pe_gr_multi_grid        ! all the grids in a multigrid
  integer pe_me_multi_grid, pe_tot_multi_grid    ! all the grids in a multigrid
  integer pe_grid, pe_gr_grid                    ! a single grid
  integer pe_me_grid, pe_tot_grid                ! a single grid
  integer pe_grid_peers, pe_gr_grid_peers        ! PEs with same rank on different grids of same multigrid
  integer pe_me_peer, pe_tot_peer                ! PEs with same rank on different grids of same multigrid
  integer pe_grid_host, pe_gr_grid_host          ! PEs on same host node (belonging to same "grid")
  integer pe_me_grid_host, pe_tot_grid_host      ! PEs on same host node (belonging to same "grid")
  integer my_colors(3)                           ! domain/multigrid/grid color
  integer, allocatable, dimension(:,:) :: pe_id
  integer, allocatable, dimension(:)   :: pe_xtab,pe_ytab
  integer, allocatable, dimension(:,:) :: pe_location   ! pe_x,pe_x,pe_ingrid,pe_insgrid,pe_indomain
  logical :: async_exch=.true.
  character *4 pe_optn(MAX_OPTN)
  integer pe_opiv(MAX_OPTN)
  real *4 pe_oprv(MAX_OPTN)
  character *4 pe_opcv(MAX_OPTN)
  integer, private :: RESTE
  parameter (RESTE=MAX_OPTN-1)
  data pe_optn/'FILL',RESTE*'    '/
  data pe_opcv/MAX_OPTN*'    '/
  data pe_oprv/MAX_OPTN*0.0/
  data pe_opiv/MAX_OPTN*0/

  integer ord_max,ord_me
  integer, allocatable, dimension(:) :: ord_tab
!
!	GLOBAL information, will be BROADCAST
!
!	WORLD_pe(1) number of PEs along x in grid
!	WORLD_pe(2) number of PEs along y in grid
!	WORLD_pe(3) deltai, size of PE blocks along x
!	WORLD_pe(4) deltaj, size of PE blocks along Y
!	WORLD_pe(5:10) provision for future expansion
!
  integer WORLD_pe(10)
!
!	TIMING information (being reworked at this time)
!
  integer TMG_CLOCK, TMG_CPU, TMG_FLOP, TMG_VEC
  parameter(TMG_CLOCK=1,TMG_CPU=2,TMG_VEC=3,TMG_FLOP=4)
  integer TMG_CLOCK_C, TMG_CPU_C
  parameter(TMG_CLOCK_C=5,TMG_CPU_C=6)
  !
  integer MAXTMG,MAXTMGELEM
  parameter (MAXTMG=1024)
  parameter (MAXTMGELEM=8)
  real *8 tmg_tbl(MAXTMGELEM,MAXTMG)
  real *8 tmg_old(MAXTMGELEM),tmg_new(MAXTMGELEM)
  integer tmg_indx
  integer, private :: MAXWDS
  parameter(MAXWDS=MAXTMG*MAXTMGELEM)
  data tmg_indx / 0 /
  data tmg_old /MAXTMGELEM * 0.0/
  data tmg_new /MAXTMGELEM * 0.0/
  data tmg_tbl /MAXWDS * 0.0/
!
!       Subgrid information
!
  logical BLOC_EXIST
  integer BLOC_SIZEX,BLOC_SIZEY, BLOC_MASTER
  integer BLOC_mybloc,BLOC_myblocx,BLOC_myblocy
  integer BLOC_me,BLOC_corner
  integer BLOC_comm_world, bloc_comm_row, bloc_comm_col
!
! Output unit information
!
  integer :: rpn_u=6
  logical :: rpn_ew_ext_l=.false.
!
! Domain type
!
  type domm
     sequence
     character(len=12) nom
     character(len=1024) path
     integer npex, npey
  end type domm
  integer domm_size, domm_num
  type(domm), allocatable, dimension(:) :: pe_domains

!
! interface to gethostid
! get the 32 bit host identifier
!
  interface
    integer(C_INT) function f_gethostid()BIND(C,name='gethostid')
      use iso_c_binding
    end function f_gethostid
  end interface


end module rpn_comm
