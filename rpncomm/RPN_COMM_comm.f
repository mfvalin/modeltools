
*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */

      integer function RPN_COMM_comm(com)
c	Luc Corbeil, 2000-11-21
c
c	lien entre chaine de caractere de communicateur
c	GRID, EW et NS et leur numero assigne par
c	MPI.
c
      use rpn_comm
      implicit none
!      include mpif.h
!        include rpn_comm.h
      character(len=*) com
      character(len=32) comm

      call rpn_comm_low2up(com,comm)

      RPN_COMM_comm = MPI_COMM_NULL

      if (trim(comm) == RPN_COMM_GRIDPEERS) then
         RPN_COMM_comm=pe_grid_peers
         return
      endif
      if (trim(comm) == RPN_COMM_GRID) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if(trim(comm) == RPN_COMM_DOMM) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if (trim(comm) == RPN_COMM_WORLD) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if(trim(comm) == RPN_COMM_ALLDOMAINS) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if (trim(comm) == RPN_COMM_ALLGRIDS) then
         RPN_COMM_comm=pe_a_domain
         return
      endif
      if (trim(comm) == RPN_COMM_MULTIGRID) then
         RPN_COMM_comm=pe_indomms ! alias pe_multi_grid
         return
      endif
      if (trim(comm) == RPN_COMM_ALL) then
        RPN_COMM_comm=pe_wcomm  ! alias pe_grid
        return
      endif
      if(trim(comm) == RPN_COMM_DEFAULT) then
        RPN_COMM_comm=pe_defcomm
        return
      endif
      if (trim(comm) == RPN_COMM_EW) then
         RPN_COMM_comm=pe_myrow
         return
      endif
      if (trim(comm) == RPN_COMM_NS) then
         RPN_COMM_comm=pe_mycol
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCMASTER) then
         RPN_COMM_comm=pe_blocmaster
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCK) then
         RPN_COMM_comm=pe_bloc
         return
      endif
      if (trim(comm) == RPN_COMM_UNIVERSE) then
         RPN_COMM_comm=MPI_COMM_WORLD
         return
      endif

      write(rpn_u,*) 'Unknown communicator ',comm,', aborting'
      stop
      return
      end
