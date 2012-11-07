!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
!=======================================================================
      integer function RPN_COMM_xch_halo_test()
!=======================================================================
      use rpn_comm
      implicit none
!
      integer, pointer, dimension(:,:,:) :: localarray
      integer*8, pointer, dimension(:,:,:) :: localarray2
      integer, parameter :: lni=3
      integer, parameter :: lnj=5
      integer, parameter :: nk=7
      integer :: gni, gnj
      integer :: i, j, k, ierr
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: minx, maxx, miny, maxy
      integer, dimension(pe_nx) :: countx, offsetx
      integer, dimension(pe_ny) :: county, offsety
      integer :: halox, haloy, npol_row
      logical :: periodx, periody

!
      gni = pe_nx*lni
      gnj = pe_ny*lnj
      halox=1
      haloy=2
      call RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      call RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      minx = lminx-halox
      maxx = lmaxx+halox
      miny = lminy-haloy
      maxy = lmaxy+haloy
      if(pe_me==8) write(rpn_u,100)'grid halo exchange test',
     %    pe_tot_grid,pe_nx,pe_ny,lminx,lmaxx,lminy,lmaxy,countx,county,
     %    minx,maxx,miny,maxy
100   format(A,25I5)
      allocate(localarray(minx:maxx,miny:maxy,nk))
      allocate(localarray2(minx:maxx,miny:maxy,nk))
!
      localarray = -1
      localarray2 = -1
      do k = 1,nk
      do j = lminy,lmaxy
      do i = lminx,lmaxx
        localarray(i,j,k) = k + 10*j + 1000*i
        localarray2(i,j,k) = k + 10*j + 1000*i
      enddo
      enddo
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
      periodx = .false.
      periody = .false.
      npol_row = 0
!      return
      call RPN_COMM_xch_halo(localarray,minx,maxx,miny,maxy,
     %             lni,lnj,nk,halox,haloy,periodx,periody,
     %             gni,npol_row)
!
      return
      end function RPN_COMM_xch_halo_test
!=======================================================================
      SUBROUTINE RPN_COMM_xch_halo(g,minx,maxx,miny,maxy,
     %             ni,nj,nk,halox,haloy,periodx,periody,
     %             gni,npol_row)
!=======================================================================
!
!  +=======================================================================================+___maxy
!  I       <halox>                                                           <halox>       I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I :     :                                    |                                  :     : I
!  I :"cTL":               "iTLCR"            haloy                                :"cTR": I
!  I :     :                                    |                                  :     : I
!  I +-----+=====+===========================================================+=====+-----: I___nj
!  I :     I     :                              |                            :     I     : I
!  I :"iTL"I" TL":            "TC"            haloy                          :"TR "I"iTR": I
!  I :     I     :                              |                            :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :"iCL"I" CL":                                                           :"CR "I"iCR": I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                              |                            :     I     : I
!  I :"iBL"I" BL":            "BC"            haloy                          :"BR "I"iBR": I
!  I :     I     :                              |                            :     I     : I
!  I +-----+=====+===========================================================+=====+-----: I___1
!  I :     :                                    |                                  :     : I
!  I :"cBL":               "iBLCR"            haloy                                :"cBR": I
!  I :     :                                    |                                  :     : I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I <halox>                                                                       <halox> I
!  +=======================================================================================+___miny
!  |       |                                                                       |       |
!  |       |                                                                       |       |
!  |       1                                                                       ni      |
!  minx                                                                                 maxx
!
!    Description of the 8 internal arrays 
!
!  BR_CR_TR array     : 3 pieces  BR[nk*halox*haloy],  CR[nk*halox*(nj-2*haloy)],  TR[nk*halox*haloy]
!                       used for West to East send  (nk*halox*nj)
!  iBR_iCR_iTR array  : 3 pieces iBR[nk*halox*haloy], iCR[nk*halox*(nj-2*haloy)], iTR[nk*halox*haloy]
!                       used for East to West receive  (nk*halox*nj)
!                       (iBR, iTR then used for North/South send)
!  BL_CL_TL array     : 3 pieces  BL[nk*halox*haloy],  CL[nk*halox*(nj-2*haloy)],  TL[nk*halox*haloy]
!                       used for East to West send  (nk*halox*nj)
!  iBL_iCL_iTL array  : 3 pieces iBL[nk*halox*haloy], iCL[nk*halox*(nj-2*haloy)], iTL[nk*halox*haloy]
!                       used for West to East receive  (nk*halox*nj)
!                       (iBL, iTL then used for North/South send)
!  iBL_BLCR_iBR array : 3 pieces iBL[nk*halox*haloy],  BLCR[nk*haloy*ni],         iBR[nk*halox*haloy]
!                       used for North to South send (nk*haloy*(ni+2*halox))
!  iTL_TLCR_iTR array : 3 pieces iTL[nk*halox*haloy],  TLCR[nk*haloy*ni],         iTR[nk*halox*haloy]
!                       used for South to North send (nk*haloy*(ni+2*halox))
!  cBL_iBLCR_cBR array: 3 pieces cBL[nk*halox*haloy],  iBLCR[nk*haloy*ni],        cBR[nk*halox*haloy]
!                       used for South to North receive (nk*haloy*(ni+2*halox))
!  cTL_iTLCR_cTR array: 3 pieces cTL[nk*halox*haloy],  iTLCR[nk*haloy*ni],        cTR[nk*halox*haloy]
!                       used for North to South receive (nk*haloy*(ni+2*halox))
!
!    Communication tags (for iSEND/Irecv)
!
!  West  -> East  : sender_pe + 100000
!  West  <- East  : sender_pe
!  North -> South : sender_pe + 100000
!  North <- South : sender_pe
!
      use rpn_comm
      implicit none
!
!     exchange a halo with N/S/E/W neighbours
!
      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer gni,npol_row
      logical periodx,periody
!     integer *8 mem_time, exch_time, ewtime
      integer g(minx:maxx,miny:maxy,nk)
!
        integer, dimension(halox*nj*nk) :: BR_CR_TR, BL_CL_TL  ! send to East/West
        integer, dimension(halox*nj*nk) :: iBR_iCR_iTR, iBL_iCL_iTL   ! recv from East/West
        integer, dimension(haloy*(ni+2*halox)*nk) ::
     %           iBL_BLCR_iBR, iTL_TLCR_iTR,    ! send to South/North
     %           cBL_iBLCR_cBR, cTL_iTLCR_cTR   ! recv from South/North
        integer :: bl,bc,br,cl,cr,tl,tc,tr
        integer :: cbl,iblcr,cbr,ibl,ibr,icl,icr
        integer :: itl,itr,ctl,itlcr,ctr
        integer :: to_east, to_west, to_north, to_south
        integer :: from_east, from_west, from_north, from_south
!     integer *8 time_base,temp_time
      integer i, j, k, m
      integer nwds_ew, nwds_ns
      integer sendtag, gettag, ierr
      integer status(MPI_STATUS_SIZE)
      logical east,west,north,south
      integer eastpe,westpe,northpe,southpe
!
      integer globalni,polarrows, nilmax, jmin,jmax
        integer RPN_COMM_topo, mini,maxi,nil,ni0
!
      integer land_fill
      real r_land_fill
      equivalence(land_fill,r_land_fill)
!
      globalni=abs(gni)
      polarrows=npol_row
      nwds_ew = size(BR_CR_TR)
      nwds_ns = size(iBL_BLCR_iBR)

1     continue
      
!     call RPN_COMM_tmg_in
      east=(bnd_east) .and. (.not.periodx)
      eastpe=pe_id(pe_mex+1,pe_mey)
      west=(bnd_west) .and. (.not.periodx)
      westpe=pe_id(pe_mex-1,pe_mey)
      north=(bnd_north) .and. (.not.periody)
      northpe=pe_id(pe_mex,pe_mey+1)
      south=(bnd_south) .and. (.not.periody)
      southpe=pe_id(pe_mex,pe_mey-1)

      jmin = 1
      jmax = nj 
      if(rpn_ew_ext_L) then
         if(north) jmax = nj+haloy
         if(south) jmin = 1-haloy
      endif
!
      if(pe_opcv(1) .ne. ' ') then !  fill halo option present
           r_land_fill=pe_oprv(1)
           if(pe_opcv(1) .eq. 'BAND') then
              if(iand(1,pe_opiv(1)) .ne. 0) then ! south band
                 do j=miny,0
                    do k=1,nk
                       do i=minx,maxx
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
              endif
            if(iand(2,pe_opiv(1)) .ne. 0) then ! east band
                 do i=ni+1,maxx
                    do j=miny,maxy
                       do k=1,nk
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
            if(iand(4,pe_opiv(1)) .ne. 0) then ! north band
                 do j=nj+1,maxy
                    do k=1,nk
                       do i=minx,maxx
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
            if(iand(8,pe_opiv(1)) .ne. 0) then ! west band
                 do i=minx,0
                    do j=miny,maxy
                       do k=1,nk
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
           endif
           if(pe_opcv(1) .eq. 'EDGE') then
              if(iand(1,pe_opiv(1)) .ne. 0) then ! south edge
                 do i=minx,0
                    do k=1,nk
                       g(i,1,k)=land_fill
                    enddo
                 enddo
                 do i=ni+1,maxx
                    do k=1,nk
                       g(i,1,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(2,pe_opiv(1)) .ne. 0) then ! east edge
                 do j=miny,0
                    do k=1,nk
                       g(ni,j,k)=land_fill
                    enddo
                 enddo
                 do j=nj+1,maxy
                    do k=1,nk
                       g(ni,j,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(4,pe_opiv(1)) .ne. 0) then ! north edge
                 do i=minx,0
                    do k=1,nk
                       g(i,nj,k)=land_fill
                    enddo
                 enddo
                 do i=ni+1,maxx
                    do k=1,nk
                       g(i,nj,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(8,pe_opiv(1)) .ne. 0) then ! west edge
                 do j=miny,0
                    do k=1,nk
                       g(1,j,k)=land_fill
                    enddo
                 enddo
                 do j=nj+1,maxy
                    do k=1,nk
                       g(1,j,k)=land_fill
                    enddo
                 enddo
            endif
           endif
        endif
!
!       use new fullly asynchronous code only if pe_nx and pe_ny both >1
!       and polarrows = 0
!
      if(pe_nx>1 .and. pe_ny>1 .and. polarrows<=0) goto 2
!
!       if no halo along x, bypass
!     call tmg_start(90,'RPN_COMM_haloew')
      if (halox .gt. 0) then
           if (.not.(min(pe_mey+1,pe_ny-pe_mey).le.polarrows)) then
              call RPN_COMM_xch_haloew(g,minx,maxx,miny,maxy,
     %             ni,jmin,jmax,nk,halox,haloy,periodx,periody)
              endif
      endif
!     call tmg_stop(90)
!     call tmg_start(91,'RPN_COMM_halons')
      if (haloy .gt. 0) then
              call RPN_COMM_xch_halons(g,minx,maxx,miny,maxy,
     %             ni,nj,nk,halox,haloy,periodx,periody)
      endif
!     call tmg_stop(91)
      if (min(pe_mey+1,pe_ny-pe_mey).le.polarrows) then
           ierr = RPN_COMM_topo(gni,mini,maxi,nil,nilmax,
     %          halox,ni0,.TRUE.,.FALSE.)
!          call tmg_start(92,'RPN_COMM_xch_halosl')
           call RPN_COMM_xch_halosl(g,minx,maxx,miny,maxy,
     %             ni,nj,nk,halox,haloy,periodx,periody,
     %             gni,npol_row,nilmax)
!          call tmg_stop(92)

      endif
!     call RPN_COMM_tmg_out
      return
!
!       version no 2 of this routine, fully asynchronous, overlapped communication /data movement
!
!       step 1 post non blocking receives in the East/West direction
!
2       continue
        if(.not. west) then
          call MPI_IRECV(iBL_iCL_iTL,nwds_ew,MPI_INTEGER,westpe, ! get from west neighbor unless i am west PE
     %         100000+westpe,PE_DEFCOMM,from_west,ierr)          ! sender was westpe, tag is westpe+100000
        endif
        if(.not. east) then
          call MPI_IRECV(iBR_iCR_iTR,nwds_ew,MPI_INTEGER,eastpe, ! get from east neighbor unless i am east PE
     %         eastpe,PE_DEFCOMM,from_east,ierr)                 ! sender was eastpe, tag is eastpe
        endif
!
!       step 2 fill East/West send buffers from inner halo
!
        bl = 1                       ; br = bl ; ibl = bl ; ibr = br
        cl = 1 + nk*halox*haloy      ; cr = cl ; icl = cl ; icr = cr
        tl = 1 + nk*halox*(ni-haloy) ; tr = tl ; itl = tl ; itr = tr
        m = bl
        do k=1,nk   ! BL/BL part
        do j=1,haloy
        do i=1,halox
           BL_CL_TL(m)=g(i,j,k)
           BR_CR_TR(m)=g(ni-halox+i,j,k)
           m=m+1
        enddo
        enddo
        enddo
        do k=1,nk  ! CL/CR part
        do j=1+haloy,nj-haloy
        do i=1,halox
           BL_CL_TL(m)=g(i,j,k)
           BR_CR_TR(m)=g(ni-halox+i,j,k)
           m=m+1
        enddo
        enddo
        enddo
        do k=1,nk  ! TL/TR part
        do j=nj-haloy+1,nj
        do i=1,halox
           BL_CL_TL(m)=g(i,j,k)
           BR_CR_TR(m)=g(ni-halox+i,j,k)
           m=m+1
        enddo
        enddo
        enddo
        if(m /= 1+halox*nj*nk) print *,'OUCH'
!      write(rpn_u,100),'DBG:',ni,nj,nk,halox,haloy,size(BL_CL_TL),
!     %      size(BR_CR_TR),m-1,halox*nj*nk
100   format(A,30I5)
!
!       step 3 non blocking send to East/West partners
!
        if(.not. east) then
          call MPI_ISEND(BR_CR_TR,nwds_ew,MPI_INTEGER,eastpe,   ! send to east neighbor unless i am east PE
     %         100000+pe_medomm,PE_DEFCOMM,to_east,ierr)         ! tag is PE grid ordinal of sender+100000
        endif
        if(.not. west) then
          call MPI_ISEND(BL_CL_TL,nwds_ew,MPI_INTEGER,westpe,  ! send to west neighbor unless i am west PE
     %         pe_medomm,PE_DEFCOMM,to_west,ierr)      ! tag is PE grid ordinal of sender
        endif
!
!       step 4 start filling the North/South send buffers
!
        cbl = 1                       ; ctl = 1
        iblcr = 1 + nk*halox*haloy    ; itlcr = iblcr
        cbr = 1 + nk*haloy*(ni+halox) ; ctr = cbr
        m = iblcr
        do k=1,nk  ! TL/TR part
        do j=1,haloy
        do i=1,ni
           iBL_BLCR_iBR(m)=g(i,j,k)
           iTL_TLCR_iTR(m)=g(ni-halox+i,j,k)
           m=m+1
        enddo
        enddo
        enddo
!
!       step 5 get East/West inbound data, then finish filling the North/South send buffers
!
        if(.not. east) then
          call MPI_wait(from_east,status,ierr)  ! wait for inbound East -> West message to complete
          do i = 0,halox*haloy*nk-1
             iBL_BLCR_iBR(ibr+i) = iBR_iCR_iTR(br+i)
             iTL_TLCR_iTR(itr+i) = iBR_iCR_iTR(tr+i)
          enddo
        endif
        if(.not. west) then
          call MPI_wait(from_west,status,ierr)  ! wait for inbound EAST <- West message to complete
          do i = 0,halox*haloy*nk-1             ! put iBL/iTL into N/S buffers
             iBL_BLCR_iBR(ibl+i) = iBL_iCL_iTL(bl+i)
             iTL_TLCR_iTR(itl+i) = iBL_iCL_iTL(tl+i)
          enddo
        endif
      write(rpn_u,100)'DBG',nwds_ns,
     %                size(cTL_iTLCR_cTR),size(cBL_iBLCR_cBR),
     %                size(iBL_BLCR_iBR),size(iTL_TLCR_iTR)
!
!       step 6 post North/South non blocking receives, post North/South non blocking sends
!
        if(.not. north) then
          call MPI_IRECV(cTL_iTLCR_cTR,nwds_ns,MPI_INTEGER,northpe,! recv from north neighbor unless I am north PE
     %         100000+northpe,PE_DEFCOMM,from_north,ierr)          ! sender was northpe therefore tag is northpe+100000
        endif
        if(.not. south) then
          call MPI_IRECV(cBL_iBLCR_cBR,nwds_ns,MPI_INTEGER,southpe,! recv from south neighbor unless I am south PE
     %         southpe,PE_DEFCOMM,from_south,ierr)                 ! sender was southpe therefore tag is southpe
        endif
        if(.not. south) then
          call MPI_ISEND(iBL_BLCR_iBR,nwds_ns,MPI_INTEGER,southpe, ! send to south neighbor unless I am south PE
     %         100000+pe_medomm,PE_DEFCOMM,to_south,ierr)          ! tag is PE grid ordinal of sender+100000
        endif
        if(.not. north) then
          call MPI_ISEND(iTL_TLCR_iTR,nwds_ns,MPI_INTEGER,northpe, ! send to north neighbor unless I am north PE
     %         pe_medomm,PE_DEFCOMM,to_north,ierr)                 ! tag is PE grid ordinal of sender
        endif
!
!       step 7  put East/West outer halos into array
!
        m = bl
        do k=1,nk   ! BL/BL part
        do j=1,haloy
        do i=1,halox
           g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iBR part
           g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iBL part
           m=m+1
        enddo
        enddo
        enddo
      if(.not. north)call MPI_wait(from_north,status,ierr)
      if(.not. south)call MPI_wait(from_south,status,ierr)
      return
        do k=1,nk  ! CL/CR part
        do j=1+haloy,nj-haloy
        do i=1,halox
           g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iCR part
           g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iCL part
           m=m+1
        enddo
        enddo
        enddo
        do k=1,nk  ! TL/TR part
        do j=nj-haloy+1,nj
        do i=1,halox
           g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iTR part
           g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iTL part
           m=m+1
        enddo
        enddo
        enddo
        if(m /= 1+halox*nj*nk) print *,'OUCH'
!
!       step 8 wait for inbound North -> South messages to complete 
!              and put into North outer halo
!
        if(.not. north) call MPI_wait(from_north,status,ierr)
        m=ctl
        do k=1,nk  ! cTL part
        do j=nj-haloy+1,nj
        do i=1,halox
           g(halox-i,j,k) = cTL_iTLCR_cTR(m)
           m=m+1
        enddo
        enddo
        enddo      
        do k=1,nk  ! iTLCR part
        do j=nj-haloy+1,nj
        do i=1,ni
           g(i,j,k) = cTL_iTLCR_cTR(m)
           m=m+1
        enddo
        enddo
        enddo      
        do k=1,nk  ! cTR part
        do j=nj-haloy+1,nj
        do i=1,halox
           g(ni+i,j,k) = cTL_iTLCR_cTR(m)
           m=m+1
        enddo
        enddo
        enddo      
        if(m /= 1+haloy*(ni+2*halox)*nk) print *,'OUCH'
!
!       step 9 wait for inbound North <- South messages to complete 
!              and put into South outer halo
!

        if(.not. south) call MPI_wait(from_south,status,ierr)
        m=cbl
        do k=1,nk  ! cBL part
        do j=nj-haloy+1,nj
        do i=1,halox
           g(halox-i,j,k) = cBL_iBLCR_cBR(m)
           m=m+1
        enddo
        enddo
        enddo      
        do k=1,nk  ! iBLCR part
        do j=nj-haloy+1,nj
        do i=1,ni
           g(i,j,k) = cBL_iBLCR_cBR(m)
           m=m+1
        enddo
        enddo
        enddo      
        do k=1,nk  ! cBR part
        do j=nj-haloy+1,nj
        do i=1,halox
           g(ni+i,j,k) = cBL_iBLCR_cBR(m)
           m=m+1
        enddo
        enddo
        enddo      
        if(m /= 1+haloy*(ni+2*halox)*nk) print *,'OUCH'
!
!       step 10 wait for all outbound messages to complete
!
        call MPI_wait(to_east,status,ierr)
        call MPI_wait(to_west,status,ierr)
        call MPI_wait(to_south,status,ierr)
        call MPI_wait(to_north,status,ierr)

        entry xch_halo(g,minx,maxx,miny,maxy,
     %             ni,nj,nk,halox,haloy,periodx,periody)
        globalni=ni
        polarrows=0

        goto 1
        end
