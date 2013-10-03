*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
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
      subroutine RPN_COMM_test_limit()
      implicit none
      integer, parameter :: NPE=6
      integer, dimension(NPE) :: count, offset
      integer :: my_id, gmin, gmax, lmini, lmaxi, status
      integer :: RPN_COMM_limit_2
      external :: RPN_COMM_limit_2

      gmin = 1
      gmax = 13
      do my_id=0,NPE-1
      status = RPN_COMM_limit_2(my_id, npe, gmin, gmax,
     &     lmini,lmaxi,count, offset,3)
      print 101, 'pe_me=',my_id, lmini,lmaxi, gmin, gmax, status
      print 101, 'count=',count
      print 101, 'offst=',offset
      print *,''
      enddo
      print *, '---------------------------'
101   format(A7,10I5)
      my_id=NPE-1
      status = RPN_COMM_limit_2(my_id, npe, gmin, gmax,
     &     lmini,lmaxi,count, offset,1)
      print 101, 'pe_me=',my_id, lmini,lmaxi, status
      print 101, 'count=',count
      print 101, 'offst=',offset
      print *,''
      gmax = 6
      status = RPN_COMM_limit_2(my_id, npe, gmin, gmax,
     &     lmini,lmaxi,count, offset,1)
      print 101, 'pe_me=',my_id, lmini,lmaxi, status
      print 101, 'count=',count
      print 101, 'offst=',offset
      print *,''
      status = RPN_COMM_limit_2(my_id, npe, gmin, gmax,
     &     lmini,lmaxi,count, offset,3)
      print 101, 'pe_me=',my_id, lmini,lmaxi, status
      print 101, 'count=',count
      print 101, 'offst=',offset
      print *,''
      end subroutine RPN_COMM_test_limit
***function  RPN_COMM_limit_2 global domain decomposition function (along one dimension)
      integer function RPN_COMM_limit_2(my_id, npe, gmin, gmax, 
     &     lmini,lmaxi,count, offset,relax)
      implicit none
*
*arguments
*  I	my_id   "tile" ordinal along decomposition axis (0 to npe-1)
*  I    npe    number of "tiles" (PEs) along this dimension
*  I    gmin,gmax
*              global index space along this dimension is gmin:gmax (usually 1:n)
*  O    lmini,lmaxi
*              this "tile" will cover index range lmini:lmaxi in global space
*  O    count(1:npe)
*              count(i) = number of points along this dimension for PE with ordinal I-1
*  O   offset(1:npe)
*              offset(i) = offset from gmin for PE with ordinal I-1
*  I    relax  decomposition mode
*          0 : strict mode, all tiles but last one must have same dimension, 
*              last tile may be shorter but may not have zero dimansion
*          1 : some tiles at end 1 shorter than tiles at beginning, zero dimension not allowed
*          2 : same as relax=1 but zero dimension tiles at end are allowed
*          3 : tiles with same length followed by a shorter tile followed by zero size tiles
*
*notes
*     this function is totally stand alone and could eventually be moved into the rmnlib library
**
      integer, intent(IN) ::  my_id, npe, gmin, gmax, relax
      integer, intent(OUT) :: lmini,lmaxi
      integer, intent(OUT) :: count(npe),offset(npe)

      integer gtot
      integer val1, val2, i

      lmini = -1
      lmaxi = -1
      gtot = gmax - gmin + 1    ! number of points to be distributed

      val1 = (gtot + npe - 1)/npe           ! ceiling(gtot/npe)
      count = val1
      val2 = val1
      offset(1) = 0
      do i=2,npe
        count(i) = min(count(i),gtot-val2)
        offset(i) = offset(i-1) + count(i-1)
        val2 = val2 + count(i)
      enddo
      if(gtot < npe .and. relax < 2) then   ! there would be zero sized tiles and relax is not 2 or 3
          goto 777
      end if
      if(relax == 3) goto 666               ! 
      val2 = gtot - (npe-1)*val1            ! potential size of last "tile" in strict mode
      if (val2 < 0  .and. relax == 0)  then ! in STRICT mode last "tile" would have a negative size
          goto 777
      end if
      if(val2 <= 0 ) val1 = gtot/npe        ! relax distribution rule, try mode 1
      count = val1
      if(val2 > 0) then                     ! strict mode will work, use it
         count(npe) = val2
      else                                  ! relaxed mode (1 or 2)
        do i=1,mod(gtot,npe)                ! add 1 to the size of some tiles at beginning
          count(i)=count(i)+1
        end do
      end if
666   do i= 2, npe   ! use count table to compute offsets from beginning
         offset(i) = offset(i-1) + count(i-1)
      end do 
      lmini = gmin + offset(my_id+1)             ! start of "tile" in global space
      lmaxi = min(gmax,lmini+count(my_id+1)-1)   ! end of "tile" in global space

      RPN_COMM_limit_2 = 0  ! SUCCESS
      return

777   continue      ! something bad happened, print error message and return error status
      print *, 'RPN_COMM_limit_2: invalid decomposition'
      print *, 'Nb of elements =', gtot
      print *, 'Nb of pe =', npe
      print *, 'relax =', relax
      RPN_COMM_limit_2 = -1
      return
      end
*     old function, calls newer RPN_COMM_limit_2 forcing STRICT decomposition mode
*     kept for compatibility with older versions of this library
      integer function RPN_COMM_limit(my_id, npe, gmin, gmax, lmini,   ! old entry point
     &     lmaxi,count, offset)
      implicit none
      integer, intent(IN) ::  my_id, npe, gmin, gmax
      integer, intent(OUT) :: lmini,lmaxi
      integer, intent(OUT) :: count(npe),offset(npe)
      external RPN_COMM_limit_2
      integer RPN_COMM_limit_2
      RPN_COMM_limit = 
     &     RPN_COMM_limit_2(my_id, npe, gmin, gmax, lmini,
     &     lmaxi,count, offset,0) ! call with relax=0 (strict mode)
      return
      end
