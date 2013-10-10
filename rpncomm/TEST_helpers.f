        subroutine TestUserInit(NX,NY) ! try to get NX,NY from file TEST.cfg if it exists
        external :: get_a_free_unit
        integer :: get_a_free_unit
        integer :: iun,ier
        iun=get_a_free_unit()
        if(iun<0)return
        open(UNIT=iun,FILE='TEST.cfg',STATUS='OLD',
     %       ACTION='READ',IOSTAT=ier)
        if(ier<0) return
        read(UNIT=iun,IOSTAT=ier,FMT=*)NX,NY
        close(UNIT=iun)
        return
        end
        integer function get_a_free_unit()
        implicit none
        integer :: i
        character (len=16) :: access_mode
          get_a_free_unit=-1
          do i = 1 , 99  ! find an available unit number
            inquire(UNIT=i,ACCESS=access_mode)
            if(trim(access_mode) == 'UNDEFINED')then ! found
              get_a_free_unit = i
              exit
            endif
          enddo
        return
        end
