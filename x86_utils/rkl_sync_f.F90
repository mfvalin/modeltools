! Copyright (C) 2021  Environnement et Changement climatique Canada
!
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! Author:
!     M. Valin,   Recherche en Prevision Numerique, 2021
module rkl_sync_mod
  use ISO_C_BINDING
  implicit none
  include 'rkl_sync.inc'

  private :: InitLock, InitLockP, SetIdLock, SetLock, ClearIdLock, ClearLock, ValidLock
  type :: rkl_lock
    private
    type(C_PTR) :: p = C_NULL_PTR  ! C pointer to address used for the lock
  contains
    procedure, PASS(this) :: init   => InitLock
    procedure, PASS(this) :: initp  => InitLockP
    procedure, PASS(this) :: lock   => SetIdLock, SetLock
    procedure, PASS(this) :: trylock   => TrySetIdLock, TrySetLock
    procedure, PASS(this) :: unlock => ClearIdLock, ClearLock
    procedure, PASS(this) :: tryunlock => TryClearIdLock, TryClearLock
    procedure, PASS(this) :: valid  => ValidLock
    procedure, PASS(this) :: owner  => OwnerOfLock
    procedure, PASS(this) :: reset  => ResetOfLock
  end type

  private :: InitBarrier, InitBarrierP, ValidBarrier
  type :: rkl_barrier
    private
    type(C_PTR) :: flag  = C_NULL_PTR   ! C pointer to address used for the barrier variable
    type(C_PTR) :: count = C_NULL_PTR   ! C pointer to address used for arrival count
  contains
    procedure, PASS(this) :: init  => InitBarrier
    procedure, PASS(this) :: initp => InitBarrierP
    procedure, PASS(this) :: valid => ValidBarrier
    procedure, PASS(this) :: set   => SetBarrier
  end type

  private :: set_to_zero_32

  contains

  subroutine set_to_zero_32(what, n)
    type(C_PTR), intent(IN), value :: what
    integer(C_INT), intent(IN), value :: n
    integer, dimension(:), pointer :: a
    call C_F_POINTER(what, a, [n])
    a(1:n) = 0
  end subroutine set_to_zero_32

  subroutine InitBarrier(this, flag, count)
     class(rkl_barrier), intent(INOUT) :: this
     integer(C_INT), intent(IN), target :: flag
     integer(C_INT), intent(IN), target :: count
     this%flag = C_LOC(flag)
     call set_to_zero_32(this%flag, 1)
     this%count = C_LOC(count)
     call set_to_zero_32(this%count, 1)
  end subroutine InitBarrier

  subroutine InitBarrierP(this, flag, count)
     class(rkl_barrier), intent(INOUT) :: this
     type(C_PTR), intent(IN), value :: flag
     type(C_PTR), intent(IN), value :: count
     this%flag = flag
     call set_to_zero_32(this%flag, 1)
     this%count = count
     call set_to_zero_32(this%count, 1)
  end subroutine InitBarrierP

  function ValidBarrier(this) result(valid)
     class(rkl_barrier), intent(IN) :: this
     logical :: valid
     valid = C_ASSOCIATED(this%flag) .and. C_ASSOCIATED(this%count)
  end function ValidBarrier

  subroutine SetBarrier(this, maxcount)
     class(rkl_barrier), intent(INOUT) :: this
     integer(C_INT), intent(IN) :: maxcount
     call BasicNodeBarrier_(this%flag, this%count, maxcount)
  end subroutine SetBarrier

  subroutine InitLock(this, what)
     class(rkl_lock), intent(INOUT) :: this
     integer(C_INT), intent(IN), target :: what
     this%p = C_LOC(what)
     call set_to_zero_32(this%p, 1)
  end subroutine InitLock

  subroutine InitLockP(this, what)
     class(rkl_lock), intent(INOUT) :: this
     type(C_PTR), intent(IN), value :: what
     this%p = what
     call set_to_zero_32(this%p, 1)
  end subroutine InitLockP

  function ValidLock(this) result(valid)
     class(rkl_lock), intent(IN) :: this
     logical :: valid
     valid = C_ASSOCIATED(this%p)
  end function ValidLock

  function OwnerOfLock(this) result(owner)
     class(rkl_lock), intent(IN) :: this
     integer(C_INT) :: owner
     owner = LockOwner_(this%p)
  end function OwnerOfLock

  function TrySetLock(this, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: fence
    integer(C_INT) :: status
    status = TryAcquireLock_(this%p, fence)
  end function TrySetLock
  function TrySetIdLock(this, id, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT), intent(IN), value :: fence
    integer(C_INT) :: status
    status = TryAcquireIdLock_(this%p, id, fence)
  end function TrySetIdLock

  subroutine SetLock(this, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: fence
    call AcquireLock_(this%p, fence)
  end subroutine SetLock
  subroutine SetIdLock(this, id, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT), intent(IN), value :: fence
    call AcquireIdLock_(this%p, id, fence)
  end subroutine SetIdLock

  subroutine ResetOfLock(this)
    class(rkl_lock), intent(IN) :: this
    call ResetLock_(this%p)
  end subroutine ResetOfLock

  function TryClearLock(this, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: fence
    integer(C_INT) :: status
    status = TryReleaseLock_(this%p, fence)
  end function TryClearLock
  function TryClearIdLock(this, id, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT), intent(IN), value :: fence
    integer(C_INT) :: status
    status = TryReleaseIdLock_(this%p, id, fence)
  end function TryClearIdLock

  subroutine ClearLock(this, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: fence
    call ReleaseLock_(this%p, fence)
  end subroutine ClearLock
  subroutine ClearIdLock(this, id, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT), intent(IN), value :: fence
    call ReleaseIdLock_(this%p, id, fence)
  end subroutine ClearIdLock
end module
