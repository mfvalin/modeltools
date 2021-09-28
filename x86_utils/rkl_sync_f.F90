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
    procedure, PASS(this) :: init       => InitLock
    procedure, PASS(this) :: initp      => InitLockP
    procedure, PASS(this) :: lock       => SetIdLock1, SetLock1              ! implies a fence
    procedure, PASS(this) :: lockf      => SetIdLock, SetLock                ! explicit fence argument
    procedure, PASS(this) :: trylock    => TrySetIdLock1, TrySetLock1        ! implies a fence
    procedure, PASS(this) :: trylockf   => TrySetIdLock, TrySetLock          ! explicit fence argument
    procedure, PASS(this) :: unlock     => ClearIdLock1, ClearLock1          ! implies a fence
    procedure, PASS(this) :: unlockf    => ClearIdLock, ClearLock            ! explicit fence argument
    procedure, PASS(this) :: tryunlock  => TryClearIdLock1, TryClearLock1    ! implies a fence
    procedure, PASS(this) :: tryunlockf => TryClearIdLock, TryClearLock      ! explicit fence argument
    procedure, PASS(this) :: valid      => ValidLock
    procedure, PASS(this) :: owner      => OwnerOfLock
    procedure, PASS(this) :: reset      => ResetOfLock
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

  function TrySetLock1(this) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT) :: status
    status = TryAcquireLock_(this%p, 1)
  end function TrySetLock1
  function TrySetLock(this, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    logical, intent(IN), value :: fence
    integer(C_INT) :: status
    if(fence) status = TryAcquireLock_(this%p, 1)
    if(.not. fence) status = TryAcquireLock_(this%p, 0)
  end function TrySetLock
  function TrySetIdLock1(this, id) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT) :: status
    status = TryAcquireIdLock_(this%p, id, 1)
  end function TrySetIdLock1
  function TrySetIdLock(this, id, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    logical, intent(IN), value :: fence
    integer(C_INT) :: status
    if(fence) status = TryAcquireIdLock_(this%p, id, 1)
    if(.not. fence) status = TryAcquireIdLock_(this%p, id, 0)
  end function TrySetIdLock

  subroutine SetLock1(this)
    class(rkl_lock), intent(IN) :: this
    call AcquireLock_(this%p, 1)
  end subroutine SetLock1
  subroutine SetLock(this, fence)
    class(rkl_lock), intent(IN) :: this
    logical, intent(IN), value :: fence
    if(fence) call AcquireLock_(this%p, 1)
    if(.not. fence) call AcquireLock_(this%p, 0)
  end subroutine SetLock
  subroutine SetIdLock1(this, id)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    call AcquireIdLock_(this%p, id, 1)
  end subroutine SetIdLock1
  subroutine SetIdLock(this, id, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    logical, intent(IN), value :: fence
    if(fence) call AcquireIdLock_(this%p, id, 1)
    if(.not. fence) call AcquireIdLock_(this%p, id, 0)
  end subroutine SetIdLock

  subroutine ResetOfLock(this)
    class(rkl_lock), intent(IN) :: this
    call ResetLock_(this%p)
  end subroutine ResetOfLock

  function TryClearLock1(this) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT) :: status
    status = TryReleaseLock_(this%p, 1)
  end function TryClearLock1
  function TryClearLock(this, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    logical, intent(IN), value :: fence
    integer(C_INT) :: status
    if(fence) status = TryReleaseLock_(this%p, 1)
    if(.not. fence) status = TryReleaseLock_(this%p, 0)
  end function TryClearLock
  function TryClearIdLock1(this, id) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    integer(C_INT) :: status
    status = TryReleaseIdLock_(this%p, id, 1)
  end function TryClearIdLock1
  function TryClearIdLock(this, id, fence) result(status)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    logical, intent(IN), value :: fence
    integer(C_INT) :: status
    if(fence) status = TryReleaseIdLock_(this%p, id, 1)
    if(.not. fence) status = TryReleaseIdLock_(this%p, id, 0)
  end function TryClearIdLock

  subroutine ClearLock1(this)
    class(rkl_lock), intent(IN) :: this
    call ReleaseLock_(this%p, 1)
  end subroutine ClearLock1
  subroutine ClearLock(this, fence)
    class(rkl_lock), intent(IN) :: this
    logical, intent(IN), value :: fence
    if(fence) call ReleaseLock_(this%p, 1)
    if(.not. fence) call ReleaseLock_(this%p, 0)
  end subroutine ClearLock
  subroutine ClearIdLock1(this, id)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    call ReleaseIdLock_(this%p, id, 1)
  end subroutine ClearIdLock1
  subroutine ClearIdLock(this, id, fence)
    class(rkl_lock), intent(IN) :: this
    integer(C_INT), intent(IN), value :: id
    logical, intent(IN), value :: fence
    if(fence) call ReleaseIdLock_(this%p, id, 1)
    if(.not. fence) call ReleaseIdLock_(this%p, id, 0)
  end subroutine ClearIdLock
end module
