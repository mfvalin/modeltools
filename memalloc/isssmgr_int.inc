      interface
!ITEM *IsmmgrInitArena(void *in, int total_size)
	function IsmmgrInitArena(in,total_size) result(p) bind(C,name='IsmmgrInitArena')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: in
	  integer(C_INT), intent(IN), value :: total_size
	  type(C_PTR) :: p
	end function IsmmgrInitArena
!int IsmmgrBlockValid(void *iap, void *ip)
	function IsmmgrBlockValid(iap,ip) result(s) bind(C,name='IsmmgrBlockValid')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrBlockValid
!int IsmmgrBlockGood(void *iap, void *ip)
	function IsmmgrBlockGood(iap,ip) result(s) bind(C,name='IsmmgrBlockGood')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrBlockGood
!int IsmmgrCheck(void *iap, int dump)
	function IsmmgrCheck(iap,dump) result(s) bind(C,name='IsmmgrCheck')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap
	  integer(C_INT), intent(IN), value :: dump
	  integer(C_INT) :: s
	end function IsmmgrCheck
!int IsmmgrSetUsed(void *iap, void *ip)
	function IsmmgrSetUsed(iap,ip) result(s) bind(C,name='IsmmgrSetUsed')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrSetUsed
!int IsmmgrBlockIndex(void *iap, void *ip)
	function IsmmgrBlockIndex(iap,ip) result(s) bind(C,name='IsmmgrBlockIndex')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrBlockIndex
!int IsmmgrBlockSize(void *iap, void *ip)
	function IsmmgrBlockSize(iap,ip) result(s) bind(C,name='IsmmgrBlockSize')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrBlockSize
!ITEM *IsmmgrNextBlock(void *iap, void *ip)
	function IsmmgrNextBlock(iap,ip) result(p) bind(C,name='IsmmgrNextBlock')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  type(C_PTR) :: p
	end function IsmmgrNextBlock
!ITEM *IsmmgrPrevBlock(void *iap, void *ip)
	function IsmmgrPrevBlock(iap,ip) result(p) bind(C,name='IsmmgrPrevBlock')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  type(C_PTR) :: p
	end function IsmmgrPrevBlock
!ITEM *IsmmgrMalloc(void *iap, int sz, int from_top)
	function IsmmgrMalloc(iap,sz,from_top) result(p) bind(C,name='IsmmgrMalloc')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap
	  integer(C_INT), intent(IN), value :: sz, from_top
	end function IsmmgrMalloc
!int IsmmgrFree(void *iap, void *ip)
	function IsmmgrFree(iap,ip) result(s) bind(C,name='IsmmgrFree')
	  import :: C_PTR, C_INT
	  type(C_PTR), intent(IN), value :: iap, ip
	  integer(C_INT) :: s
	end function IsmmgrFree
      end interface
