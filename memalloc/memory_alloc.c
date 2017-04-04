#include <stdio.h>
#include <stdlib.h>

/*
 * memory arena structure
 * 
 *                |<--- data blocks --->|
 * +------+-------+------         ------+-------+------+
 * |  SZA |   0   | SZA          items  |   0   |  SZA |
 * +------+-------+------         ------+-------+------+
 *   [-2]   [-1]    [0]                   [SZA]  [SZA+1]
 * 
 * 
 * 2 consecutive data blocks
 * |<------------------------ block 1 ------------------------->|<------------------------- block 2 ------------------------>|
 * +---------+----------+------    --------+----------+---------+---------+----------+------    --------+----------+---------+
 * |   NW1   |  marker  |  NW1 data items  |  marker  |   NW1   |   NW2   |  marker  |  NW2 data items  |  marker  |   NW2   |
 * +---------+----------+------    --------+----------+---------+---------+----------+------    --------+----------+---------+
 *    [-2]       [-1]    [0]                  [NW1]     [NW1+1]    [-2]       [-1]    [0]                  [NW1]     [NW1+1]
 * 
 * same data blocks, fused (NW = NW1 + nW2 +4)
 * 
 * +---------+----------+------------------------------------        -----------------------------------+----------+---------+
 * |   NW    |  marker  |  NW data items                                                                |  marker  |   NW    |
 * +---------+----------+------------------------------------        -----------------------------------+----------+---------+
 *    [-2]       [-1]    [0]                                                                                [NW]      [NW+1]
 * 
 * sza    : size of useful memory arena
 * marker : 0xCAFEDECA   (block contains data)
 *          0xBEBEFADA   (block is empty)
 * NW     : size of block
 * item   : a 32 bit quantity (markers and NW are of type ITEM)
 * []     : indices used for addressing [0] at start of data
 * 
 */

typedef unsigned int ITEM ;
#define MINARENA 1024
#define MINBLOCK 16

// Idiot simple memory mamager (I s s mgr)

// initialize memory arena and create one data block to fill it entirely
// total_size is in ITEM units (see typedef above)
// returns arena pointer to be passed to other routines
ITEM *IsmmgrInitArena(void *in, int total_size){
  ITEM *ap = (ITEM *)in ;   // pointer to arena
  ITEM *hp ;                // pointer to data block
  int sz = total_size - 4;  // useful size of arena

  if(sz < MINARENA) return(NULL) ;                 // minimum size not met
  sz &= 0x7FFFFFFE ;                               // force even size

  // initialize arena
  ap += 2;  // from now on use arena indexing (data relative)
  ap[-2]   = sz ;          // arena size
  ap[-1]   = 0 ;           // pseudo size (like a data bkock)
  ap[sz]   = 0 ;           // pseudo size (like a data bkock)
  ap[sz+1] = sz ;          // arena size
  
  // create first data block in arena (empty block filling the whole arena)
  hp = ap + 2 ;            // a normal data block pointer
  sz = sz - 4 ;            // there will be sz -4 data items in single block occupying the whole arena
  hp[-1] = 0xBEBEFADA ;    // empty block marker
  hp[-2] = sz ;            // data block size (head)
  hp[sz] = 0xBEBEFADA ;    // empty block marker
  hp[sz+1] = sz ;          // data block size (tail)

  return(ap) ;                  // return data relative arena pointer that will be used from now on
}

// is arena minimally valid ? (return arena size if OK, 0 if not)
// iap  : arena pointer from IsmmgrInitArena
// returns arena size if OK, 0 otherwise
static int IsmmgrArenaValid(void *iap){
  ITEM *ap = (ITEM *)iap ;
  int sza ;

  sza = ap[-2] ;                   // arena size
  if(sza < MINARENA) return(0);    // unbelievably small, error
  if(ap[-1] != 0 || ap[sza] != 0 || ap[sza+1] != sza) return(0);   // inconsistent markers
  return(sza);                     // looks like an arena, return size
}

// are block and arena valid ?
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns block size if OK, 0 otherwise
static int IsmmgrBlockValid(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sza, sz ;

  sza = IsmmgrArenaValid(iap) ;
  if(sza == 0) return(0);                    // invalid arena

  sz = p[-2] ;                               // data block size
  if(sz == 0) return(sz);                    // invalid size
  if( (p - ap) < 2 ) return(0) ;             // data block starts before arena beginning
  if( (p - ap) + sz > sza - 2 ) return(0) ;  // end of data block beyond arena end

  if(p[sz+1] != sz)  return(0);              // inconsistent size
  if(p[-1] != p[sz]) return(0);              // inconsistent free/used marker
  if(p[-1] != 0xBEBEFADA && p[-1] != 0xCAFEDECA ) return(0); // invalid marker

  return(sz);                                // everything OK, return data block size
}

// is block valid and has size > 0 ? (arena is assumed valid)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns block size if OK, 0 otherwise
int IsmmgrBlockGood(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sza, sz ;

  sza = ap[-2] ;
  sz = p[-2] ;
  if(sz == 0) return(0);                      // invalid block size
  if( (p - ap) < 2 ) return(0) ;              // data block starts before arena
  if( (p - ap) + sz > sza - 2 ) return(0) ;   // end of data block beyond arena end

  if(p[sz+1] != sz)  return(0);               // inconsistent size
  if(p[-1] != p[sz]) return(0);               // inconsistent free/used marker
  if(p[-1] != 0xBEBEFADA && p[-1] != 0xCAFEDECA ) return(0); // invalid marker

  return(sz);
}

// get free block with best matching size
// from_top == 0 , scan from bottom end of memory arena
// from_top != 0 , scan from top end of memory arena
// iap  : arena pointer from IsmmgrInitArena
// size : desired minimum free size
// returns pointer to data block if found, NULL otherwise
static ITEM *IsmmgrBestMatch(void *iap, int size, int from_top){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;
  ITEM *px = NULL;
  int sza, sz ;
  int i0, match ;
  char *msg;

  sza = IsmmgrArenaValid(iap) ;
  if(sza <= 0) return (NULL);     // invalid arena

  match = 0x7FFFFFFF ;   // huge as initial best match size
  if(from_top){
    i0 = sza - ap[sza-1] -2 ; // start of top block
  }else{
    i0 = 2 ;                  // first block at ap[2]
  }
  sz = ap[i0-2] ;         // size of block ;
  while(sz > 0){          // until end of block chain reached
    p = &ap[i0] ;              // pointer to block ;
    sz = p[-2] ;
    if(p[-1] == 0xBEBEFADA){  // free block, try to match
      if( sz > size ){        // possible match
	if(sz <= size + MINBLOCK){  // we have a good enough match
	  px = p ;
	  break ;
	}
	if(sz < match){   // best match so far
	  match = sz ;
	  px = p;
	}
      }  // ( sz > size )
    }
    if(from_top){
      sz = p[-3] ;         // size of previous block
      i0 = i0 - sz - 4 ;   // start of previous block
    }else{
      i0 = i0 + sz + 4 ;   // start of next block
    }
// printf("i0 , sz : %d %d\n",i0,sz);
  }

  return(px);
}

// check integrity of arena (block chains)
// iap  : arena pointer from IsmmgrInitArena
// dump : if nonzero, print block metadata
// returns 0 if OK, nonzero otherwise
int IsmmgrCheck(void *iap, int dump){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;
  int sza, sz ;
  int i0 ;
  char *msg;

  sza = IsmmgrArenaValid(iap) ;
  if(sza == 0) return (-1);     // invalid arena
  i0 = 2 ;  // first block starts at ap[2], block occupies sz+4 items where sz is block data size
  while(i0 < sza){
    p = &ap[i0] ;
    sz = IsmmgrBlockGood(ap,p); // valid block ?
    if(dump) {
      msg = (p[-1] == 0xCAFEDECA) ? "USED" : "FREE" ;  // msg will be USED|FREE|INVALID
      if(sz == 0) msg = "INVALID" ;
      printf("[%10d] sz=%10d (%s)\n",i0,sz,msg);       // block index, block length, block status
    }
    if(sz == 0) break ;                                // end or bad block
    i0=i0+sz+4 ;
  }

  if(i0 < sza) printf("ERROR: premature end of block chain\n");
  return((i0 < sza) ? -1 : 0);
}

// set data block status to used 
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns 0 if OK, nonzero otherwise ( block invalid or alredy used)
int IsmmgrSetUsed(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz;

  if(IsmmgrBlockValid(ap,p) <= 0) return(-1);   // bad block
  if(p[-1] == 0xCAFEDECA) return(-1) ;          // block already used

  sz = p[-2] ;
  p[-1] = 0xCAFEDECA ;
  p[sz] = 0xCAFEDECA ;
  return(0);
}

// get block index relative to arena (0 if invalid)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns block index (in ITEM units) relative to ap if OK, 0 otherwise
// the minimum possible value is 2
int IsmmgrBlockIndex(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz = 0;

  if(IsmmgrBlockValid(ap,p) == 0) return(0);

  return(p - ap);
}

// get block size (negative if free block, 0 if invalid block, positice if used block)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns block size (in ITEM units) if OK, 0 otherwise
// if block is FREE, a negative size will be returned
// if block is USED, a positive size will be returned
int IsmmgrBlockSize(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz = 0;

  if(IsmmgrBlockValid(ap,p) == 0) return(0);

  if(p[-1] == 0xCAFEDECA) sz =  p[-2]  ;            // size of current used block
  if(p[-1] == 0xBEBEFADA) sz = -p[-2] ;             // size of current free block
  return(sz);
}

// get pointer to next block (NULL if no valid next block)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns pointer to block following this block (NULL if at end of chain or following block invalid)
ITEM *IsmmgrNextBlock(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz ;

  if( (sz = IsmmgrBlockGood(ap,p)) <= 0) return(NULL);

  p = p + sz + 4 ;         // point to next block
  return( (p[-2] > 0) ? p : NULL);   // return NULL if next block has size 0
}

// get pointer to previous block (NULL if no valid previous block)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns pointer to block preceding this block (NULL if at beginning of chain or preceding block invalid)
ITEM *IsmmgrPrevBlock(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz ;

  if(IsmmgrBlockValid(ap,p) == 0) return(NULL);

  sz = p[-3] ;       // size of previous block
  p = p - sz -4 ;    // point to previous block
  return( (sz > 0) ? p : NULL);   // return NULL if next block has size 0
}

// fuse adjacent data blocks if free (fuse next block then previous block if they are free)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns 0 if OK, nonzero otherwise
static int IsmmgrFuse(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p1 = (ITEM *)ip ;
  ITEM *p2 ;
  int s1, s2 ;

  s1 = IsmmgrBlockValid(ap,p1);   // get size of block (+validity check)
  if(s1 == 0) return (-1) ;       // invalid block
  if(p1[-1] != 0xBEBEFADA) return(-1) ; // this block is not free
  if(p1[s1 + 2] == 0) goto fuseprev ; // end of block chain, nothing to fuse forward
  p2 = &p1[s1 + 4] ;
  s2 = IsmmgrBlockValid(ap,p2);   // get size of next block
  if(s2 == 0) return (-1) ;       // invalid block

  if(p2[-1] == 0xBEBEFADA){       // next block is free, fuse with current
    p1[s1] = 0 ;                  // wipe upper p1 metadata
    p1[s1 + 1] = 0 ;
    p2[-1] = 0 ;                  // wipe lower p2 metadata
    p2[-2] = 0 ;
    s1 = s1 + s2 + 4;             // new size for p1
    p1[-2] = s1 ;                 // update lower p1 size
    p1[s1 + 1] = s1 ;             // update upper p1 size
  }
fuseprev:
  s2 = p1[-3] ;                   // size of brevious block ;
  if(s2 == 0) ;                   // beginning of block chain, nothing to fuse
  if(p1[-4] != 0xBEBEFADA) return(0) ; // previous block is not free
  p2 = &p1[-4 -s2] ;              // previous block
  s2 = IsmmgrBlockValid(ap,p2);   // get size of block (+validity check)
  if(s2 == 0) return (-1) ;       // invalid block

  p1[-1] = 0 ;                    // wipe lower p1 metadata
  p1[-2] = 0 ;
  p2[s2] = 0 ;                    // wipe upper p2 metadata
  p2[s2 + 1] = 0 ;
  s2 = s2 + s1 + 4 ;              // new size for p2
  p2[-2] = s2 ;                   // update lower p2 size
  p2[s2 + 1] = s2 ;               // update upper p1 size

  return(0) ;
}

// split a data block
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// size : >0 split with size ITEMs at beginning of block
// size : <0 split with size ITEMs at end of block
// returns pointer to leftovers data block if enough leftovers
// returns pointer to original block if not worth splitting
// returns NULL in case of error (invalid block, size larger than block size)
// splitting costs 4 ITEMs for new block metadata
static ITEM *IsmmgrSplit(void *iap, void *ip, int size){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p1 = (ITEM *)ip ;
  ITEM *p2 ;
  int sz, s1, s2 ;

  sz = IsmmgrBlockValid(ap,p1);   // get size

  if(sz ==0) return(NULL) ;      // bad block / arena, ERROR
  if(p1[-1] != 0xBEBEFADA) return(NULL) ;   // not a free blck, ERROR

  if(size == 0) return(p1) ;     // no split

  if(size >0) {   // size is low part
    s1 = (size+1) & 0x7FFFFFFE ; // force even size
    s2 = sz - s1 - 4 ;           // leftovers
    if(s2 < 0) return(NULL) ;    // sz too small, cannot split, ERROR
    if(s2 < MINBLOCK) return(p1) ;     // no split
  }else{          // size is high part
    size = -size ;
    s2 = (size+1) & 0x7FFFFFFE ; // force even size
    s1 = sz - s2 - 4 ;           // leftovers
    if(s1 < 0) return(NULL) ;    // sz too small, cannot split, ERROR
    if(s1 < MINBLOCK) return(p1) ;     // no split
  }

  p1[-2]   = s1 ;  // new size for block 1
  p1[s1+1] = s1 ;  // new size for block 1
  p1[s1]   = 0xBEBEFADA ; // free block flag

  p2 = &p1[s1+4] ;        // new block
  p2[-1]   = 0xBEBEFADA ; // free block flag
  p2[-2]   = s2 ;  // new size for block 2
  p2[s2+1] = s2 ;  // new size for block 2

  return(p2) ; // return pointer to block 2
}

// allocate a data block of size sz in arena ap (returns pointer to data block)
// iap  : arena pointer from IsmmgrInitArena
// sz   : size  (in ITEM units) of desired new block
// from_top : if 0 allocate as close as possible to beginning of arena
//            if nonzero allocate as close as possible to end of arena
// returns pointer to new block (or NULL if request cannot be satisfied)
ITEM *IsmmgrMalloc(void *iap, int sz, int from_top){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;
  ITEM *p2 = NULL;
  int status ;

  if(IsmmgrArenaValid(ap) == 0) return(NULL) ; // bad arena

  p = IsmmgrBestMatch(ap, sz, from_top) ;  // locate a free block of size at least sz (try a best fit)
  if(p == NULL) return(p) ;                // match failed 

  if(from_top) {
    p2 = IsmmgrSplit(ap, p, -sz) ;         // split to get a block of size at least sz in upper part
  }else{
    p2 = IsmmgrSplit(ap, p,  sz) ;         // split to get a block of size at least sz in lower part
    p2 = p ;
  }
  status = IsmmgrSetUsed(ap,p2) ;          // set block status to used
  return(p2);
}

// free data block p from arena ap (returns 0 upon success, != 0 otherwise)
// iap  : arena pointer from IsmmgrInitArena
// ip   : pointer to data block
// returns 0 if OK, nonzero otherwise
int IsmmgrFree(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz ;

  if(IsmmgrBlockValid(ap, p) == 0) return(-1) ; // bad block or arena
  sz = p[-2] ;
  if(p[-1] == 0xBEBEFADA) return(-1) ; // already free
  p[-1] = 0xBEBEFADA ;   // set marker to free
  p[sz] = 0xBEBEFADA ;   // set marker to free
  
  return (IsmmgrFuse(ap,p) );  // fuse with free blocks around freed block
}
#if defined(SELF_TEST)
main(int argc,char**argv){
  int arena[2048];
  int status, sza;
  ITEM *ap, *p1, *p2, *p3, *p4, *p5, *p6, *p0;
  int i1, i2, i3, i4, i5, i6;

  ap = IsmmgrInitArena(&arena[0], 2048);
  sza = ap[-2] ;
  printf("sza = %d\n",sza);
  printf("arena initialized, %8.8x %8.8x | %8.8x %8.8x \n",ap[-2],ap[-1],ap[sza],ap[sza+1]);
  printf("block initialized, %8.8x %8.8x | %8.8x %8.8x \n",ap[-0],ap[ 1],ap[sza-2],ap[sza-1]);
  
  status = IsmmgrArenaValid(ap);
  printf("status IsmmgrArenaValid = %d\n",status);

  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);

  p1 = &ap[2] ;
  printf("p1 = ap[%ld], size=%d\n",p1-ap,IsmmgrBlockSize(ap,p1)) ;
  printf("=============== split p1 at bottom (512) ===============\n");
  p2 = IsmmgrSplit(ap, p1, 512) ;
  printf("p1 = ap[%ld], size=%d\n",p1-ap,IsmmgrBlockSize(ap,p1)) ;
  printf("p2 = ap[%ld], size=%d\n",p2-ap,IsmmgrBlockSize(ap,p2)) ;
  printf("p2 does %s follow p1\n",(p2 == IsmmgrNextBlock(ap,p1)) ? "" : "not");
  printf("=============== split p2 at top (-128) ===============\n");
  p3 = IsmmgrSplit(ap, p2, -128) ;
  printf("p1 = ap[%ld], size=%d\n",p1-ap,IsmmgrBlockSize(ap,p1)) ;
  printf("p2 = ap[%ld], size=%d\n",p2-ap,IsmmgrBlockSize(ap,p2)) ;
  printf("p3 = ap[%ld], size=%d\n",p3-ap,IsmmgrBlockSize(ap,p3)) ;
  printf("p2 does %s precede p3\n",( p2 == IsmmgrPrevBlock(ap,p3) ) ? "" : "not");

  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  printf("IsmmgrSetUsed(ap,p2) = %d\n",IsmmgrSetUsed(ap,p2));
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);

  p0 = IsmmgrBestMatch(ap, 100, 0) ;   // from bottom
  printf("IsmmgrBestMatch(ap,100,0) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 100, 1) ;   // from top
  printf("IsmmgrBestMatch(ap,100,1) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 400, 0) ;   // from bottom
  printf("IsmmgrBestMatch(ap,400,0) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 400, 1) ;   // from top
  printf("IsmmgrBestMatch(ap,400,1) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 514, 0) ;   // from bottom
  printf("IsmmgrBestMatch(ap,514,0) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 514, 1) ;   // from top
  printf("IsmmgrBestMatch(ap,514,1) = %ld\n",(p0 ? p0 : ap)-ap);
  p2[-1] = 0xBEBEFADA ;  // artificially make p2 free
  p3[-4] = 0xBEBEFADA ;
  printf("artificially made p2(1392) free\n");
  p0 = IsmmgrBestMatch(ap, 514, 0) ;   // from bottom
  printf("IsmmgrBestMatch(ap,514,0) = %ld\n",(p0 ? p0 : ap)-ap);
  p0 = IsmmgrBestMatch(ap, 514, 1) ;   // from top
  printf("IsmmgrBestMatch(ap,514,1) = %ld\n",(p0 ? p0 : ap)-ap);

  printf("p2 bottom marker inconsistent\n");
  p2[-1] = 0xCAFEDECA;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  p2[-1] = 0xBEBEFADA;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  printf("p2 size inconsistent\n");
  p2[-2] = p2[-2] - 1;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);

  p2[-2] = p2[-2] + 1;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  printf("fuse p2 with neighbours\n");
  status = IsmmgrFuse(ap,p2);
  printf("status IsmmgrFuse = %d\n",status);
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);

  p1 = IsmmgrMalloc(ap, 512, 0) ; i1 = IsmmgrBlockIndex(ap,p1) ;
  p2 = IsmmgrMalloc(ap, 256, 0) ; i2 = IsmmgrBlockIndex(ap,p2) ;
  p3 = IsmmgrMalloc(ap, 128, 0) ; i3 = IsmmgrBlockIndex(ap,p3) ;
  p4 = IsmmgrMalloc(ap, 512, 1) ; i4 = IsmmgrBlockIndex(ap,p4) ;
  p5 = IsmmgrMalloc(ap, 256, 1) ; i5 = IsmmgrBlockIndex(ap,p5) ;
  p6 = IsmmgrMalloc(ap, 128, 1) ; i6 = IsmmgrBlockIndex(ap,p6) ;
  p0 = IsmmgrMalloc(ap, 224, 0) ;
  if(p0 == NULL) printf("allocate p0 failed as expected\n");
  if(p0 != NULL) printf("allocate p0 did not fail as expected\n");
  printf("p1 -> p6 = %5d %5d %5d %5d %5d %5d\n",i1, i2, i3, i4, i5, i6) ;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  status = IsmmgrFree(ap,p2) ; printf("free p2 \n");
  if(status) printf("free p2 failed\n");
  status = IsmmgrFree(ap,p5) ; printf("free p5 \n");
  if(status) printf("free p5 failed\n");
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  status = IsmmgrFree(ap,p3) ; printf("free p3 \n");
  if(status) printf("free p3 failed\n");
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  status = IsmmgrFree(ap,p4) ; printf("free p4 \n");
  if(status) printf("free p4 failed\n");
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  status = IsmmgrFree(ap,p6) ; printf("free p6 \n");
  if(status) printf("free p6 failed\n");
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
  status = IsmmgrFree(ap,p1) ; printf("free p1 \n");
  if(status) printf("free p1 failed\n");
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
}
#endif
