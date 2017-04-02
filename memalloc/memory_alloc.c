#include <stdio.h>
#include <stdlib.h>

/*
 * memory arena structure
 * 
 * +------+-------+------    ------+-------+------+
 * |  SZ  |   0   | SZ      items  |   0   |  SZ  |
 * +------+-------+------    ------+-------+------+
 *   [-2]   [-1]    [0]              [SZ]   [SZ+1]
 * 
 * 
 * split data blocks (2 consecutive blocks)
 * |<------------------------ block 1 ------------------------->|<------------------------- block 2 ------------------------>|
 * +---------+----------+------    --------+----------+---------+---------+----------+------    --------+----------+---------+
 * |   NW1   |  marker  |  NW1 data items  |  marker  |   NW1   |   NW2   |  marker  |  NW2 data items  |  marker  |   NW2   |
 * +---------+----------+------    --------+----------+---------+---------+----------+------    --------+----------+---------+
 *    [-2]       [-1]    [0]                  [NW1]     [NW1+1]    [-2]       [-1]    [0]                  [NW1]     [NW1+1]
 * 
 * fused data blocks (NW = NW1 + nW2 +4)
 * 
 * +---------+----------+------------------------------------        -----------------------------------+----------+---------+
 * |   NW    |  marker  |  NW data items                                                                |  marker  |   NW    |
 * +---------+----------+------------------------------------        -----------------------------------+----------+---------+
 *    [-2]       [-1]    [0]                                                                                [NW]      [NW+1]
 * 
 * marker : 0xCAFEDECA   (block contains data)
 *          0xBEBEFADA   (block is empty)
 * NW     : size of block
 * item   : a 32 bit quantity (markers and NW are unsigned integers)
 * []     : indices used for addressing [0] at start of data
 * 
 */

#define ITEM unsigned int
#define MINARENA 1024
#define MINBLOCK 16

// idiot simple memory mamager

// initialize memory arena
ITEM *IsmmgrInitArena(void *in, int total_size){
  ITEM *ap = (ITEM *)in ;   // pointer to arena
  ITEM *hp ;                // pointer to data block
  int sz = total_size - 4;  // size of arena

  if(sz < MINARENA) return(NULL) ;                 // minimum size not met
  sz &= 0x7FFFFFFE ;                               // force even size

  // initialize arena
  ap += 2;  // from now on use arena indexing
  ap[-2]   = sz ;
  ap[-1]   = 0 ;
  ap[sz]   = 0 ;
  ap[sz+1] = sz ;
  
  // create first data block in arena (empty block filling the whole arena)
  hp = ap + 2 ;            // a normal data block pointer
  sz = sz - 4 ;            // there will be sz -4 data items in single block occupying the whole arena
  hp[-1] = 0xBEBEFADA ;    // empty block marker
  hp[-2] = sz ;            // data block size (head)
  hp[sz] = 0xBEBEFADA ;    // empty block marker
  hp[sz+1] = sz ;          // data block size (tail)

  return(ap) ;                  // return arena pointer that will be used from now on
}

// is arena minimally valid ?
static int IsmmgrArenaValid(void *iap){
  ITEM *ap = (ITEM *)iap ;
  int sz ;

  sz = ap[-2] ;
  if(sz < MINARENA) return(0);
  if(ap[-1] != 0 || ap[sz] != 0 || ap[sz+1] != sz) return(0);
  return(sz);
}

// are block and arena valid ? (return size of block or 0)
static int IsmmgrBlockValid(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sza, sz ;

  sza = IsmmgrArenaValid(iap) ;
  if(sza <= 0) return(0);   // invalid arena

  sz = p[-2] ;
  if(sz == 0) return(sz);                    // invalid size
  if( (p - ap) < 2 ) return(0) ;             // data block starts before arena
  if( (p - ap) + sz > sza - 2 ) return(0) ;  // end of data block beyond arena end

  if(p[sz+1] != sz)  return(0);              // inconsistent size
  if(p[-1] != p[sz]) return(0);              // inconsistent free/used marker
  if(p[-1] != 0xBEBEFADA && p[-1] != 0xCAFEDECA ) return(0); // invalid marker

  return(sz);
}

// is block valid and has size > 0 ? (arena assumed valid) (return size of block or 0)
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

// check integrity of arena (print block metadata id dump != 0)
int IsmmgrCheck(void *iap, int dump){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;
  int sza, sz ;
  int i0 ;
  char *msg;

  sza = IsmmgrArenaValid(iap) ;
  if(sza <= 0) return (-1);     // invalid arena
  i0 = 2 ;  // first block starts at ap[2], bblocx occupies sz+4 items where sz is block data size
  while(i0 < sza){
    p = &ap[i0] ;
    sz = IsmmgrBlockGood(ap,p);
    if(dump) {
      msg = (p[-1] == 0xCAFEDECA) ? "USED" : "FREE" ;
      if(sz == 0) msg = "INV " ;
      printf("[%10d] sz=%10d (%s)\n",i0,sz,msg);    // block index, block length, block status
    }
    if(sz == 0) break ;
    i0=i0+sz+4 ;
  }

  if(i0 < sza) printf("ERROR: premature end of block chain\n");
  return((i0 < sza) ? -1 : 0);
}

// get block size (negative if free block, 0 if invalid block)
static int IsmmgrBlockSize(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz = 0;

  if(IsmmgrBlockValid(ap,p) <= 0) return(0);

  if(p[-1] == 0xCAFEDECA) sz =  p[-2]  ;            // size of current used block
  if(p[-1] == 0xBEBEFADA) sz = -p[-2] ;             // size of current free block
  return(sz);
}
// get pointer to next block (NULL if no next block)
static ITEM *IsmmgrNextBlock(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz ;

  if( (sz = IsmmgrBlockGood(ap,p)) <= 0) return(NULL);

  p = p + sz + 4 ;         // point to next block
  return( (p[-2] > 0) ? p : NULL);   // return NULL if next block has size 0
}

// get pointer to previous block (NULL if no previous block)
static ITEM *IsmmgrPrevBlock(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sz ;

  if(IsmmgrBlockValid(ap,p) <= 0) return(NULL);

  sz = p[-3] ;       // size of previous block
  p = p - sz -4 ;    // point to previous block
  return( (sz > 0) ? p : NULL);   // return NULL if next block has size 0
}

// split a data block
static ITEM *IsmmgrSplit(void *iap, void *ip, int size){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p1 = (ITEM *)ip ;
  ITEM *p2 ;
  int sz, s1, s2 ;

  sz = IsmmgrBlockValid(ap,p1);   // get size

  if(sz <=0) return(NULL) ;      // bad block / arena, ERROR
  if(p1[-1] != 0xBEBEFADA) return(NULL) ;   // not a free blck, ERROR

  if(size == 0) return(p1) ;     // no split

  if(size >0) {   // size is low part
    s1 = (size+1) & 0x7FFFFFFE ; // force even
    s2 = sz - s1 - 4 ;           // leftovers
    if(s2 < 0) return(NULL) ;    // sz too small, cannot split, ERROR
    if(s2 < MINBLOCK) return(p1) ;     // no split
  }else{          // size is high part
    size = -size ;
    s2 = (size+1) & 0x7FFFFFFE ; // force even
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

// allocate a data block of size sz in arena ap
ITEM *IsmmgrMalloc(void *iap, int sz){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;

  if(IsmmgrArenaValid(ap) <= 0) return(NULL) ; // bad arena
  // locate a hole of size at least sz (try a best fit)
  return(p);
}

// free data block p from arena ap
int IsmmgrFree(void *iap, void *ip, int sz){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;

  if(IsmmgrBlockValid(ap, p) <= 0) return(-1) ; // bad block or arena
  return(0);
}

main(int argc,char**argv){
  int arena[2048];
  int status, sza;
  ITEM *ap, *p1, *p2, *p3;

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
  p2 = IsmmgrSplit(ap, p1, 512) ;
  printf("p2 = ap[%ld], size=%d\n",p2-ap,IsmmgrBlockSize(ap,p2)) ;
  printf("p2 does %s follow p1\n",(p2 == IsmmgrNextBlock(ap,p1)) ? "" : "not");
  p3 = IsmmgrSplit(ap, p2, -128) ;
  printf("p3 = ap[%ld], size=%d\n",p3-ap,IsmmgrBlockSize(ap,p3)) ;
  printf("p2 does %s precede p3\n",( p2 == IsmmgrPrevBlock(ap,p3) ) ? "" : "not");

  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);

  p3[-2] = 0;
  status = IsmmgrCheck(ap,1);
  printf("status IsmmgrCheck = %d\n",status);
}
