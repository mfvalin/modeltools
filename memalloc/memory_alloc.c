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

// idiot simple memory mamager

// initialize memory arena
ITEM *IsmmgrInitArena(void *in, int total_size){
  ITEM *ap = (ITEM *)in ;   // pointer to arena
  ITEM *hp ;                // pointer to data block
  int sz = total_size - 4;  // size of arena

  if(sz < 1024) return(NULL) ;                     // minimum size not met
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
  if(sz < 0) return(0);
  if(ap[-1] != 0 || ap[sz] != 0 || ap[sz+1] != sz) return(0);
  return(sz);
}

// is arena minimally valid ?
static int IsmmgrBlockValid(void *iap, void *ip){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = (ITEM *)ip ;
  int sza, sz ;

  sza = IsmmgrArenaValid(iap) ;
  if(sza <= 0) return(0);   // invalid arena

  sz = ap[-2] ;
  if(sz <= 0) return(0);                     // invalid size
  if( (p - ap) < 2 ) return(0) ;             // data block starts before arena
  if( (p - ap) + sz > sza - 4 ) return(0) ;  // end of data block beyond arena end

  if(p[sz+1] != sz)  return(0);              // inconsistent size
  if(p[-1] != p[sz]) return(0);              // inconsistent free/used marker
  if(p[-1] != 0xBEBEFADA && p[-1] != 0xCAFEDECA ) // invalid marker
  return(sz);
}

// check integrity of arena
int IsmmgrCheck(void *iap){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;

  return(0);
}

// debug dump of arena structure
void IsmmgrDump(void *iap){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;

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
