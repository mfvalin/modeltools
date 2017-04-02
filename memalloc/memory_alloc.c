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
 * item   : a 32 bit quantity (markers and NW are unsigned integers)
 * []     : indices used for addressing [0] at start of data
 * 
 */

#define ITEM unsigned int

// idiot simple memory mamager

// create first data block in arena (empty block filling the whole arena)
static void ismmgr_create_heap(void *in){
  ITEM *ap = (ITEM *)in ;
  ITEM *hp ;
  int sz ;

  sz = ap[-2] ;            // arena size
  hp = ap + 2 ;            // a normal data block pointer
  sz = sz - 4 ;            // there will be sz -4 data items in single block occupying the whole arena
  hp[-1] = 0xBEBEFADA ;    // empty block marker
  hp[-2] = sz ;            // block data size
  hp[sz] = 0xBEBEFADA ;    // empty block marker
  hp[sz+1] = sz ;
}

// initialize memory arena
ITEM *ismmgr_init_arena(void *in, int total_size){
  ITEM *ap = (ITEM *)in;
  int sz = total_size - 4;

  if(sz < 1024) return(NULL) ;                     // minimum size not met

  ap += 2;  // from now on use arena indexing
  ap[-2]   = sz ;
  ap[-1]   = 0 ;
  ap[sz]   = 0 ;
  ap[sz+1] = sz ;
  
  init_ismmgr_create_heap(ap) ; // create heap

  return(ap) ;                  // return arena pointer that will be used from now on
}

ITEM *ismmgr_malloc(void *iap, int sz){
  ITEM *ap = (ITEM *)iap ;
  ITEM *p = NULL;

  // locate a hole of size at least sz (try a best fit)
  return(p);
}
