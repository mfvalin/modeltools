#define MAX_NAME     256
#define MAXWAFILES  1024
#define MAXPAGES      10


#define new_age_rd(age) (age+256)
#define new_age_wr(age) (age+512)
#define decay(age)      (age - (age >> 2))
#define LLSK long long
#define LSEEK lseek64
#define WSEEK(fd,offst,posi)\
 {\
  LLSK local_off;\
  local_off = offst;\
  LSEEK(fd,local_off * sizeof(int32_t),posi);\
 }

#define CMCARC_SIGN "CMCARCHS"  /* signature du debut d'un fichier cmcarc */
#define CMCARC_SIGN_V5 "CMCARCH5"  /* signature du debut d'un fichier cmcarc version 5 */

typedef struct {
   uint32_t *page_adr;
   uint64_t wa0;
   uint64_t walast;
   int access_count;
   int last_access;
   int touch_flag;
   int not_used_pad_for_word_alignment;
   } PAGEINFO;

typedef struct {
   int file_desc;
   int nb_page_in_use;
//   PAGEINFO page[MAXPAGES];
   PAGEINFO *page;             // ready for dynamic number of pages
   uint64_t offset;
   uint64_t segments[3];       // segment table, segment i addresses go from segments[i] -> segments[i+1]-1 0 means inactive
   } FILEINFO;

typedef struct {
   unsigned char ntc[4];        /* nt (longueur totale du fichier) en unites 64 bits */
   unsigned char ndc[4];        /* nd (longueur des donnees) en unites 64 bits */
   char cmcarc_name[MAX_NAME];
   } ENTETE_CMCARC;

typedef struct {
   unsigned char ntc[8];        /* nt (64 bits) (longueur totale du fichier) en unites 64 bits */
   unsigned char ndc[8];        /* nd (64 bits) (longueur des donnees) en unites 64 bits */
   char cmcarc_name[MAX_NAME];
   } ENTETE_CMCARC_V5;
   
