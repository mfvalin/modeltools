#include <stdint.h>

typedef struct{
  uint32_t lock;
  uint32_t flags;
  uint32_t max_entries;
  uint32_t first_free;
  uint32_t n_entries;
  uint32_t arena_size;
} arena_header;
#define ArenaHeaderSize64 (sizeof(arena_header) / sizeof(uint64_t))

typedef struct{
  uint32_t lock;
  uint32_t flags;
  uint32_t data_index;
  uint32_t data_size;
  uint64_t data_name;
} symtab_entry;
#define SymtabEntrySize64 (sizeof(symtab_entry) / sizeof(uint64_t))

typedef struct{
  uint32_t fwd;
  uint32_t ix;
  uint32_t nwd;
  uint32_t sign;
}block_header;
#define BlockHeaderSize64 (sizeof(block_header) / sizeof(uint64_t))

typedef struct{
  uint32_t sign;
  uint32_t bwd;
}block_tail;
#define BlockTailSize64 (sizeof(block_tail) / sizeof(uint64_t))
