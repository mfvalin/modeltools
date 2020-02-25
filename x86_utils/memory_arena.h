#include <stdint.h>

typedef struct{              // only used to calculate size of constant part of arena header (see memory_arena)
  uint32_t lock;
  uint32_t owner;
  uint32_t max_entries;
  uint32_t first_free;
  uint32_t n_entries;
  uint32_t arena_size;
} arena_header;
#define ArenaHeaderSize64 (sizeof(arena_header) / sizeof(uint64_t))

typedef struct{
  uint32_t lock;             // to lock this memory block
  uint32_t flags;
  uint32_t data_index;       // index relative to start of memory arena
  uint32_t data_size;        // size of data portion of block (64 bit units)
  uint64_t data_name;        // block name (max 8 characters)
} symtab_entry;
#define SymtabEntrySize64 (sizeof(symtab_entry) / sizeof(uint64_t))

#define MAX_SYMS 1000000
typedef struct{              // MUST BE CONSISTENT WITH arena_header
  uint32_t lock;             // to lock this memory arena
  uint32_t owner;            // MPI rank or PID of owner process
  uint32_t max_entries;      // max number of entries in t[]
  uint32_t first_free;       // index of first free location in arena
  uint32_t n_entries;        // number of entries in use in t[]
  uint32_t arena_size;       // size of memory arena (data + metadata) (64 bit units)
  symtab_entry t[MAX_SYMS];  // will never get allocated at that size, only used for indexing
} memory_arena;

typedef struct{
  uint64_t arena_name;       // name of segment (max 8 characters)
  size_t   arena_sz;         // size of segment
  int      arena_id;         // shared memory id of shared memory segment
  int      owner_id;
}master_entry;               // one entry per memory arena

#define MAX_MASTER 256
typedef struct{
  uint32_t lock;             // to lock master arena
  int      arena_id;         // shared memory id of master arena
  uint64_t arena_name;       // name of master arena  (max 8 characters)
  size_t   arena_sz;         // size of master arena segment
  master_entry me[MAX_MASTER];
} master_header;
#define MasterHeaderSize64 (sizeof(master_header) / sizeof(uint64_t))

typedef struct{              // MUST BE CONSISTENT WITH master_header
  uint32_t lock;             // to lock master arena
  int      arena_id;         // shared memory id of master arena
  uint64_t arena_name;       // name of master arena
  size_t   arena_sz;         // size of master arena segment
  master_entry me[MAX_MASTER];
  memory_arena ma;           // a normal memory arena
} master_arena;              // master arena contains the master table, followed by a normal arena

typedef struct{
  uint64_t arena_name;       // same as in associated master arena table
  size_t   arena_sz;         // same as in associated master arena table
  memory_arena *ma;          // pointer to memory arena in process space
}local_entry;                // one entry per memory arena

typedef struct{
  uint32_t     lock;         // should not be necessary
  int          master_id;    // shared memory id of master arena
  size_t       master_sz;    // size of segment
  master_arena *MA;          // pointer to master arena
  local_entry  le[MAX_MASTER];
}local_arena;                // copy in local process memory pointing to memory arenas

typedef struct{
  uint32_t fwd;              // forward index to next block (64 bit units) (0 for last block)
  uint32_t ix;               // index to this block (64 bit units)
  uint32_t nwd;              // length of data portion of block (in 64 bit units)
  uint32_t sign;             // low marker signature
}block_header;
#define BlockHeaderSize64 (sizeof(block_header) / sizeof(uint64_t))

typedef struct{
  uint32_t sign;             // high marker signature
  uint32_t bwd;              // backward index to start of this block (in 64 bit units)
}block_tail;
#define BlockTailSize64 (sizeof(block_tail) / sizeof(uint64_t))
