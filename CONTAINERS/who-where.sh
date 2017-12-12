#!/bin/bash
export MP_CHILD=${MP_CHILD:-${PMI_RANK:-${OMPI_COMM_WORLD_RANK:-${ALPS_APP_PE}}}}
echo "CHILD = $MP_CHILD, $(taskset -cp $$) on $(hostname)"
[[ -d "$FLAGDIR" ]] && touch $FLAGDIR/FLAG_$MP_CHILD
true
