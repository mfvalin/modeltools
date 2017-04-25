#include <stdio.h>

#pragma weak flush_listing_c_stream__=flush_listing_c_stream
#pragma weak flush_listing_c_stream_=flush_listing_c_stream
void flush_listing_c_stream__() ;
void flush_listing_c_stream_() ;
void flush_listing_c_stream() {   // flush listing streams, stderr and stdout
  fflush(stdout);
  fflush(stderr);
}