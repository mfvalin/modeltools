#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define MAX_LINE 2048

static int header_written = 0 ;
static FILE *outfile = NULL ;
char *filename = NULL ;

void write_header(char *header_tag){   // write listing header (once)
  char header[256] ;
  size_t nb, one ;
  FILE *out = NULL;

  if(header_written) return ;
  outfile = stdout ;

  if(filename != NULL) out = fopen(filename, "w") ;
  if(out      != NULL) outfile = out ; // cannot open file, will send output to stdout

  header_written = 1 ;
  if(header_tag[0] == '\0') return ;  // no title

  nb = snprintf(header, sizeof(header), "\n================================================== %s ==================================================\n\n", header_tag) ;
  one = 1 ;
  fwrite(header, nb, one, outfile) ;
}

// line "tagger" for r.run_in_parallel / r.launch_mpi
//
// Usage: program  min_lines line_tag title [file_name]
//
//        min_lines    : suppression threshold, no output when less than min_lines lines
//        line_tag     : string that will prefix each output line
//        title        : ======= title =======   header
//        file_name    : use this file instead of stdout for output  [OPTIONAL argument]

int main(int argc, char **argv){
  char text[MAX_LINE] ;  // max length of expected output lines
  char *buf ;
  char *p ;
  int nb ;
  int maxbuf ;
  size_t one = 1 ;
  int minlines = 0 ;    // if less than minlines in input, there will be no output
  int nlines   = 0 ;    // lines counter in input
  char *line_tag = "" ;
  char *header_tag = "" ;

  if(argc > 1) minlines = atoi(argv[1]) ;   // minimum number of lines to avoid suppression [0]
  if(argc > 2) line_tag = argv[2] ;         // tag at start of lines [none]
  if(argc > 3) header_tag = argv[3] ;       // title header
  if(argc > 4) filename = argv[4] ;         // file to receive output                       [stdout]

  maxbuf = (((minlines > 0) ? minlines : 1)+1) * MAX_LINE ;
  p = buf = malloc( maxbuf * sizeof(char) ) ;  // allocate enough for (minlines + 1) lines of max length characters
  while( fgets(text, (MAX_LINE), stdin) != NULL ) {                 // read 1 line from stdin, until EOF
    nb = snprintf(p, (maxbuf) - (p - buf), "%s%s", line_tag, text) ; // copy it into temporary buffer
    p += nb ;
    nlines ++ ;
    if( (p - buf) > ((maxbuf) - (MAX_LINE)) || (nlines > minlines)) {     // less than MAX_LINE chars free in buffer
      write_header(header_tag) ;                                              // or min lines exceeded, write buffer to stdout
      fwrite(buf, p - buf, one, outfile);
      p = buf ;
    }
  }
}
