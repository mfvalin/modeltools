/* functions for C and FORTRAN programming
 * Copyright (C) 2019  M.Valin (UQAM)
 *
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <errno.h>

static int new_stdout, new_stderr, save_stderr, save_stdout;
static int saved = 0;

// redirect stdout to file newout and stderr to file newerr
// if flush is nonzero, flush stdout and stderr before redirecting
// providing the same name for newstdout and stderr is O.K.
// override environment variables REDIRECT_STDOUT and REDIRECT_STDERR have precendence over newout and newerr
int redirect_stdout_stderr(char *newout, char *newerr, int flush, int append){
    char *newstdout = getenv("REDIRECT_STDOUT");
    char *newstderr = getenv("REDIRECT_STDERR");
    int mode = O_RDWR|O_CREAT;

    if(newstdout == NULL) newstdout = newout;
    if(newstderr == NULL) newstderr = newerr;
    if(append) mode |= O_APPEND;

    new_stdout = open(newstdout, mode, 0600);        // new descriptor for stdout
    if (-1 == new_stdout) { perror("problem opening new stdout"); return 255; }

    new_stderr = open(newstderr, mode, 0600);        // new descriptor for stderr
    if (-1 == new_stderr) { perror("problem opening new_stderr"); return 255; }

    if(flush) fflush(stdout);      // flush buffers before redirecting if requested
    if(flush) fflush(stderr);      // flush buffers before redirecting if requested
    save_stdout = dup(fileno(stdout));
    save_stderr = dup(fileno(stderr));

    if (-1 == dup2(new_stdout, fileno(stdout))) { perror("problem redirecting stdout"); return 255; }
    if (-1 == dup2(new_stderr, fileno(stderr))) { perror("problem redirecting stderr"); return 255; }

    saved = 1;
    return 0;
}

// copy contents of fd onto stream, prefixing each text line with prefix string
static void prefix_lines(int fd, FILE *stream, char *prefix){
    off_t offset = 0;
    char line[4097];
    int nc, start, end, lastc;

    lseek(fd, offset, SEEK_SET);            // rewind text source
    nc = read(fd, line, sizeof(line)-1) ;   // read first buffer
    line[nc] = '\0';
    lastc = '\n';
    while(nc > 0){                          // while not end of file
      start = 0; end = 0;
      while(end < nc){                      // do not go beyond end of buffer
	while(end < nc && line[end] != '\n') end++;
	if(lastc == '\n') fprintf(stream,"\n%s",prefix);  // last character seen was newline, print prefix
	lastc = line[end];                  // null of lf (lf if end points ot end of line)
	line[end] = '\0';
	fprintf(stream,"%s",line+start);    // write line or portion of line
	end++;                              // next 
	start = end;
      }
      nc = read(fd, line, sizeof(line)) ;   // next buffer
    }
    fprintf(stream,"\n\n");
}

// restore stdout and stderr to their original destination
// print error message if no stdout/stderr was saved
// if prefixo is not NULL, copy captured stdout onto original stdout with prefixo at start of lines
// if prefixe is not NULL, copy captured stderr onto original stderr with prefixe at start of lines
// formato end formate are used for the title line when copying stdout/stderr onto original
// in an MPI context, formato/formate would be expected to contain info such as host name, MPI rank, ...
void restore_stdout_stderr(char *prefixo, char *formato, char *prefixe, char *formate){

    if(! saved) {
      errno = ENOTSUP;
      perror("no original stdout/stderr saved");
    }

    fflush(stdout); 
    fflush(stderr); 

    dup2(save_stdout, fileno(stdout));
    close(save_stdout);
    dup2(save_stderr, fileno(stderr));
    close(save_stderr);

    if(prefixo != NULL){
      fprintf(stdout, "%s%s", prefixo, formato);  // put out header line
      prefix_lines(new_stdout, stdout, prefixo);  // copy captured onto original
    }
    if(prefixe != NULL){
      fprintf(stderr, "%s%s", prefixe, formate);  // put out header line
      prefix_lines(new_stderr, stderr, prefixe);  // copy captured onto original
    }

    close(new_stdout);
    close(new_stderr);

    saved = 0;
}
#if defined(SELF_TEST)
int main(int argc, const char *argv[])   // idiot test program
{
    fprintf(stdout,"This should appear(1) argc = %d (stdout)\n", argc);
    fprintf(stderr,"This should appear(1) argc = %d (stderr)\n", argc);

    restore_stdout_stderr(NULL, "", NULL, "");   // deliberate error

    if( redirect_stdout_stderr("stdout.log", "stderr.log", argc > 1 ? atoi(argv[1]) : 1, 1 ) ) perror("redirection of stdout/stderr") ;

    fprintf(stdout,"1 some output after redirection (stdout)\n");
    fprintf(stderr,"1 some output after redirection (stderr)\n");
    fprintf(stdout,"2 some output after redirection (stdout)\n");
    fprintf(stderr,"2 some output after redirection (stderr)\n");

    restore_stdout_stderr("0000o: ", "==== captured stdout (0000) ===", "0000e: ", "==== captured stderr (0000) ===");

    fprintf(stdout,"This should appear(2) (stdout)\n");
    fprintf(stderr,"This should appear(2) (stderr)\n");
    fprintf(stdout,"back to original (stdout)\n");
    fprintf(stderr,"back to original (stderr)\n");

    return 0;
}
#endif
