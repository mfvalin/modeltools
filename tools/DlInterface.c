#include <stdlib.h>
#include <dlfcn.h>
/*
  rmnlib interface to the dynamic loading functions
  if this module is compiled without -DLIVE
  it provides stubs that return a failure code when called
  the purpose is to have a default version in the library that does not
  necessitate -ldl at link time for applications
  if this module is compiled with -DLIVE
  it becomes a direct interface do dlopen/dlsym/dlerror/dlclose
*/
#define ERR_NOT_ACTIVE "ERROR: this is the dummy dynamic loader\n"

void *DlOpen(const char *filename, int flag)
{
#ifdef LIVE
  return(dlopen(filename,flag));
#else
  return(NULL);
#endif
}
void *DlSym(void *handle, const char *symbol)
{
#ifdef LIVE
  return (dlsym(handle,symbol));
#else
  return(NULL);
#endif
}
char *DlError(void)
{
#ifdef LIVE
  return (dlerror());
#else
  return(ERR_NOT_ACTIVE);
#endif
}
int DlClose(void *handle)
{
#ifdef LIVE
  return (dlclose(handle));
#else
  return(-1);
#endif
}
