#include "smdr_internal.h"
#include <string.h>

/* Error and warning messages. */

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Warn (char *func, char *msg)
{
  fprintf (stderr, "WARN (%s): %s\n", func, msg);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Prints to stdout so that true errors will always be seen. */

#include <execinfo.h>

void SMDR_Error (char *func, char *msg, int code)
{
  void *callstack[128];
  int i, frames = backtrace(callstack, 128);
  char **strs = backtrace_symbols(callstack, frames);

  fprintf(stdout, "ERROR (%s): %s Exiting...\n", func, msg);
  for (i=0; i<frames; i++)
    fprintf(stdout, "%s\n", strs[i]);

  free(strs);
  exit(code);
}
