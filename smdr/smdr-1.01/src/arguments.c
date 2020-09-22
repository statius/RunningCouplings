/*
   This is a simple system for processing a set of command line
   arguments. The arguments and their types are defined in the calling
   program, and passed to SMDR_Process_Arguments () via the variables:

   nargs
   -----
   The total number of defined arguments, equal to the sum of the
   number of required arguments and the number of possible optional
   arguments.

   arglist 
   -------
   A list of strings (of dimension nargs) specifying the
   arguments. *Required* arguments must come first, and are given in a
   specified order that cannot be modified. A required argument is
   indicated with "req". After all required args are given, the
   arglist entries should define a set of flags that begin with
   '-'. Each of these specifies an *optional* argument. As an example,
   if there are three required arguments and two optional ones, then
   arglist might look like this:

   char arglist[] = {"req","req","req","-a","-i"};

   In this case the two optional argw would be specified on the
   command line as

   [-a <arg1> [-i <arg2>]]

   Optional arguments can be given in any order on the commad line,
   but the flags specified here must match up with the correspondng
   elements in the following arrays that define them.
   
   argtype
   -------
   An array of strings (of dimension nargs) defining the types of the
   arguments, in order. Allowed type strings are currently "real",
   "int", "string", and "toggle". As an example,

   char argtype[] = {"string","int","real","real","string"};

   Type "toggle" is used for an argument that does not have an
   associated value. If specified, the associated value is set to 1
   (TRUE).

   argvar
   ------
   This is an array of pointers-to-void (of dimension nargs), to data
   objects in the calling program to which the read values should be
   assigned. As an example,

   void argvar[] = {inFile, &i, &x, &Q, outFile};

   would assign the first argument in the variable inFile, the second
   in the variable i, and so on. Note that no ampersand is required
   for "string"-type variables.

   For arguments of type "toggle", the assigned variable should
   normally be of type int.  The effect of including the toggle will
   be to set this int to 1 (YES).

   As another example, for three optional arguments corresponding to
   an input filename, and output filename, and a real error tolerance,
   one could specify:

   int nargs = 3;
   char *arglist[] = {"-e","-i","-o"};
   char *argtype[] = {"real","string","string"};
   void *argvar[] = {&ERROR_TOLERANCE, inputFile, outputFile};

   Users should be sure to specify default values for all optional
   arguments in the calling program.

   Some simple tests for errors are included, but these are probably
   not bulletproof. Caveat emptor!
*/

#include "smdr_internal.h"

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* The assignment of an argument to the specified variable, in one
   convenient location to allow for extensions if desired.
*/

void SMDR_Assign_Argument (const char *arg, const char *atype, void *avar)
{
  char funcname[] = "SMDR_Assign_Argument";

  if (!strcmp (atype, "real"))
    *((SMDR_REAL *) avar) = (SMDR_REAL) strtold (arg, NULL);
  else if (!strcmp (atype, "int"))
    *((int *) avar) = (int) atoi (arg);
  else if (!strcmp (atype, "string"))
    strcpy (avar, arg);
  else if (!strcmp (atype, "toggle"))
    *((int *) avar) = TRUE;
  else {
    printf ("Type \"%s\" not recognized!\n", atype);
    SMDR_Error (funcname, "Error in argument type specification.", 876);
  }

  return;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Simplified errorr reporting without stack trace or calling function
   name, for use in argument processing.

   Prints to stdout so that true errors will always be seen. 
*/

void SMDR_ArgError (char *msg, int code)
{
  fprintf(stdout, "ERROR: %s Exiting...\n", msg);
  exit(code);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

int SMDR_Process_Arguments (int   argc, 
			    char *argv[], 
			    int   nargs,
			    char *arglist[],
			    char *argtype[],
			    void *argvar[])
{
  int i, j, k;
  /* char funcname[] = "SMDR_Process_Arguments"; */
  int numRequiredArgs;
  char errmsg[100], foo[5];

  /* First, find the number of required arguments: */
  numRequiredArgs = 0;
  for (i=0; i<nargs; i++)
    if (!strcmp (arglist[i], "req"))
      numRequiredArgs++;

  /* Then, a simple check: */
  if (numRequiredArgs > argc - 1)
    SMDR_ArgError ("Some required argument(s) is/are missing.", 876);

  i = 1;
  while (i < argc) {

    j = i - 1;
    
    if (!strcmp (arglist[j], "req")) {
      
      /* Simple test for an optional argument. */
      if (!strncmp (argv[i], "-", 1)) {
	sprintf (foo, "%d", i);
	strcpy (errmsg, "Argument ");
	strcat (errmsg, foo);
	strcat (errmsg, " is required but has been specified as optional.");
	SMDR_ArgError (errmsg, 876);
      }

      SMDR_Assign_Argument (argv[i], argtype[j], argvar[j]);
      i++;
    }
    else
      break;
  }

  /* Done with required args, now for the optional ones... */
  while (i < argc) {

    if (!strncmp (argv[i], "-", 1))
      ;
    else {
      strcpy (errmsg, "\"");
      strcat (errmsg, argv[i]);
      strcat (errmsg, "\" is not an optional argument flag, something is wrong.");
      SMDR_ArgError (errmsg, 876);
    }

    for (j=numRequiredArgs; j<nargs; j++) {
      if (!strcmp (argv[i], arglist[j])) {
	SMDR_Assign_Argument (argv[i+1], argtype[j], argvar[j]);
	break;
      }
    }

    /* If we get here, the flag was not recognized: */
    if (j == nargs) {
      strcpy (errmsg, "Optional argument \"");
      strcat (errmsg, argv[i]);
      strcat (errmsg, "\" not recognized.");
      SMDR_ArgError (errmsg, 876);
    }

    /* Go to the next argument */
    i++;

    /* If the flag was *not* a toggle, also step over the argument: */
    if (strcmp(argtype[j],"toggle"))
      i++;
  }
  
  return 0;
}
