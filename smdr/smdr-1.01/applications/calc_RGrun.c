/* 
   This program implements renormalization group running of the
   Standard Model MSbar parameters from one specified scale to
   another.  It takes one required argument Qfinal and up to three
   optional arguments:

     Q_final               This required command line argument is the 
                           destination RG scale for the running.

     -i <input_filename>   Reads model data from <input_filename>; if
                           not specified, data are read from the file
                           ReferenceModel.dat, which should be in the
                           current directory.

     -int                  Runs the program in interactive mode. For each
                           input parameter, the user is prompted to either
                           accept the default value, or enter a new value.
                           The default values are either the ones from
                           the input file specified by the -i flag, or else
                           the ones found in ReferenceModel.dat.

     -l <loopOrder>        Runs the RG equations at the specified loop order.
                           Allowed values for loopOrder are 1,2,3,4,5.
                           If not specified, by default loopOrder = 5.
     
   As output, the program prints the values of the MSbar quantities at
   the scales Q = Q_in and Q_final. General and timing information for
   the computation are also printed.

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

      gcc -o calc_RGrun calc_RGrun.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and lib3vil.a and libtsil.a
   and the header file smdr.h are all in the current directory.

   Run as, for example:

      ./calc_RGrun 2.4e18

   for full 5-loop running and input data in ReferenceModel.dat, or

      ./calc_RGrun 2.4e18 -l 2 -i myInFile -int

   for 2-loop running and model data in myInFile, or

      ./calc_RGrun 2.4e18 -i myInFile -int

   for prompts that allow interactive modification of the initial
   parameter values.
*/

#include "smdr.h"

int main (int argc, char *argv[])
{
  char inputFile[50];
  SMDR_REAL Qfinal;
  int loopOrder;
  int interactiveMode;
  /* char funcname[] = "calc_RGrun"; */

  int nVars = 17;
  char *varList[] = {
    "SMDR_Q_in",
    "SMDR_lambda_in",
    "SMDR_g3_in",
    "SMDR_g_in",
    "SMDR_gp_in",
    "SMDR_yt_in",
    "SMDR_yb_in",
    "SMDR_yc_in",
    "SMDR_ys_in",
    "SMDR_yu_in",
    "SMDR_yd_in",
    "SMDR_ytau_in",
    "SMDR_ymu_in",
    "SMDR_ye_in",
    "SMDR_v_in",
    "SMDR_m2_in",
    "SMDR_Lambda_in"
  };

  /* Suppress warning messages */
  /* fclose (stderr); */

  /* Define arguments: */
  int nargs = 4;
  char *arglist[] = {"req", "-l","-i","-int"};
  char *argtype[] = {"real","int","string","toggle"};
  void *argvar[] = {&Qfinal, &loopOrder, inputFile, &interactiveMode};

  /* Set default values for optional args: */
  loopOrder = 5;
  strcpy (inputFile, "ReferenceModel.dat");
  interactiveMode = 0;

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  SMDR_Read_Values (inputFile, nVars, varList);

  if (interactiveMode)
    SMDR_Set_Values_Interactively (nVars, varList);

  SMDR_Load_Inputs ();

  printf("\nINITIAL MODEL PARAMETERS");
  if (0 == interactiveMode) {
    printf(" from \"%s\"", inputFile);
  }
  printf(":\n");

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();
  SMDR_Display_Lambda ();

  SMDR_Start_Timer();

  SMDR_RGeval_SM (Qfinal, loopOrder);

  SMDR_Timer ();

  printf("\nFINAL MODEL PARAMETERS:\n");
  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();
  SMDR_Display_Lambda ();

  printf("\nTotal time: %.2f seconds\n\n", SMDR_Time_Total);

  return 0;
}
