/* 
   This program calculates the Higgs VEV, obtained by minimizing the
   Landau gauge effective potential with resummed Goldston boson
   contributions, using, by default, the full 3-loop results of
   arXiv:1709.02397 plus the 4-loop QCD contribution from 1508.00912.
   It takes one optional command-line inputs:

   -l <loopOrder>

   The allowed values for loopOrder are:

   0   [tree-level]
   1   [1-loop]
   2   [2-loop, from 1406.2355 eq. (4.20), or SMDeltas.anc of 1709.02397]
   2.5 [adds the leading 3-loop contributions from 1406.2355 eq. (4.21)]
   3   [full 3-loop, from ancillary file SMDeltas.anc of 1709.02397]
   3.5 [full 3-loop plus the 4-loop contribution at leading order in QCD,
        from 1508.00912 equation (5.5)]

   If -l <loopOrder> is omitted, it is taken to be 3.5 by default.  

   The program then prompts the user for the quantities:
     Q, m2, lambda, g3, g, gp, yt, yb, ytau 
   with default values taken from the file ReferenceModel.dat, which
   should be in the current directory. It then calculates and prints
   the full set of MSbar quantities including the value of the VEV v.

   Dimensionful quantities are in appropriate powers of GeV.  Timing
   information for the computation is also printed.

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o calc_vev calc_vev.c -L. -lsmdr -ltsil -l3vil -lm

   assuming the necessary archive files libsmdr.a and lib3vil.a and
   libtsil.a and the header file smdr.h are all in the current
   directory.

   Run as, for example:

   ./calc_vev 

   or 

   ./calc_vev -l 2 

   for the result including only the 2-loop contributions.
*/

#include "smdr.h"

int main (int argc, char *argv[])
{
  SMDR_REAL loopOrder;
  char inputFile[] = "ReferenceModel.dat";
  /* char funcname[] = "calc_vev"; */

  /* Input variables to be read in: */
  int nVars = 9;
  char *varList[] = {
     "SMDR_Q_in",
     "SMDR_m2_in",
     "SMDR_lambda_in",
     "SMDR_g3_in",
     "SMDR_g_in",
     "SMDR_gp_in",
     "SMDR_yt_in",
     "SMDR_yb_in",
     "SMDR_ytau_in"};

  /* Define arguments: */
  int nargs = 1;
  char *arglist[] = {"-l"};
  char *argtype[] = {"real"};
  void *argvar[] = {&loopOrder};

  /* Set default values for optional args: */
  loopOrder = 3.5;

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  SMDR_Read_Values (inputFile, nVars, varList);

  /* Set all 1st and 2nd family Yukawas to zero. */
  SMDR_yc_in = SMDR_ys_in = SMDR_yu_in = SMDR_yd_in = 0;
  SMDR_ymu_in = SMDR_ye_in = 0;

  SMDR_Set_Values_Interactively (nVars, varList);

  SMDR_Load_Inputs ();

  if(loopOrder > 2.99999999) {
    printf("\nPlease be patient, need to compute some 3-loop integrals...\n"); 
  }

  SMDR_Start_Timer ();

  SMDR_v = SMDR_Eval_vev (-1, loopOrder);

  SMDR_Timer ();

  printf("\n"); 
  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_m2 ();
  SMDR_Display_v ();

  printf("\nTotal calculation time: %.2f seconds\n", SMDR_Time_Total);

  return 0;
}
