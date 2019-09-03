/*
   Finds the Standard Model MSbar Higgs parameters lambda and m^2 at
   their best fit to current experimental observables, and also the
   envelope of values that lambda and m^2 take as the observables are
   varied within 1 sigma and 2 sigma of their central values, as a
   function of the renormalization scale from Q = MZ to Q = 10^19 GeV.
   The current experimental values are all taken from the file
   "smdr_pdg.h".

   Output data are given by default in the file "FIG_lambda_m2_highQ.dat",
   in 12 columns:

   1  Q
   2  Log[10,Q]
   3  lambda (lower bound for observables within 2 sigma)  
   4  lambda (lower bound for observables within 1 sigma)  
   5  lambda (central value)  
   6  lambda (upper bound for observables within 1 sigma)  
   7  lambda (upper bound for observables within 2 sigma)  
   8  sqrt(-m2) (lower bound for observables within 2 sigma)  
   9  sqrt(-m2) (lower bound for observables within 1 sigma)  
   10 sqrt(-m2) (central value)  
   11 sqrt(-m2) (upper bound for observables within 1 sigma)  
   12 sqrt(-m2) (upper bound for observables within 2 sigma)  

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 4.3a graphs columns 3, 4, 5, 6, and 7 as a function of column 2.
   Figure 4.3b graphs columns 8, 9, 10, 11, and 12 as a function of column 2.

   The program takes an optional command line argument:

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_lambda_m2_highQ.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_lambda_m2_vs_Q fig_lambda_m2_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   lib3vil.a and the header file smdr.h are all in the current
   directory.

   Run as:

   ./fig_lambda_m2_vs_Q

   Running in the background is highly recommended, because the running
   time can be of order an hour, depending on your hardware.
*/

#include "smdr.h"

#ifndef SMDR_PI
#define SMDR_PI     3.1415926535897932385L
#endif

/* Default value: */
#define OUTFILENAME_DEF "FIG_lambda_m2_highQ.dat"

int main (int argc, char *argv[])
{
  char outFileName[50];
  FILE *outfile;
  SMDR_REAL ERROR_TOLERANCE = 1.0e-8;
  SMDR_REAL Q_target;
  SMDR_REAL Q, log10Q;
  int fivelog10Q;
  int iii, jjj, kkk, nnn;
  SMDR_REAL m2_hi[101][3];
  SMDR_REAL m2_lo[101][3];
  SMDR_REAL lambda_hi[101][3];
  SMDR_REAL lambda_lo[101][3];
  char funcname[] = "fig_lambda_m2_vs_Q";

  /* Define arguments: */
  int nargs = 1;
  char *arglist[] = {"-o"};
  char *argtype[] = {"string"};
  void *argvar[] = {outFileName};

  char *columnDescriptor[] = {
    "Fits from 2019 update of the Review of Particle Properties.",
    "",
    "1  Q",
    "2  Log[10,Q]",
    "3  lambda (lower bound for observables within 2 sigma)",
    "4  lambda (lower bound for observables within 1 sigma)",
    "5  lambda (central value)",
    "6  lambda (upper bound for observables within 1 sigma)",
    "7  lambda (upper bound for observables within 2 sigma)",
    "8  sqrt(-m2) (lower bound for observables within 2 sigma)",
    "9  sqrt(-m2) (lower bound for observables within 1 sigma)",
    "10 sqrt(-m2) (central value)",
    "11 sqrt(-m2) (upper bound for observables within 1 sigma)",
    "12 sqrt(-m2) (upper bound for observables within 2 sigma)",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 4.3a graphs columns 3, 4, 5, 6, and 7 as a function of column 2, and",
    "Figure 4.3b graphs columns 8, 9, 10, 11, and 12 as a function of column 2."
  };
  int nDescriptors = 20;

  /* Set default values for optional args: */
  strcpy (outFileName, OUTFILENAME_DEF);

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  Q_target = SMDR_Mt_EXPT;

  /* Initialize bins to store maximum and minimum values. */
  for (fivelog10Q = 8; fivelog10Q <= 100; fivelog10Q = fivelog10Q + 1) {
    for (nnn = 0; nnn <= 2; nnn = nnn + 1) {
      m2_hi[fivelog10Q][nnn] = -1000000; 
      m2_lo[fivelog10Q][nnn] = 1000000; 
      lambda_hi[fivelog10Q][nnn] = -100; 
      lambda_lo[fivelog10Q][nnn] = 100; 
    }
  }

  printf("\nThis may take of order 1 hour, depending on your hardware.\n");

  for (nnn = 0; nnn <= 2; nnn = nnn + 1) {
  for (iii = -nnn; iii <= nnn; iii = iii + 1) {
  for (jjj = -nnn; jjj <= nnn; jjj = jjj + 1) {
  for (kkk = -nnn; kkk <= nnn; kkk = kkk + 1) {

    printf("\nDoing %d sigma envelope: (%d %d %d)\n",nnn,iii,jjj,kkk);

    SMDR_Fit_Inputs (Q_target,
                     SMDR_alphaS_MZ_EXPT + iii * SMDR_alphaS_MZ_EXPT_UNC,
                     SMDR_alpha_EXPT,
                     SMDR_GFermi_EXPT,
                     SMDR_MZ_EXPT,
                     SMDR_Mh_EXPT + jjj * SMDR_Mh_EXPT_UNC,
                     SMDR_Mt_EXPT + kkk * SMDR_Mt_EXPT_UNC,
                     SMDR_mbmb_EXPT,
                     SMDR_mcmc_EXPT,
                     SMDR_ms_2GeV_EXPT,
                     SMDR_md_2GeV_EXPT,
                     SMDR_mu_2GeV_EXPT,
                     SMDR_Mtau_EXPT,
                     SMDR_Mmuon_EXPT,
                     SMDR_Melectron_EXPT,
                     SMDR_Delta_alpha_had_5_MZ_EXPT,
                     ERROR_TOLERANCE);

    SMDR_Display_MSbar_Parameters ();
    SMDR_Display_v ();
    SMDR_Display_m2 ();

    for (fivelog10Q = 8; fivelog10Q <= 100; fivelog10Q += 1) {
      log10Q = fivelog10Q/5.;
      Q = SMDR_POW (10.0, log10Q);
      SMDR_RGeval_SM (Q, 5);

      if (SMDR_lambda < lambda_lo[fivelog10Q][nnn]) 
          lambda_lo[fivelog10Q][nnn] = SMDR_lambda;

      if (SMDR_lambda > lambda_hi[fivelog10Q][nnn]) 
          lambda_hi[fivelog10Q][nnn] = SMDR_lambda;

      if (SMDR_m2 < m2_lo[fivelog10Q][nnn]) 
          m2_lo[fivelog10Q][nnn] = SMDR_m2;

      if (SMDR_m2 > m2_hi[fivelog10Q][nnn]) 
          m2_hi[fivelog10Q][nnn] = SMDR_m2;
    }

  }}}}

  /* Open the output file: */
  if ((outfile = fopen (outFileName, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  SMDR_Write_Version (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  for (fivelog10Q = 8; fivelog10Q <= 100; fivelog10Q += 1) {
      log10Q = fivelog10Q/5.;
      Q = SMDR_POW (10.0, log10Q);
      fprintf (outfile, "%.8Le", Q);
      fprintf (outfile, "\t%.8Lf", log10Q);

      fprintf (outfile, "  %.8Lf", lambda_lo[fivelog10Q][2]);
      fprintf (outfile, "  %.8Lf", lambda_lo[fivelog10Q][1]);
      fprintf (outfile, "  %.8Lf", lambda_hi[fivelog10Q][0]);
      fprintf (outfile, "  %.8Lf", lambda_hi[fivelog10Q][1]);
      fprintf (outfile, "  %.8Lf", lambda_hi[fivelog10Q][2]);

      fprintf (outfile, "  %.8Lf", SMDR_SQRT(-m2_lo[fivelog10Q][2]));
      fprintf (outfile, "  %.8Lf", SMDR_SQRT(-m2_lo[fivelog10Q][1]));
      fprintf (outfile, "  %.8Lf", SMDR_SQRT(-m2_hi[fivelog10Q][0]));
      fprintf (outfile, "  %.8Lf", SMDR_SQRT(-m2_hi[fivelog10Q][1]));
      fprintf (outfile, "  %.8Lf", SMDR_SQRT(-m2_hi[fivelog10Q][2]));       
      fprintf (outfile, "\n");
  } 

  fflush (outfile);
  fclose (outfile);
  return 0;
}
