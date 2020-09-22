/* 
   Prints the running MSbar parameters in the QCD+QED effective theory
   as a function of the renormalization scale from Q = MZ down to Q =
   1 GeV, decoupling the bottom, tau, and charm quarks at appropriate
   thresholds.  The input parameters come from the file
   "ReferenceModel.dat" unless a different file is specified using the
   "-i" option; see below. The state-of-the-art beta functions and
   decoupling relations are used.

   Data are output in the following files: 

   "FIG_RGrun_53.dat" for the 5-quark, 3-lepton theory, in 11 columns, 
     1  Q (from SMDR_MZ_EXPT down to SMDR_mbmb_EXPT)
     2  alphaS 
     3  alpha 
     4  mb 
     5  mc 
     6  ms 
     7  md
     8  mu
     9  mtau 
     10  mmuon 
     11  melectron 

   "FIG_RGrun_43.dat" for the 4-quark, 3-lepton theory, in 10 columns,
     1  Q (from SMDR_mbmb_EXPT down to SMDR_Mtau_EXPT)
     2  alphaS 
     3  alpha 
     4  mc 
     5  ms 
     6  md
     7  mu
     8  mtau 
     9  mmuon 
     10  melectron 

   "FIG_RGrun_42.dat" for the 4-quark, 2-lepton theory, in 9 columns,
     1  Q (from SMDR_Mtau_EXPT down to SMDR_mcmc_EXPT)
     2  alphaS 
     3  alpha 
     4  mc 
     5  ms 
     6  md
     7  mu
     8  mmuon 
     9  melectron 

   "FIG_RGrun_32.dat" for the 3-quark, 2-lepton theory, in 8 columns,
     1  Q (from SMDR_mcmc_EXPT down to 1 GeV)
     2  alphaS 
     3  alpha 
     4  ms 
     5  md
     6  mu
     7  mmuon 
     8  melectron 

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 2.2a graphs columns 2 and 3 as a function of column 1 from
     FIG_RGrun_53.dat FIG_RGrun_43.dat FIG_RGrun_42.dat FIG_RGrun_32.dat
   Figure 2.2b graphs 
     columns 4, 5, 6, 7, 8, 9, 10, and 11 as a function of column 1 
       from FIG_RGrun_53.dat
     columns 4, 5, 6, 7, 8, 9, and 10 as a function of column 1 
       from FIG_RGrun_43.dat
     columns 4, 5, 6, 7, 8, and 9 as a function of column 1 
       from FIG_RGrun_42.dat
     columns 4, 5, 6, 7, and 8 as a function of column 1 
       from FIG_RGrun_32.dat

   The program takes an optional command line argument:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_RGrun_QCDQED fig_RGrun_QCDQED.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.

   Run as: 
   ./fig_RGrun_QCDQED
*/

#include "smdr.h"

#define INFILENAME_DEF "ReferenceModel.dat"

int main (int argc, char *argv[])
{
  char inFileName[50];
  int j;
  SMDR_REAL Q, deltalogQ;
  FILE *outfile;
  char outFileName53[] = "FIG_RGrun_53.dat";
  char outFileName43[] = "FIG_RGrun_43.dat";
  char outFileName42[] = "FIG_RGrun_42.dat";
  char outFileName32[] = "FIG_RGrun_32.dat";
  char funcname[] = "fig_RGrun_QCDQED";

  /* Define arguments: */
  int nargs = 1;
  char *arglist[] = {"-i"};
  char *argtype[] = {"string"};
  void *argvar[] = {inFileName};

  char *columnDescriptor53[] = {
    "1  Q (from SMDR_MZ_EXPT down to SMDR_mbmb_EXPT)",
    "2  alphaS",
    "3  alpha",
    "4  mb",
    "5  mc",
    "6  ms",
    "7  md",
    "8  mu",
    "9  mtau",
    "10  mmuon",
    "11  melectron",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 2.2a graphs columns 2 and 3 as a function of column 1, and",
    "Figure 2.2b graphs columns 4, 5, 6, 7, 8, 9, 10, and 11 as a", 
    "function of column 1."
  };
  int nDescriptors53 = 18;

  char *columnDescriptor43[] = {
    "1  Q (from SMDR_mbmb_EXPT down to SMDR_Mtau_EXPT)",
    "2  alphaS",
    "3  alpha",
    "4  mc",
    "5  ms",
    "6  md",
    "7  mu",
    "8  mtau",
    "9  mmuon",
    "10  melectron",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 2.2a graphs columns 2 and 3 as a function of column 1, and",
    "Figure 2.2b graphs columns 4, 5, 6, 7, 8, 9, and 10 as a function", 
    "of column 1."
  };
  int nDescriptors43 = 17;

  char *columnDescriptor42[] = {
    "1  Q (from SMDR_Mtau_EXPT down to SMDR_mcmc_EXPT)",
    "2  alphaS",
    "3  alpha",
    "4  mc",
    "5  ms",
    "6  md",
    "7  mu",
    "8  mmuon",
    "9  melectron",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 2.2a graphs columns 2 and 3 as a function of column 1, and",
    "Figure 2.2b graphs columns 4, 5, 6, 7, 8, and 9 as a function of column 1."
  };
  int nDescriptors42 = 15;

  char *columnDescriptor32[] = {
    "1  Q (from SMDR_mcmc_EXPT down to 1 GeV)",
    "2  alphaS",
    "3  alpha",
    "4  ms",
    "5  md",
    "6  mu",
    "7  mmuon",
    "8  melectron",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 2.2a graphs columns 2 and 3 as a function of column 1, and",
    "Figure 2.2b graphs columns 4, 5, 6, 7, and 8 as a function of column 1."
  };
  int nDescriptors32 = 14;

  /* Set default values for optional args: */
  strcpy (inFileName, INFILENAME_DEF);

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  /* Set benchmark parameters: */
  SMDR_Read_MSbar_Inputs (inFileName);
  SMDR_Read_Value (inFileName, "SMDR_v_in");
  SMDR_Read_Value (inFileName, "SMDR_Delta_alpha_had_5_MZ_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n",inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_Delta_alpha_had5 ();

  /* The following shouldn't matter at all, but are run by SMDR_RGeval_SM() */
  SMDR_m2_in = SMDR_m2 = 0;
  SMDR_Lambda_in = SMDR_Lambda = 0;

  /* Open output data file for 5-quark, 3-lepton theory: */
  if ((outfile = fopen (outFileName53, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  SMDR_Write_Version (outfile, "# ");

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_Delta_alpha_had5 (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors53, columnDescriptor53, "# ");

  deltalogQ = 0.05 * SMDR_LOG(SMDR_mbmb_EXPT/SMDR_MZ_EXPT);

  for (j = 0; j <= 20; j++) {
    Q = SMDR_MZ_EXPT * SMDR_EXP (j * deltalogQ);
    SMDR_RGeval_QCDQED_53 (Q, SMDR_MZ_EXPT, 4);

    fprintf (outfile, "%Lf", Q);
    fprintf (outfile, "  %.8Lf", SMDR_alphaS_53);
    fprintf (outfile, "  %.8Lf", SMDR_alpha_53);
    fprintf (outfile, "  %.8Lf", SMDR_mb_53);
    fprintf (outfile, "  %.8Lf", SMDR_mc_53);
    fprintf (outfile, "  %.8Lf", SMDR_ms_53);
    fprintf (outfile, "  %.8Lf", SMDR_md_53);
    fprintf (outfile, "  %.8Lf", SMDR_mu_53);
    fprintf (outfile, "  %.8Lf", SMDR_mtau_53);
    fprintf (outfile, "  %.8Lf", SMDR_mmuon_53);
    fprintf (outfile, "  %.8Lf", SMDR_melectron_53);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  /* Open output data file for 4-quark, 3-lepton theory: */
  if ((outfile = fopen (outFileName43, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  SMDR_Load_Inputs ();

  SMDR_Write_Version (outfile, "# ");
  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_Delta_alpha_had5 (outfile, "# ");
  SMDR_Write_Column_Data (outfile, nDescriptors43, columnDescriptor43, "# ");

  deltalogQ = 0.05 * SMDR_LOG(SMDR_Mtau_EXPT/SMDR_mbmb_EXPT);

  for (j = 0; j <= 20; j++) {
    Q = SMDR_mbmb_EXPT * SMDR_EXP (j * deltalogQ);
    SMDR_RGeval_QCDQED_43 (Q, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 4);

    fprintf (outfile, "%Lf", Q);
    fprintf (outfile, "  %.8Lf", SMDR_alphaS_43);
    fprintf (outfile, "  %.8Lf", SMDR_alpha_43);
    fprintf (outfile, "  %.8Lf", SMDR_mc_43);
    fprintf (outfile, "  %.8Lf", SMDR_ms_43);
    fprintf (outfile, "  %.8Lf", SMDR_md_43);
    fprintf (outfile, "  %.8Lf", SMDR_mu_43);
    fprintf (outfile, "  %.8Lf", SMDR_mtau_43);
    fprintf (outfile, "  %.8Lf", SMDR_mmuon_43);
    fprintf (outfile, "  %.8Lf", SMDR_melectron_43);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  /* Open output data file for 4-quark, 2-lepton theory: */
  if ((outfile = fopen (outFileName42, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  deltalogQ = 0.05 * SMDR_LOG(SMDR_mcmc_EXPT/SMDR_Mtau_EXPT);

  SMDR_Load_Inputs ();

  SMDR_Write_Version (outfile, "# ");
  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_Delta_alpha_had5 (outfile, "# ");
  SMDR_Write_Column_Data (outfile, nDescriptors42, columnDescriptor42, "# ");

  for (j = 0; j <= 20; j++) {
    Q = SMDR_Mtau_EXPT * SMDR_EXP (j * deltalogQ);
    SMDR_RGeval_QCDQED_42 (Q, SMDR_Mtau_EXPT,SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 4);

    fprintf (outfile, "%Lf", Q);
    fprintf (outfile, "  %.8Lf", SMDR_alphaS_42);
    fprintf (outfile, "  %.8Lf", SMDR_alpha_42);
    fprintf (outfile, "  %.8Lf", SMDR_mc_42);
    fprintf (outfile, "  %.8Lf", SMDR_ms_42);
    fprintf (outfile, "  %.8Lf", SMDR_md_42);
    fprintf (outfile, "  %.8Lf", SMDR_mu_42);
    fprintf (outfile, "  %.8Lf", SMDR_mmuon_42);
    fprintf (outfile, "  %.8Lf", SMDR_melectron_42);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  /* Open output data file for 3-quark, 2-lepton theory: */
  if ((outfile = fopen (outFileName32, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  deltalogQ = 0.005 * SMDR_LOG(0.6/SMDR_mcmc_EXPT);

  SMDR_Load_Inputs ();

  SMDR_Write_Version (outfile, "# ");
  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_Delta_alpha_had5 (outfile, "# ");
  SMDR_Write_Column_Data (outfile, nDescriptors32, columnDescriptor32, "# ");

  for (j = 0; j <= 200; j++) {
    Q = SMDR_mcmc_EXPT * SMDR_EXP (j * deltalogQ);
    SMDR_RGeval_QCDQED_32 (Q, SMDR_mcmc_EXPT, SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, 
                           SMDR_MZ_EXPT, 4);

    fprintf (outfile, "%Lf", Q);
    fprintf (outfile, "  %.8Lf", SMDR_alphaS_32);
    fprintf (outfile, "  %.8Lf", SMDR_alpha_32);
    fprintf (outfile, "  %.8Lf", SMDR_ms_32);
    fprintf (outfile, "  %.8Lf", SMDR_md_32);
    fprintf (outfile, "  %.8Lf", SMDR_mu_32);
    fprintf (outfile, "  %.8Lf", SMDR_mmuon_32);
    fprintf (outfile, "  %.8Lf", SMDR_melectron_32);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  return 0;
}
