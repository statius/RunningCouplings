/* 
   Calculates the Higgs self coupling lambda as a function of the
   MSbar renormalization scale Q. The inputs are the experimental
   Higgs pole mass SMDR_Mh_EXPT from smdr_pdg.h and the other MSbar
   input parameters obtained from the file "ReferenceModel.dat" unless
   a different file is specified using the "-i" option; see
   below. This program produces (by default) an output file
   "FIG_lambda_vs_Q.dat" containing data in 13 columns:

   1  Q (MSbar renormalization scale)
   2  lambda (tree-level)
   3  lambda (1-loop)
   4  lambda (1-loop plus 2-loop QCD)
   5  lambda (2-loop)
   6  lambda (2-loop plus leading 3-loop QCD)
   7  lambda (2-loop plus leading 3-loop QCD and yt)
   8  lambda/lambdarun (tree-level)
   9  lambda/lambdarun (1-loop)
   10 lambda/lambdarun (1-loop plus 2-loop QCD)
   11 lambda/lambdarun (2-loop)
   12 lambda/lambdarun (2-loop plus leading 3-loop QCD)
   13 lambda/lambdarun (2-loop plus leading 3-loop QCD and yt)
 
   where lambda is the value obtained at Q from the input Higgs pole
   mass with all other pertinent MSbar parameters having been run from
   their benchmark values, and lambdarun is instead obtained by
   directly RG running its benchmark input value to the scale Q.

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme 
   by Stephen P. Martin and David G. Robertson,
   Figure 4.2a graphs columns 3, 4, 5, and 7 as a function of column 1.
   Figure 4.2b graphs columns 9, 10, 11, and 13 as a function of column 1.

   The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_lambda_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_lambda_vs_Q fig_lambda_vs_Q.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.
 
   Run as:

   ./fig_lambda_vs_Q

   The running time is typically of order 3 minutes, depending on your
   hardware.
*/

#include "smdr.h"

#define QSTART 50
#define QEND   271
#define QSTEP  2

#define Mhpoletarget SMDR_Mh_EXPT

/* These are the defaults: */
#define INFILENAME_DEF "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_lambda_vs_Q.dat"


int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  FILE *outfile;
  SMDR_REAL Q;
  SMDR_REAL lambda_result0, lambda_result1, lambda_result15;
  SMDR_REAL lambda_result2, lambda_result23, lambda_result25;
  char funcname[] = "fig_lambda_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q (MSbar renormalization scale)",
    "2  lambda (tree-level)",
    "3  lambda (1-loop)",
    "4  lambda (1-loop plus 2-loop QCD)",
    "5  lambda (2-loop)",
    "6  lambda (2-loop plus leading 3-loop QCD)",
    "7  lambda (2-loop plus leading 3-loop QCD and yt)",
    "8  lambda/lambdarun (tree-level)",
    "9  lambda/lambdarun (1-loop)",
    "10 lambda/lambdarun (1-loop plus 2-loop QCD)",
    "11 lambda/lambdarun (2-loop)",
    "12 lambda/lambdarun (2-loop plus leading 3-loop QCD)",
    "13 lambda/lambdarun (2-loop plus leading 3-loop QCD and yt)",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme", 
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 4.2a graphs columns 3, 4, 5, and 7 as a function of column 1, and",
    "Figure 4.2b graphs columns 9, 10, 11, and 13 as a function of column 1."
  };
  int nDescriptors = 19;

  /* Set default values for optional args: */
  strcpy (inFileName, INFILENAME_DEF);
  strcpy (outFileName, OUTFILENAME_DEF);

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Open the output file: */
  if ((outfile = fopen (outFileName, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  /* Print version information: */
  SMDR_Display_Version ();
  SMDR_Write_Version (outfile, "# ");

  /* Set benchmark parameters */
  SMDR_Read_MSbar_Inputs (inFileName);
  SMDR_Read_Value (inFileName, "SMDR_v_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n", inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  printf ("Mh = %.8Lf (pole mass)\n", Mhpoletarget);

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  fprintf (outfile, "# Mh = %.8Lf (pole mass)\n", Mhpoletarget);

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  /* The following shouldn't matter at all, but are run by SMDR_RGeval_SM() */
  SMDR_m2_in = SMDR_m2 = 0;
  SMDR_Lambda_in = SMDR_Lambda = 0;

  printf("\nThis may take of order 3 minutes, depending on your hardware.\n");
  printf("Output data is going into %s...\n", outFileName);

  for (Q = QSTART; Q <= QEND; Q += QSTEP) {

    SMDR_RGeval_SM (Q, 5);
    lambda_result0 = SMDR_Eval_lambda (-1, Mhpoletarget, 0);
    lambda_result1 = SMDR_Eval_lambda (-1, Mhpoletarget, 1);
    lambda_result15 = SMDR_Eval_lambda (-1, Mhpoletarget, 1.5);
    lambda_result2 = SMDR_Eval_lambda (-1, Mhpoletarget, 2);
    lambda_result23 = SMDR_Eval_lambda (-1, Mhpoletarget, 2.3);
    lambda_result25 = SMDR_Eval_lambda (-1, Mhpoletarget, 2.5);

    fprintf (outfile, "%.1Lf", Q);
    fprintf (outfile, "  %.8Lf", lambda_result0);
    fprintf (outfile, "  %.8Lf", lambda_result1);
    fprintf (outfile, "  %.8Lf", lambda_result15);
    fprintf (outfile, "  %.8Lf", lambda_result2);
    fprintf (outfile, "  %.8Lf", lambda_result23);
    fprintf (outfile, "  %.8Lf", lambda_result25);
    fprintf (outfile, "  %.8Lf", lambda_result0/SMDR_lambda);
    fprintf (outfile, "  %.8Lf", lambda_result1/SMDR_lambda);
    fprintf (outfile, "  %.8Lf", lambda_result15/SMDR_lambda);
    fprintf (outfile, "  %.8Lf", lambda_result2/SMDR_lambda);
    fprintf (outfile, "  %.8Lf", lambda_result23/SMDR_lambda);
    fprintf (outfile, "  %.8Lf", lambda_result25/SMDR_lambda);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  return 0;
}
