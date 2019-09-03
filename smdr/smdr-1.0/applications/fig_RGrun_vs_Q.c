/* 
   Calculates the Standard Model MSbar parameters as a function of the
   renormalization scale from Q = 100 GeV to Q = 10^19 GeV, given
   their input values at a benchmark point defined in the file
   "ReferenceModel.dat" unless a different file is specified using the
   "-i" option; see below.  The full state-of-the art beta functions
   are used.

   The data is output by default in the file "FIG_RGrun_vs_Q.dat" in
   24 columns, respectively:

   1  Q
   2  Log[Q]
   3  g3 
   4  g 
   5  gp 
   6  alpha3 
   7  alpha2
   8  alpha1 
   9  1/alpha3 
   10 1/alpha2
   11 1/alpha1 
   12  yt 
   13  yb 
   14  yc 
   15  ys 
   16  yd 
   17  yu 
   18  ytau 
   19  ymu 
   20  ye 
   21  lambda 
   22  sqrt(-m2) 
   23  v
   24  Lambda

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 2.1a graphs columns 9, 10, and 11 as a function of column 2.
   Figure 2.1b graphs columns 12, 13, 14, 15, 16, 17, 18, 19, and 20 
   as a function of column 2.

   The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_RGrun_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_RGrun_vs_Q fig_RGrun_vs_Q.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.

   Run as:
   ./fig_RGrun_vs_Q
*/

#include "smdr.h"

#ifndef SMDR_PI
#define SMDR_PI     3.1415926535897932385L
#endif

/* These are the defaults: */
#define INFILENAME_DEF  "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_RGrun_vs_Q.dat"

int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  SMDR_REAL Q, log10Q;
  FILE *outfile;
  char funcname[] = "fig_RGrun_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q",
    "2  Log[Q]",
    "3  g3",
    "4  g",
    "5  gp",
    "6  alpha3",
    "7  alpha2",
    "8  alpha1",
    "9  1/alpha3",
    "10 1/alpha2",
    "11 1/alpha1",
    "12  yt",
    "13  yb",
    "14  yc",
    "15  ys",
    "16  yd",
    "17  yu",
    "18  ytau",
    "19  ymu",
    "20  ye",
    "21  lambda",
    "22  sqrt(-m2)",
    "23  v",
    "24  Lambda",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 2.1a graphs columns 9, 10, and 11 as a function of column 2, and",
    "Figure 2.1b graphs columns 12, 13, 14, 15, 16, 17, 18, 19, and 20", 
    "as a function of column 2."
  };
  int nDescriptors = 31;

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

  /* Set benchmark parameters: */
  SMDR_Read_MSbar_Inputs (inFileName);
  SMDR_Read_Value (inFileName, "SMDR_v_in");
  SMDR_Read_Value (inFileName, "SMDR_m2_in");
  SMDR_Read_Value (inFileName, "SMDR_Lambda_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n", inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();
  SMDR_Display_Lambda ();
  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_m2 (outfile, "# ");
  SMDR_Write_Lambda (outfile, "# ");
  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  printf("\nOutput data is going into \"%s\"...\n", outFileName);
  
  for (log10Q = 1.75; log10Q <= 19; log10Q += 0.25) {
    Q = SMDR_POW (10.0, log10Q);
    SMDR_RGeval_SM (Q, 5);

    fprintf (outfile, "%.6Le", Q);
    fprintf (outfile, "\t%.8Lf", log10Q);
    fprintf (outfile, "  %.8Lf", SMDR_g3);
    fprintf (outfile, "  %.8Lf", SMDR_g);
    fprintf (outfile, "  %.8Lf", SMDR_gp);
    fprintf (outfile, "  %.8Lf", SMDR_g3 * SMDR_g3/(4. * SMDR_PI));
    fprintf (outfile, "  %.8Lf", SMDR_g * SMDR_g/(4. * SMDR_PI));
    fprintf (outfile, "  %.8Lf", (5./3.) * SMDR_gp * SMDR_gp/(4. * SMDR_PI));
    fprintf (outfile, "  %.8Lf", 1./(SMDR_g3 * SMDR_g3/(4. * SMDR_PI)));
    fprintf (outfile, "  %.8Lf", 1./(SMDR_g * SMDR_g/(4. * SMDR_PI)));
    fprintf (outfile, "  %.8Lf", 1./((5./3.) * SMDR_gp * SMDR_gp/(4. * SMDR_PI)));
    fprintf (outfile, "  %.8Lf", SMDR_yt);
    fprintf (outfile, "  %.8Lf", SMDR_yb);
    fprintf (outfile, "  %.12Lf", SMDR_yc);
    fprintf (outfile, "  %.12Lf", SMDR_ys);
    fprintf (outfile, "  %.12Lf", SMDR_yd);
    fprintf (outfile, "  %.12Lf", SMDR_yu);
    fprintf (outfile, "  %.8Lf", SMDR_ytau);
    fprintf (outfile, "  %.12Lf", SMDR_ymu);
    fprintf (outfile, "  %.12Lf", SMDR_ye);
    fprintf (outfile, "  %.10Lf", SMDR_lambda);
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-SMDR_m2));
    fprintf (outfile, "  %.8Lf", SMDR_v);
    fprintf (outfile, "  %.8Lf", SMDR_Lambda);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  return 0;
}
