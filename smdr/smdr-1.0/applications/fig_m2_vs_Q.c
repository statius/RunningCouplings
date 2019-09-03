/*
   Calculates the Higgs lagrangian mass^2 parameter m2 as a function
   of the MSbar renormalization scale Q. Here m2 is determined by
   minimization of the effective potential, in various approximations
   from:
      hep-ph/0111190 (full 1-loop and 2-loop).
      1310.7553 (3-loop at leading order in g3 and yt), 
      1709.02397 (full 3-loop), and 
      1508.00912 (4-loop at leading order in g3). 

   The input parameters are obtained from the file
   "ReferenceModel.dat" unless a different file is specified using the
   "-i" option; see below.  This program produces by default an output
   file "FIG_m2_vs_Q.dat", with data in 13 columns:

   1  Q (MSbar renormalization scale)
   2  sqrt(-m2(0))   (loopOrder = 0, tree-level),
   3  sqrt(-m2(1))   (loopOrder = 1, 1-loop),
   4  sqrt(-m2(2))   (loopOrder = 2, 2-loop),
   5  sqrt(-m2(2.5)) (loopOrder = 2.5, includes leading 3-loop),
   6  sqrt(-m2(3))   (loopOrder = 3, full 3-loop),
   7  sqrt(-m2(3.5)) (loopOrder = 3, full 3-loop plus QCD 4-loop),
   8   m2(0)/m2run,
   9   m2(1)/m2run,
   10  m2(2)/m2run,
   11  m2(2.5)/m2run,
   12  m2(3)/m2run,
   13  m2(3.5)/m2run

   where m2(n) is the value of m2 obtained by requiring the effective
   potential to be minimized at the scale Q at loopOrder=n, and m2run
   is obtained by directly running the input value of m2 to the scale
   Q using 3-loop RG equations.

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 3.1a graphs columns 3, 4, 5, 6, and 7 as a function of column 1.
   Figure 3.1b graphs columns 10, 11, 12, and 13 as a function of column 1.

   The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes data to <output_filename>; if not
                          specified, the results appear in
                          "FIG_m2_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_m2_vs_Q fig_m2_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   lib3vil.a and the header file smdr.h are all in the current
   directory.

   Run as:

   ./fig_m2_vs_Q

   Running time is typically of order 15 minutes, depending on your
   hardware.
*/

#include "smdr.h"

/* These are the defaults: */
#define INFILENAME_DEF "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_m2_vs_Q.dat"

#define QSTART 60
#define QEND   271
#define QSTEP  2

int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  FILE *outfile;
  SMDR_REAL Q, result0, result1, result2, result2plus, result3, result3plus;
  char funcname[] = "fig_m2_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q (MSbar renormalization scale)",
    "2  sqrt(-m2(0))   (loopOrder = 0, tree-level)",
    "3  sqrt(-m2(1))   (loopOrder = 1, 1-loop)",
    "4  sqrt(-m2(2))   (loopOrder = 2, 2-loop)",
    "5  sqrt(-m2(2.5)) (loopOrder = 2.5, includes leading 3-loop)",
    "6  sqrt(-m2(3))   (loopOrder = 3, full 3-loop)",
    "7  sqrt(-m2(3.5)) (loopOrder = 3, full 3-loop plus QCD 4-loop)",
    "8   m2(0)/m2run",
    "9   m2(1)/m2run",
    "10  m2(2)/m2run",
    "11  m2(2.5)/m2run",
    "12  m2(3)/m2run",
    "13  m2(3.5)/m2run",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 3.1a graphs columns 3, 4, 5, 6, and 7 as a function of column 1, and",
    "Figure 3.1b graphs columns 10, 11, 12, and 13 as a function of column 1."
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
  SMDR_Read_Value (inFileName, "SMDR_m2_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n", inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_m2 (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  /* The following shouldn't matter at all, but is run by SMDR_RGeval_SM() */
  SMDR_Lambda_in = SMDR_Lambda = 0;

  
  printf("\nThis may take of order 15 minutes, depending on your hardware.\n");
  printf("Output data is going into %s...\n", outFileName);
  
  for (Q = QSTART; Q <= QEND; Q += QSTEP) {

    SMDR_RGeval_SM (Q, 5);

    result0 = SMDR_Eval_m2 (-1, 0);
    result1 = SMDR_Eval_m2 (-1, 1);
    result2 = SMDR_Eval_m2 (-1, 2);
    result2plus = SMDR_Eval_m2 (-1, 2.5);
    result3 = SMDR_Eval_m2 (-1, 3);
    result3plus = SMDR_Eval_m2 (-1, 3.5);

    fprintf (outfile, "%.1Lf", Q);
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result0));
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result1));
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result2));
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result2plus));
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result3));
    fprintf (outfile, "  %.8Lf", SMDR_SQRT(-result3plus));
    fprintf (outfile, "  %.8Lf", result0/SMDR_m2);
    fprintf (outfile, "  %.8Lf", result1/SMDR_m2);
    fprintf (outfile, "  %.8Lf", result2/SMDR_m2);
    fprintf (outfile, "  %.8Lf", result2plus/SMDR_m2);
    fprintf (outfile, "  %.8Lf", result3/SMDR_m2);
    fprintf (outfile, "  %.8Lf", result3plus/SMDR_m2);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  return 0;
}
