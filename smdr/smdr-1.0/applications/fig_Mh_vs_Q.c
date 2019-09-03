/* 
   Calculates the Higgs pole mass as a function of the MSbar
   renormalization scale Q at which it is computed, in various
   approximations, using the results in 1407.4336. The input
   parameters are obtained from the file "ReferenceModel.dat" unless a
   different file is specified using the "-i" option; see below.

   This program produces by default an output file: "FIG_Mh_vs_Q.dat"
   with data in 7 columns, respectively
 
   1  Q  (MSbar renormalization scale)
   2  Mh (loopOrder = 0, tree-level)
   3  Mh (loopOrder = 1, 1-loop) 
   4  Mh (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   5  Mh (loopOrder = 2, full 2-loop) 
   6  Mh (loopOrder = 2.3, includes 3-loop QCD) 
   7  Mh (loopOrder = 2.5, includes 3-loop QCD and yt corrections) 

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 4.1a graphs columns 3, 4, 5, and 7 as a function of column 1.

   The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_Mh_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_Mh_vs_Q fig_Mh_vs_Q.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.

   Run as:

   ./fig_Mh_vs_Q

   Running time is typically of order 4 minutes, depending on your hardware,
   so you might want to run it in the background.
*/

#include "smdr.h"

/* These are the defaults: */
#define INFILENAME_DEF "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_Mh_vs_Q.dat"

#define QSTART 60
#define QEND   271
#define QSTEP  2

int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  FILE *outfile;
  SMDR_REAL Q;
  char funcname[] = "fig_Mh_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q  (MSbar renormalization scale)",
    "2  Mh (loopOrder = 0, tree-level)",
    "3  Mh (loopOrder = 1, 1-loop)",
    "4  Mh (loopOrder = 1.5, 1-loop plus 2-loop QCD)",
    "5  Mh (loopOrder = 2, full 2-loop)",
    "6  Mh (loopOrder = 2.3, includes 3-loop QCD)",
    "7  Mh (loopOrder = 2.5, includes 3-loop QCD and yt corrections)",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 4.1a graphs columns 3, 4, 5, and 7 as a function of column 1."
  };
  int nDescriptors = 12;

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

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n", inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  /* The following shouldn't matter at all, but are run by SMDR_RGeval_SM() */
  SMDR_m2_in = SMDR_m2 = 0;
  SMDR_Lambda_in = SMDR_Lambda = 0;

  printf("\nThis may take of order 4 minutes, depending on your hardware.\n");
  printf("Output data is going into \"%s\"...\n", outFileName);

  for (Q = QSTART; Q <= QEND; Q += QSTEP) {

    SMDR_RGeval_SM (Q, 5);
    fprintf (outfile, "%.1Lf", Q);

    SMDR_Eval_Mh_pole (-1, 0, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    SMDR_Eval_Mh_pole (-1, 1, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    SMDR_Eval_Mh_pole (-1, 1.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    SMDR_Eval_Mh_pole (-1, 2, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    SMDR_Eval_Mh_pole (-1, 2.3, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    SMDR_Eval_Mh_pole (-1, 2.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    fprintf (outfile, "  %.8Lf", SMDR_Mh_pole);

    fprintf (outfile, "\n");
    fflush (outfile);
  }
  fclose (outfile);

  return 0;
}
