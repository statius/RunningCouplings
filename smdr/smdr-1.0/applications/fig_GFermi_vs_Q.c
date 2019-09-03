/* 
   Calculates GFermi as a function of the MSbar renormalization scale
   Q at which it is computed, in various approximations. The input
   parameters are obtained from the file "ReferenceModel.dat". This
   program produces by default an output file "FIG_GFermi_vs_Q.dat",
   with data in 11 columns, respectively

   1  Q  (MSbar renormalization scale) 
   2  GFermi 10^5 (loopOrder = 0, tree-level) 
   3  GFermi 10^5 (loopOrder = 1, 1-loop) 
   4  GFermi 10^5 (loopOrder = 1.5, 1-loop plus 2-loop LO in QCD) 
   5  GFermi 10^5 (loopOrder = 2, full 2-loop) 
   6, 7, 8, 9 = same, but divided by GFermi(experiment)

   In the paper: arXiv:1907.02500 
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson, 
   Figure 3.3 graphs columns 3, 4, and 5 as a function of column 1.

   This program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_GFermi_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_GFermi_vs_Q fig_GFermi_vs_Q.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.

   Run as:
   ./fig_GFermi_vs_Q 
*/

#include "smdr.h"

/* These are the defaults: */
#define INFILENAME_DEF "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_GFermi_vs_Q.dat"

#define QSTART 50
#define QEND   271
#define QSTEP  1

int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  FILE *outfile;
  SMDR_REAL Q, result0, result1, result15, result2;
  char funcname[] = "fig_GFermi_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q           (MSbar renormalization scale)",
    "2  GFermi 10^5 (loopOrder = 0, tree-level)",
    "3  GFermi 10^5 (loopOrder = 1, 1-loop)",
    "4  GFermi 10^5 (loopOrder = 1.5, 1-loop plus 2-loop LO in QCD)",
    "5  GFermi 10^5 (loopOrder = 2, full 2-loop)",
    "6  GFermi/GFermi(experiment) (loopOrder = 0, tree-level)", 
    "7  GFermi/GFermi(experiment) (loopOrder = 1, 1-loop)", 
    "8  GFermi/GFermi(experiment) (loopOrder = 1.5, 1-loop plus 2-loop LO in QCD)", 
    "9  GFermi/GFermi(experiment) (loopOrder = 2, full 2-loop)", 
    "",
    "In the paper: arXiv:1907.02500", 
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,", 
    "Figure 3.3 graphs columns 3, 4, and 5 as a function of column 1."
  };
  int nDescriptors = 14;

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

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  /* The following shouldn't matter at all, but are run by SMDR_RGeval_SM() */
  SMDR_m2_in = SMDR_m2 = 0;
  SMDR_Lambda_in = SMDR_Lambda = 0;

  printf("\nOutput data is going into %s...\n", outFileName);

  for (Q = QSTART; Q <= QEND; Q += QSTEP) {

    SMDR_RGeval_SM (Q, 5);

    result0 = SMDR_Eval_GFermi (-1, 0);
    result1 = SMDR_Eval_GFermi (-1, 1);
    result15 = SMDR_Eval_GFermi (-1, 1.5);
    result2 = SMDR_Eval_GFermi (-1, 2);

    fprintf (outfile, "%.1Lf", Q);  
    fprintf (outfile, "  %.9Lf", 100000 * result0 );
    fprintf (outfile, "  %.9Lf", 100000 * result1 );
    fprintf (outfile, "  %.9Lf", 100000 * result15 );
    fprintf (outfile, "  %.9Lf", 100000 * result2 );
    fprintf (outfile, "  %.9Lf", result0/SMDR_GFermi_EXPT);
    fprintf (outfile, "  %.9Lf", result1/SMDR_GFermi_EXPT);
    fprintf (outfile, "  %.9Lf", result15/SMDR_GFermi_EXPT);
    fprintf (outfile, "  %.9Lf", result2/SMDR_GFermi_EXPT);
    
    fprintf (outfile, "\n");
    fflush (outfile);
  }
  fclose (outfile);

  return 0;
}
