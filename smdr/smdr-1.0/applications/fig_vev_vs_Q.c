/*
   Calculates the Higgs VEV as a function of the MSbar renormalization
   scale Q.  Here v is determined by minimization of the effective
   potential, in various approximations from:
      hep-ph/0111190 (full 1-loop and 2-loop).
      1310.7553 (3-loop at leading order in g3 and yt),
      1709.02397 (full 3-loop), and
      1508.00912 (4-loop at leading order in g3).

   The input parameters are obtained from the file
   "ReferenceModel.dat" unless a different file is specified using the
   "-i" option; see below. This program produces by default an output
   file "FIG_vev_vs_Q.dat", with data in 13 columns:
   
   1  Q (MSbar renormalization scale) 
   2  v(0)   (loopOrder = 0,   tree-level), 
   3  v(1)   (loopOrder = 1,   1-loop), 
   4  v(2)   (loopOrder = 2,   2-loop), 
   5  v(2.5) (loopOrder = 2.5, includes leading 3-loop), 
   6  v(3)   (loopOrder = 3,   full 3-loop), 
   7  v(3.5) (loopOrder = 3.5, full 3-loop + QCD 4-loop), 
   8  v(0)/vrun, 
   9  v(1)/vrun, 
   10  v(2)/vrun, 
   11  v(2.5)/vrun, 
   12  v(3)/vrun
   13  v(3.5)/vrun

   where v(n) is the value of v obtained by minimizing the effective
   potential at the scale Q in the given approximation, and vrun is
   obtained by directly running the input value of v to the scale Q
   using 3-loop RG equations.

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 3.2a graphs columns 3, 4, 5, 6, and 7 as a function of column 1.
   Figure 3.2b graphs columns 10, 11, 12, and 13 as a function of column 1.

  The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes data to <output_filename>; if not
                          specified, the results appear in
                          "FIG_vev_vs_Q.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_vev_vs_Q fig_vev_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   lib3vil.a and the header file smdr.h are all in the current
   directory.

   Run as:

   ./fig_vev_vs_Q

   Running in the background is highly recommended, as the running time is
   typically of order 2 hours, depending on your hardware.
*/

#include "smdr.h"

/* These are the defaults: */
#define INFILENAME_DEF "ReferenceModel.dat"
#define OUTFILENAME_DEF "FIG_vev_vs_Q.dat"

#define QSTART 60
#define QEND   241
#define QSTEP  2

int main (int argc, char *argv[])
{
  char inFileName[50], outFileName[50];
  FILE *outfile;
  SMDR_REAL Q, resultvtree, resultv1, resultv2, resultv2plus; 
  SMDR_REAL resultv3, resultv3plus;
  char funcname[] = "fig_vev_vs_Q";

  /* Define arguments: */
  int nargs = 2;
  char *arglist[] = {"-i","-o"};
  char *argtype[] = {"string","string"};
  void *argvar[] = {inFileName, outFileName};

  char *columnDescriptor[] = {
    "1  Q (MSbar renormalization scale)",
    "2  v(0)   (loopOrder = 0,   tree-level)",
    "3  v(1)   (loopOrder = 1,   1-loop)",
    "4  v(2)   (loopOrder = 2,   2-loop)",
    "5  v(2.5) (loopOrder = 2.5, includes leading 3-loop)",
    "6  v(3)   (loopOrder = 3,   full 3-loop)",
    "7  v(3.5) (loopOrder = 3.5, full 3-loop + QCD 4-loop)",
    "8  v(0)/vrun",
    "9  v(1)/vrun",
    "10  v(2)/vrun",
    "11  v(2.5)/vrun",
    "12  v(3)/vrun",
    "13  v(3.5)/vrun",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 3.2a graphs columns 3, 4, 5, 6, and 7 as a function of column 1, and",
    "Figure 3.2b graphs columns 10, 11, 12, and 13 as a function of column 1."
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

  /* Set benchmark parameters: */
  SMDR_Read_MSbar_Inputs (inFileName);
  SMDR_Read_Value (inFileName, "SMDR_v_in");
  SMDR_Read_Value (inFileName, "SMDR_m2_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n",inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();

  SMDR_Write_MSbar_Parameters (outfile, "# ");
  SMDR_Write_v (outfile, "# ");
  SMDR_Write_m2 (outfile, "# ");

  SMDR_Write_Column_Data (outfile, nDescriptors, columnDescriptor, "# ");

  /* The following shouldn't matter at all, but is run by SMDR_RGeval_SM() */
  SMDR_Lambda_in = SMDR_Lambda = 0;

  printf("\nThis may take of order 2 hours, depending on your hardware.");
  printf("\nOutput data is going into %s...\n", outFileName);
  
  for (Q = QSTART; Q <= QEND; Q += QSTEP) {

    SMDR_RGeval_SM (Q, 5);

    resultvtree = SMDR_Eval_vev (-1, 0);
    resultv1 = SMDR_Eval_vev (-1, 1);
    resultv2 = SMDR_Eval_vev (-1, 2);
    resultv2plus = SMDR_Eval_vev (-1, 2.5);
    resultv3 = SMDR_Eval_vev (-1, 3);
    resultv3plus = SMDR_Eval_vev (-1, 3.5);

    fprintf (outfile, "%.1Lf", Q);
    fprintf (outfile, "  %.8Lf", resultvtree);
    fprintf (outfile, "  %.8Lf", resultv1);
    fprintf (outfile, "  %.8Lf", resultv2);
    fprintf (outfile, "  %.8Lf", resultv2plus);
    fprintf (outfile, "  %.8Lf", resultv3);
    fprintf (outfile, "  %.8Lf", resultv3plus);
    fprintf (outfile, "  %.8Lf", resultvtree/SMDR_v);
    fprintf (outfile, "  %.8Lf", resultv1/SMDR_v);
    fprintf (outfile, "  %.8Lf", resultv2/SMDR_v);
    fprintf (outfile, "  %.8Lf", resultv2plus/SMDR_v);
    fprintf (outfile, "  %.8Lf", resultv3/SMDR_v);
    fprintf (outfile, "  %.8Lf", resultv3plus/SMDR_v);
    fprintf (outfile, "\n");
    fflush (outfile);
  }

  fclose (outfile);

  return 0;
}
