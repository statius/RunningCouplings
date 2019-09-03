/* 
   Calculates the top-quark pole mass as a function of the MSbar
   renormalization scale Q at which it is computed, in various
   approximations, using the results of 1604.01134 and references
   therein. The input parameters are obtained from the file
   "ReferenceModel.dat" unless a different file is specified using the
   "-i" option; see below. This program produces output data in 12
   columns:

   1    Q 
   2    Mt tree
   3    Mt 1-loop QCD
   4    Mt 2-loop QCD
   5    Mt 3-loop QCD
   6    Mt 4-loop QCD
   7    Mt     4-loop QCD + 1-loop nonQCD
   8    Gammat 4-loop QCD + 1-loop nonQCD
   9    Mt     4-loop QCD + 1-loop nonQCD + 2-loop mixed
   10   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop mixed
   11   Mt     4-loop QCD + 1-loop nonQCD + 2-loop full  
   12   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop full

   Two output files are produced: "FIG_Mt_vs_Q_0.dat" and
   "FIG_Mt_vs_Q_1.dat", using expansions around the tree-level and
   real pole masses, respectively.

   In the paper: arXiv:1907.02500
   Standard Model parameters in the tadpole-free pure MSbar scheme
   by Stephen P. Martin and David G. Robertson,
   Figure 4.1b graphs columns 6, 7, 9, and 11 as a
   function of column 1 from FIG_Mt_vs_Q_1.dat.

   The program takes an optional command line argument:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

   gcc -o fig_Mt_vs_Q fig_Mt_vs_Q.c -L. -lsmdr -ltsil -lm

   assuming the necessary archive files libsmdr.a and libtsil.a and
   the header file smdr.h are all in the current directory.

   Run as:

   ./fig_Mt_vs_Q
   
   The running time is typically of order 6 minutes, depending on your
   hardware, so you might want to run it in the background.
*/

#include "smdr.h"

#define INFILENAME_DEF "ReferenceModel.dat"

#define QSTART 50
#define QEND   301
#define QSTEP  2

int main (int argc, char *argv[])
{
  char inFileName[50];
  FILE *outfile;
  char outFileName0[] = "FIG_Mt_vs_Q_0.dat";
  char outFileName1[] = "FIG_Mt_vs_Q_1.dat";
  char funcname[] = "fig_Mt_vs_Q";

  SMDR_REAL Q;
  /* SMDR_COMPLEX CM2T; */
  /* SMDR_REAL MTreal, GammaT; */
  int method;

  /* Define arguments: */
  int nargs = 1;
  char *arglist[] = {"-i"};
  char *argtype[] = {"string"};
  void *argvar[] = {inFileName};

  char *columnDescriptor0[] = {
    "1    Q",
    "2    Mt tree",
    "3    Mt 1-loop QCD",
    "4    Mt 2-loop QCD",
    "5    Mt 3-loop QCD",
    "6    Mt 4-loop QCD",
    "7    Mt     4-loop QCD + 1-loop nonQCD",
    "8    Gammat 4-loop QCD + 1-loop nonQCD",
    "9    Mt     4-loop QCD + 1-loop nonQCD + 2-loop mixed",
    "10   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop mixed",
    "11   Mt     4-loop QCD + 1-loop nonQCD + 2-loop full",
    "12   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop full",
  };
  int nDescriptors0 = 12;

  char *columnDescriptor1[] = {
    "1    Q",
    "2    Mt tree",
    "3    Mt 1-loop QCD",
    "4    Mt 2-loop QCD",
    "5    Mt 3-loop QCD",
    "6    Mt 4-loop QCD",
    "7    Mt     4-loop QCD + 1-loop nonQCD",
    "8    Gammat 4-loop QCD + 1-loop nonQCD",
    "9    Mt     4-loop QCD + 1-loop nonQCD + 2-loop mixed",
    "10   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop mixed",
    "11   Mt     4-loop QCD + 1-loop nonQCD + 2-loop full",
    "12   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop full",
    "",
    "In the paper: arXiv:1907.02500",
    "Standard Model parameters in the tadpole-free pure MSbar scheme",
    "by Stephen P. Martin and David G. Robertson,",
    "Figure 4.1b graphs columns 6, 7, 9, and 11 as a function of column 1." 
  };
  int nDescriptors1 = 17;

  /* Set default values for optional args: */
  strcpy (inFileName, INFILENAME_DEF);

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  /* Set benchmark parameters */
  SMDR_Read_MSbar_Inputs (inFileName);
  SMDR_Read_Value (inFileName, "SMDR_v_in");

  SMDR_Load_Inputs ();

  printf("\nInput MSbar parameters read from \"%s\":\n\n", inFileName);

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  
  /* The following shouldn't matter at all, but are run by SMDR_RGeval_SM() */
  SMDR_m2_in = SMDR_m2 = 0;
  SMDR_Lambda_in = SMDR_Lambda = 0;

  printf("\nThis may take of order 6 minutes, depending on your hardware.\n");
  printf("Output data is going into %s and %s...\n", outFileName0, outFileName1);

  for (method = 0; method < 2; method++) {

    /* Open output data file: */
    if (0 == method) 
      if ((outfile = fopen (outFileName0, "w")) == NULL)
	SMDR_Error (funcname, "Output file cannot be opened.", 246);
    if (1 == method)
      if ((outfile = fopen (outFileName1, "w")) == NULL)
	SMDR_Error (funcname, "Output file cannot be opened.", 246);

    SMDR_Load_Inputs ();

    SMDR_Write_Version (outfile, "# ");
    SMDR_Write_MSbar_Parameters (outfile, "# ");
    SMDR_Write_v (outfile, "# ");
    if (0 == method) 
      SMDR_Write_Column_Data (outfile, nDescriptors0, columnDescriptor0, "# ");
    else
      SMDR_Write_Column_Data (outfile, nDescriptors1, columnDescriptor1, "# ");

    for (Q = QSTART; Q <= QEND; Q += QSTEP) {

      SMDR_RGeval_SM (Q, 5);

      fprintf (outfile, "%.1Lf", Q); 

      SMDR_Eval_Mt_pole (-1, method, 0, 0, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);

      SMDR_Eval_Mt_pole (-1, method, 1, 0, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);

      SMDR_Eval_Mt_pole (-1, method, 2, 0, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);

      SMDR_Eval_Mt_pole (-1, method, 3, 0, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);

      SMDR_Eval_Mt_pole (-1, method, 4, 0, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);

      SMDR_Eval_Mt_pole (-1, method, 4, 1, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Gammat_pole);

      SMDR_Eval_Mt_pole (-1, method, 4, 1.5, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Gammat_pole);
    
      SMDR_Eval_Mt_pole (-1, method, 4, 2, &SMDR_Mt_pole, &SMDR_Gammat_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Mt_pole);
      fprintf (outfile, "\t%.8Lf", SMDR_Gammat_pole);

      fprintf (outfile, "\n");
      fflush (outfile);
    } /* End of Q loop */

    fclose (outfile);
  } /* End of method loop */

  return 0;
}
