#include "smdr_internal.h"
#include "smdr_ioparams.h"

/* General I/O Routines */

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Skips everything until the next EOL; leaves the file stream
   pointing at the first character of the following line. */

int SkipRestOfLine (FILE *fp)
{
  char tmp[160];
  fgets (tmp, 160, fp);

  return 0;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Skips spaces and tabs and leaves the file stream pointing at the
   first character that is not a space or tab (this could be EOL,
   i.e., it does *not* proceed to the next line). */

int SkipWhiteSpace (FILE *fp)
{
  char a;
  while (isspace(a=getc(fp)) && a != '\n')
    ;
  fseek(fp, -1, SEEK_CUR);

  return 0;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* This can be called even if we are not at the beginning of the word;
   any leading white space will be skipped.

   Returns 0 if a word was read with either a terminating EOL or ';',
   -1 if a word was read with a terminating '/', and 1 if only white
   space was found up to the next EOL.

   If EOL was encountered, file pointer is at the EOL; if a word was
   read, the pointer is at the first character after the word.
*/

int GetWord (FILE *fp, char *wd)
{
  char a;
  int i = 0;
  int retval = 0;

  SkipWhiteSpace (fp);

  if ((a=getc(fp)) == '\n') {
    fseek(fp, -1, SEEK_CUR);
    return 1;
  }
  else
    wd[i++] = a;

  while ( !isspace(a=getc(fp)) && a != ';' && a != '/')
    wd[i++] = a;

  wd[i] = '\0';

  if (a == ';')
    ;
  else if (a == '/')
    retval = -1;
  else
    fseek(fp, -1, SEEK_CUR);

  return retval;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int SMDR_Get_Value (FILE *fp, char *whichVar)
{
  SMDR_REAL value;
  int nVals = 35;
  char a, tmpword[30];
  int i, check;
  int success = 0;
  char funcname[] = "SMDR_Get_Value";
  char warnmsg[100];

  /* Which one are we looking for? */
  for (i=0; i<nVals; i++)
    if (strcmp(whichVar, varName[i]) == 0)
      break;

  /* Start at beginning of file always: */
  rewind (fp);

  /* Skip comments and look for the next line of data: */
  while ((a = getc(fp)) != EOF) {

    if (a == '#')
      SkipRestOfLine (fp);

    else if (a=='\n')
      ;

    else if (isspace(a)) {  /* Possible blank line or comment,
			       otherwise a format error. */
      SkipWhiteSpace (fp);
      if ((a=getc(fp)) == '\n')
	;
      else if (a == '#')
	SkipRestOfLine (fp);
    }
    else {
      fseek(fp, -1, SEEK_CUR);
      GetWord (fp, tmpword);
      
      if (strcmp(tmpword, whichVar) == 0) {
	check = GetWord (fp, tmpword);
	if (strcmp(tmpword, "=") == 0)
	  check = GetWord (fp, tmpword);
	value = (SMDR_REAL) strtold (tmpword, NULL);
	/* value = (SMDR_REAL) atof (tmpword); */
	if (check == -1) {
	  check = GetWord (fp, tmpword);
	  value /= (SMDR_REAL) strtold (tmpword, NULL);
	}
	success = 1;
	break;
      }
      else
	SkipRestOfLine (fp);
    }
  }

  /* Set value: */
  if (success == 1)
    *varValue[i] = value;
  else {
    strcpy (warnmsg, "Variable not found: ");
    strcat (warnmsg, whichVar);
    SMDR_Warn (funcname, warnmsg);
  }

  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Reads a single specified parameter from a file.                     */

int SMDR_Read_Value (char *infile, char *whichVar)
{
  FILE *fp;
  int success;
  char funcname[] = "SMDR_Read_Value";

  if ((fp = fopen(infile, "r")) == NULL)
    SMDR_Error (funcname, "Specified input file does not exist.", 123);

  success = SMDR_Get_Value (fp, whichVar);

  fclose(fp);

  if (success == 0)
    SMDR_Error (funcname, "Missing input datum.", 111);

  return (success);
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Test routine: takes a list of values and reads them.  */

int SMDR_Read_Values (char *file, int n, char *whichVar[])
{
  int j, success = 1;
  FILE *fp;
  char funcname[] = "SMDR_Read_Values";

  if ((fp = fopen (file, "r")) == NULL)
    SMDR_Error (funcname, "Specified input file does not exist.", 123);

  for (j=0; j<n; j++)
    success *= SMDR_Get_Value (fp, whichVar[j]);

  fclose (fp);

  if (success == 0)
    SMDR_Error (funcname, "Some input data is missing.", 111);

  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Allows for interactive setting of a list of values. */

int SMDR_Set_Values_Interactively (int n, char *whichVar[])
{
  int i, j, success = 1;
  char tmpword[50];
  char funcname[] = "SMDR_Set_Values_Interactively";

  printf ("Entering Interactive Mode.\n");
  printf ("At each prompt, enter a new value or [Return] to accept.\n");

  for (j=0; j<n; j++) {

    /* Which one is it? */
    for (i=0; i<nVals; i++)
      if (strcmp(whichVar[j], varName[i]) == 0)
	break;

    if (i == nVals)
      SMDR_Error (funcname, "Variable not found!", 777);

    printf ("%s = ? [%.10Lf] ", whichVar[j], *varValue[i]);

    fgets (tmpword, sizeof(tmpword), stdin);

    if (!strcmp (tmpword, "\n"))
      ;
    else {
      *varValue[i] = (SMDR_REAL) strtold (tmpword, NULL);

    }
    /* Just to check: */
    /* printf("%s = %Lf\n", whichVar[j], *varValue[i]); */
  }

  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Reads all values from an input file. */

int SMDR_Read_Model_File (char *infile)
{
  FILE *fp;
  int nVars = 35;
  int success = 1;
  char funcname[] = "SMDR_Read_Model_File";

  char *varList[] = {
    "SMDR_Q_in",
    "SMDR_g3_in",
    "SMDR_gp_in",
    "SMDR_g_in",
    "SMDR_yt_in",
    "SMDR_yb_in",
    "SMDR_yc_in",
    "SMDR_ys_in",
    "SMDR_yu_in",
    "SMDR_yd_in",
    "SMDR_ytau_in",
    "SMDR_ymu_in",
    "SMDR_ye_in",
    "SMDR_lambda_in",
    "SMDR_m2_in",
    "SMDR_v_in",
    "SMDR_Lambda_in",
    "SMDR_Delta_alpha_had_5_MZ_in",
    "SMDR_Mt_pole",
    "SMDR_Mh_pole",
    "SMDR_MZ_BreitWigner",
    "SMDR_MZ_pole",
    "SMDR_MW_BreitWigner",
    "SMDR_MW_pole",
    "SMDR_GFermi",
    "SMDR_alpha",
    "SMDR_alphaS_5_MZ",
    "SMDR_mbmb",
    "SMDR_mcmc",
    "SMDR_ms_2GeV",
    "SMDR_md_2GeV",
    "SMDR_mu_2GeV",
    "SMDR_Mtau_pole",
    "SMDR_Mmuon_pole",
    "SMDR_Melectron_pole"};

  if ((fp = fopen(infile, "r")) == NULL)
    SMDR_Error (funcname, "Specified input file does not exist.", 123);

  success = SMDR_Read_Values (infile, nVars, varList);

  fclose(fp);

  if (success == 0)
    SMDR_Error (funcname, "Some input data is missing.", 111);

  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Reads MSbar input parameters from a file. */

int SMDR_Read_MSbar_Inputs (char *infile)
{
  FILE *fp;
  int nVars = 14;
  int success;
  char funcname[] = "SMDR_Read_MSbar_Inputs";

  char *varList[] = {
    "SMDR_Q_in",
    "SMDR_g3_in",
    "SMDR_gp_in",
    "SMDR_g_in",
    "SMDR_yt_in",
    "SMDR_yb_in",
    "SMDR_yc_in",
    "SMDR_ys_in",
    "SMDR_yu_in",
    "SMDR_yd_in",
    "SMDR_ytau_in",
    "SMDR_ymu_in",
    "SMDR_ye_in",
    "SMDR_lambda_in"
  };

  if ((fp = fopen (infile, "r")) == NULL)
    SMDR_Error (funcname, "Specified input file does not exist.", 123);

  success = SMDR_Read_Values (infile, nVars, varList);

  fclose(fp);

  if (success == 0)
    SMDR_Error (funcname, "Some input data is missing.", 111);

  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Reads on-shell input parameters from a file. */

int SMDR_Read_OS_Inputs (char *infile)
{
  FILE *fp;
  int nVars = 15;
  int success;
  char funcname[] = "SMDR_Read_OS_Inputs";

  char *varList[] = {
    "SMDR_Delta_alpha_had_5_MZ_in",
    "SMDR_Mt_pole",
    "SMDR_Mh_pole",
    "SMDR_MZ_BreitWigner",
    "SMDR_GFermi",
    "SMDR_alpha",
    "SMDR_alphaS_5_MZ",
    "SMDR_mbmb",
    "SMDR_mcmc",
    "SMDR_ms_2GeV",
    "SMDR_md_2GeV",
    "SMDR_mu_2GeV",
    "SMDR_Mtau_pole",
    "SMDR_Mmuon_pole",
    "SMDR_Melectron_pole"};

  if ((fp = fopen(infile, "r")) == NULL)
    SMDR_Error (funcname, "Specified input file does not exist.", 123);

  success = SMDR_Read_Values (infile, nVars, varList);

  fclose(fp);

  if (success == 0)
    SMDR_Error (funcname, "Some input data is missing.", 111);
  
  return success;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int SMDR_Write_Model_File (char *outfile)
{
  FILE *fp;
  char funcname[] = "SMDR_Write_Model_File";

  if ((fp = fopen(outfile, "w")) == NULL)
    SMDR_Error (funcname, "Output file cannot be opened.", 246);

  fprintf(fp, "# %s Version %s\n\n", SMDR_NAME, SMDR_VERSION);

  fprintf(fp, "SMDR_Q_in  = %.8Lf;\n", SMDR_Q_in);
  fprintf(fp, "SMDR_g3_in = %.16Lf;\n", SMDR_g3_in);
  fprintf(fp, "SMDR_g_in  = %.16Lf;\n", SMDR_g_in);
  fprintf(fp, "SMDR_gp_in = %.16Lf;\n", SMDR_gp_in);
  fprintf(fp, "SMDR_v_in  = %.14Lf;\n", SMDR_v_in);
  fprintf(fp, "SMDR_lambda_in  = %.16Lf;\n", SMDR_lambda_in);
  fprintf(fp, "SMDR_yt_in = %.16Lf;\n", SMDR_yt_in);
  fprintf(fp, "SMDR_yb_in = %.16Lf;\n", SMDR_yb_in);
  fprintf(fp, "SMDR_yc_in = %.18Lf;\n", SMDR_yc_in);
  fprintf(fp, "SMDR_ys_in = %.18Lf;\n", SMDR_ys_in);
  fprintf(fp, "SMDR_yd_in = %.20Lf;\n", SMDR_yd_in);
  fprintf(fp, "SMDR_yu_in = %.20Lf;\n", SMDR_yu_in);
  fprintf(fp, "SMDR_ytau_in = %.16Lf;\n", SMDR_ytau_in);
  fprintf(fp, "SMDR_ymu_in  = %.17Lf;\n", SMDR_ymu_in);
  fprintf(fp, "SMDR_ye_in   = %.20Lf;\n", SMDR_ye_in);
  fprintf(fp, "\n");
  fprintf(fp, "SMDR_Delta_alpha_had_5_MZ_in = %.12Lf;\n", SMDR_Delta_alpha_had_5_MZ_in);
  fprintf(fp, "\n");
  fprintf(fp, "SMDR_m2_in = %.12Lf;\n",SMDR_m2_in);
  fprintf(fp, "SMDR_Lambda_in = %.6Lf;\n",SMDR_Lambda_in);
  fprintf(fp, "\n");
  fprintf(fp, "SMDR_Mt_pole = %.14Lf;\n", SMDR_Mt_pole);
  fprintf(fp, "SMDR_Mh_pole = %.14Lf;\n", SMDR_Mh_pole);
  fprintf(fp, "SMDR_MZ_BreitWigner = %.14Lf;\n", SMDR_MZ_BreitWigner);
  fprintf(fp, "SMDR_MZ_pole = %.14Lf;\n", SMDR_MZ_pole);
  fprintf(fp, "SMDR_MW_BreitWigner = %.14Lf;\n", SMDR_MW_BreitWigner);
  fprintf(fp, "SMDR_MW_pole = %.14Lf;\n", SMDR_MW_pole);
  fprintf(fp, "SMDR_GFermi = %.20Lf;\n", SMDR_GFermi);
  fprintf(fp, "SMDR_alpha = 1/%.16Lf;\n", 1./SMDR_alpha);
  fprintf(fp, "SMDR_alphaS_5_MZ = %.16Lf;\n", SMDR_alphaS_5_MZ);
  fprintf(fp, "SMDR_mbmb = %.16Lf;\n", SMDR_mbmb);
  fprintf(fp, "SMDR_mcmc = %.16Lf;\n", SMDR_mcmc);
  fprintf(fp, "SMDR_ms_2GeV = %.16Lf;\n", SMDR_ms_2GeV);
  fprintf(fp, "SMDR_md_2GeV = %.16Lf;\n", SMDR_md_2GeV);
  fprintf(fp, "SMDR_mu_2GeV = %.16Lf;\n", SMDR_mu_2GeV);
  fprintf(fp, "SMDR_Mtau_pole = %.16Lf;\n", SMDR_Mtau_pole);
  fprintf(fp, "SMDR_Mmuon_pole = %.17Lf;\n", SMDR_Mmuon_pole);
  fprintf(fp, "SMDR_Melectron_pole = %.19Lf;\n", SMDR_Melectron_pole);

  /* fprintf(fp, "SMDR_Q_in = %.15Le;\n", (long double) SMDR_Q_in);  */
  /* fprintf(fp, "SMDR_g3_in = %.15Le;\n", (long double) SMDR_g3_in);  */
  /* fprintf(fp, "SMDR_gp_in = %.15Le;\n", (long double) SMDR_gp_in);  */
  /* fprintf(fp, "SMDR_g_in = %.15Le;\n", (long double) SMDR_g_in);  */
  /* fprintf(fp, "SMDR_yt_in = %.15Le;\n", (long double) SMDR_yt_in);  */
  /* fprintf(fp, "SMDR_yb_in = %.15Le;\n", (long double) SMDR_yb_in);  */
  /* fprintf(fp, "SMDR_yc_in = %.15Le;\n", (long double) SMDR_yc_in);  */
  /* fprintf(fp, "SMDR_ys_in = %.15Le;\n", (long double) SMDR_ys_in);  */
  /* fprintf(fp, "SMDR_yu_in = %.15Le;\n", (long double) SMDR_yu_in);  */
  /* fprintf(fp, "SMDR_yd_in = %.15Le;\n", (long double) SMDR_yd_in);  */
  /* fprintf(fp, "SMDR_ytau_in = %.15Le;\n", (long double) SMDR_ytau_in);  */
  /* fprintf(fp, "SMDR_ymu_in = %.15Le;\n", (long double) SMDR_ymu_in);  */
  /* fprintf(fp, "SMDR_ye_in = %.15Le;\n", (long double) SMDR_ye_in);  */
  /* fprintf(fp, "SMDR_lambda_in = %.15Le;\n", (long double) SMDR_lambda_in);  */
  /* fprintf(fp, "SMDR_m2_in = %.15Le;\n", (long double) SMDR_m2_in);  */
  /* fprintf(fp, "SMDR_v_in = %.15Le;\n", (long double) SMDR_v_in);  */
  /* fprintf(fp, "SMDR_Lambda_in = %.15Le;\n", (long double) SMDR_Lambda_in);  */
  /* fprintf(fp, "SMDR_Delta_alpha_had_5_MZ_in = %.15Le;\n", (long double) SMDR_Delta_alpha_had_5_MZ_in);  */
  /* fprintf(fp, "SMDR_Mt_pole = %.15Le;\n", (long double) SMDR_Mt_pole);  */
  /* fprintf(fp, "SMDR_Mh_pole = %.15Le;\n", (long double) SMDR_Mh_pole);  */
  /* fprintf(fp, "SMDR_MZ_BreitWigner = %.15Le;\n", (long double) SMDR_MZ_BreitWigner);  */
  /* fprintf(fp, "SMDR_MZ_pole = %.15Le;\n", (long double) SMDR_MZ_pole);  */
  /* fprintf(fp, "SMDR_MW_BreitWigner = %.15Le;\n", (long double) SMDR_MW_BreitWigner);  */
  /* fprintf(fp, "SMDR_MW_pole = %.15Le;\n", (long double) SMDR_MW_pole);  */
  /* fprintf(fp, "SMDR_GFermi = %.15Le;\n", (long double) SMDR_GFermi);  */
  /* fprintf(fp, "SMDR_alpha = %.15Le;\n", (long double) SMDR_alpha);  */
  /* fprintf(fp, "SMDR_alphaS_5_MZ = %.15Le;\n", (long double) SMDR_alphaS_5_MZ);  */
  /* fprintf(fp, "SMDR_mbmb = %.15Le;\n", (long double) SMDR_mbmb);  */
  /* fprintf(fp, "SMDR_mcmc = %.15Le;\n", (long double) SMDR_mcmc);  */
  /* fprintf(fp, "SMDR_ms_2GeV = %.15Le;\n", (long double) SMDR_ms_2GeV);  */
  /* fprintf(fp, "SMDR_md_2GeV = %.15Le;\n", (long double) SMDR_md_2GeV);  */
  /* fprintf(fp, "SMDR_mu_2GeV = %.15Le;\n", (long double) SMDR_mu_2GeV);  */
  /* fprintf(fp, "SMDR_Mtau_pole = %.15Le;\n", (long double) SMDR_Mtau_pole);  */
  /* fprintf(fp, "SMDR_Mmuon_pole = %.15Le;\n", (long double) SMDR_Mmuon_pole);  */
  /* fprintf(fp, "SMDR_Melectron_pole = %.15Le;\n", (long double) SMDR_Melectron_pole);  */

  fclose (fp);

  return 0;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Displays version information. */

void SMDR_Write_Version (FILE *fp, char *prepend)
{
  fprintf (fp, "%s%s v%s\n", prepend, SMDR_NAME, SMDR_VERSION);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Displays version information to stdout. */

void SMDR_Display_Version ()
{
  SMDR_Write_Version (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Displays version information. */

void SMDR_Write_Column_Data (FILE *fp, int n, char *desc[], char *prepend)
{
  int i;

  fprintf (fp, "%s\n", prepend);
  fprintf (fp, "%sCOLUMN DESCRIPTIONS:\n", prepend);

  for (i=0; i<n; i++)
    fprintf (fp, "%s%s\n", prepend, desc[i]);

  fprintf (fp, "%s\n", prepend);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Displays basic MSbar parameter values. */

void SMDR_Write_MSbar_Parameters (FILE *fp, char *prepend)
{
  fprintf (fp, "%sQ = ", prepend);
  if (Q < 100000)
    fprintf(fp, "%.6Lf;\n", Q);
  else
    fprintf(fp, "%.6e;\n", (double) Q);

  fprintf(fp, "%sgauge  couplings: g3 = %Lf;      g = %Lf;      gp = %Lf;\n",
	  prepend, SMDR_g3, SMDR_g, SMDR_gp);
  fprintf(fp, "%sYukawa couplings: yt = %Lf;     yb = %Lf;    ytau = %.8Lf;\n",
	  prepend, SMDR_yt,SMDR_yb, SMDR_ytau);
  fprintf(fp, "%s                  yc = %.7Lf;    ys = %.8Lf;   ymu = %.9Lf;\n",
	  prepend, SMDR_yc,SMDR_ys, SMDR_ymu);
  fprintf(fp, "%s                  yu = %.10Lf; yd = %.9Lf;   ye = %.11Lf;\n",
	  prepend, SMDR_yu,SMDR_yd, SMDR_ye);
  fprintf(fp, "%sHiggs self-coupling lambda = %Lf;\n", prepend, SMDR_lambda);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Displays basic MSbar parameter values. */

void SMDR_Display_MSbar_Parameters (void)
{
  SMDR_Write_MSbar_Parameters (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Write_v (FILE *fp, char *prepend)
{
  fprintf(fp, "%sHiggs VEV = %Lf;\n", prepend, SMDR_v);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Display_v (void)
{
  SMDR_Write_v (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Write_m2 (FILE *fp, char *prepend)
{
  /* int sign = 1; */
  /* if (SMDR_m2 < 0) sign = -1; */

  fprintf(fp, "%sHiggs mass^2 parameter m2 = %.6Lf = ", prepend, SMDR_m2);
  if (SMDR_m2 < 0) 
    fprintf(fp, "-(%.6Lf)^2;\n",SMDR_SQRT(-SMDR_m2));
  else
    fprintf(fp, "(%.6Lf)^2;\n",SMDR_SQRT(SMDR_m2));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Display_m2 (void)
{
  SMDR_Write_m2 (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Write_Lambda (FILE *fp, char *prepend)
{
  fprintf(fp, "%sTree-level Vacuum Energy Lambda = %.6Lf = ",
	  prepend, SMDR_Lambda);

  if (SMDR_Lambda < 0) 
    fprintf(fp, "-(%.6Lf)^4;\n",SMDR_SQRT(SMDR_SQRT(-SMDR_Lambda)));
  else
    fprintf(fp, "(%.6Lf)^4;\n",SMDR_SQRT(SMDR_SQRT(SMDR_Lambda)));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Display_Lambda (void)
{
  SMDR_Write_Lambda (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Write_Delta_alpha_had5 (FILE *fp, char *prepend)
{
  fprintf(fp, "%sDelta_hadronic^(5) alpha(MZ) = %Lf\n",
	  prepend, SMDR_Delta_alpha_had_5_MZ);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Display_Delta_alpha_had5 (void)
{
  SMDR_Write_Delta_alpha_had5 (stdout, "");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Write_OS_Inputs (FILE *fp, char *prepend)
{
  fprintf(fp, "%sMt = %.3Lf;\n", prepend, SMDR_Mt_pole);
  fprintf(fp, "%sMh = %.3Lf;\n", prepend, SMDR_Mh_pole);
  fprintf(fp, "%sMZ = %.5Lf;\n", prepend, SMDR_MZ_BreitWigner);
  fprintf(fp, "%salpha_S_5_MZ = %Lf;\n", prepend, SMDR_alphaS_5_MZ);
  fprintf(fp, "%sGFermi = %.8Lf 10^-5;\n", prepend, 100000 * SMDR_GFermi);
  fprintf(fp, "%salpha = 1/%.8Lf;\n", prepend, 1./SMDR_alpha);
  fprintf(fp, "%sDelta_hadronic^(5) alpha(MZ) = %.6Lf;\n",
	  prepend, SMDR_Delta_alpha_had_5_MZ_in);

  fprintf(fp, "%smb(mb) = %.5Lf;      ", prepend, SMDR_mbmb);
  fprintf(fp, "mc(mc) = %.5Lf;\n", SMDR_mcmc);

  fprintf(fp, "%sms(2 GeV) = %.5Lf;   ", prepend, SMDR_ms_2GeV);
  fprintf(fp, "mu(2 GeV) = %.6Lf;  ", SMDR_mu_2GeV);
  fprintf(fp, "md(2 GeV) = %.6Lf;\n", SMDR_md_2GeV);

  fprintf(fp, "%sMtau = %.6Lf;       ", prepend, SMDR_Mtau_pole);
  fprintf(fp, "Mmuon = %.8Lf;    ", SMDR_Mmuon_pole);
  fprintf(fp, "Melectron = %.12Lf;\n", SMDR_Melectron_pole);
  
  return;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Display_OS_Inputs (void)
{
  SMDR_Write_OS_Inputs (stdout, "");
}
