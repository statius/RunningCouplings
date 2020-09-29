/* RG running of the basic model paramters. */

#include "smdr_internal.h"


/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Loads the MSbar input parameters and runs them to Q.             */
int SMDR_RGeval_SM (SMDR_REAL Q_final, int loopOrder)
{
  SMDR_Load_Inputs();
  return(SMDR_RGrun_SM (Q_final, loopOrder));
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Runs all parameters from Q to Q_final at loop order specified.   */

int SMDR_RGrun_SM (SMDR_REAL Q_final, int loopOrder)
{
  SMDR_REAL t_final, t, dt;
  int min_steps, max_steps;
  int force_step;
  int goodsteps, badsteps; 
  int rk6status; /* 1 for success or forced; 0 for need retry. */
  SMDR_REAL Q_init = Q;
  char funcname[] = "RGrun";

  if ( (loopOrder < 1) || (loopOrder > 5) )
    SMDR_Error (funcname,
    "Invalid loop order specified, should be 0, 1, 2, 3, 4, or 5", 3);

  t_final = TSIL_LOG(Q_final/Q_init);

  if (TSIL_FABS(t_final) < TSIL_TOL) return 0;
  min_steps = 1 + ((int) (5.0L * TSIL_FABS(t_final)));
  max_steps = 1 + ((int) (125.0L * TSIL_FABS(t_final)));

  dt = t_final/(1 + ((int) (124.0L * TSIL_FABS(t_final))));
  t = 0.0L;

  SMDR_Betas (loopOrder);

  goodsteps = badsteps = 0;

  while ( TSIL_FABS(dt) < 0.5 * TSIL_FABS(t_final - t) )
    {
      if ( TSIL_FABS(dt) < TSIL_FABS(t_final/max_steps) )
	{
	  force_step = 1;
	  dt = t_final/((SMDR_REAL) max_steps);
	}
      else force_step = 0;
      
      if ( TSIL_FABS(dt) > TSIL_FABS(t_final)/min_steps )
	dt = (t_final)/min_steps;

      rk6status = SMDR_RG_rk6 (&dt, force_step, loopOrder);
      t = TSIL_LOG (Q/Q_init);
      
      if (1 == rk6status) {
	goodsteps += (1 - force_step);
	badsteps += force_step;
      }
    }

  /* The remaining distance is less than 2.0 times the step size.  So,
     take exactly two more steps. */
  dt = 0.5L*(t_final - t);
  SMDR_RG_rk6 (&dt, 1, loopOrder);
  t = TSIL_LOG (Q/Q_init);

  /* Arrange final step to land exactly on t_final, and force it. */
  dt = t_final - t;
  SMDR_RG_rk6 (&dt, 1, loopOrder);

  SMDR_Update ();

  /*
  printf("goodsteps = %d    badsteps  = %d\n",goodsteps, badsteps);
  */

  /* Return a status code eventually */
  return 0;
}
