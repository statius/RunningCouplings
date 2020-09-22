#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "smdr.h"
#include "smdr_internal.h"

/* 
  This program just calls the function SMDR_Make_Accelcoeffs(). 
  After running it, there will be a new file accelcoeffs.h, overwriting
  the old one, and corresponding to the file ReferenceModel.dat. 
  After running it, one should issue make again. 
*/

int main ()
{
  SMDR_Make_Accelcoeffs();
  return 0;
}

