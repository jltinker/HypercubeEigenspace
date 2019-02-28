#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"

#define POWELL 0

float chi2(float *a);
float gasdev(long *idum);

int NP;

int main(int argc, char **argv)
{
  int i,j, ncf=6, niter, ntrial = 0, naccept = 0;
  float *a,**pp,*yy,FTOL=1.0E-5,chi2min,s1,dlogm,m, d[20], *aprev,chi2prev, chi2i, var[100], hotfac=10;
  long IDUM=555;

  NP = ncf;

  a = vector(1,ncf);
  aprev = vector(1,ncf);

  if(argc>1)
    {
      for(i=1;i<=ncf;++i)
	a[i] = atof(argv[i]);
    }
  else
    {
      for(i=1;i<=ncf;++i)
	a[i] = 1;
    }

  for(i=1;i<=NP;++i)
    var[i] = 0.1;

  for(i=1;i<=NP;++i)
    aprev[i] = a[i];

  chi2prev = chi2(a);

  while(1)
    {
      for(i=1;i<=NP;++i)
	a[i] = aprev[i] + gasdev(&IDUM)*var[i];
      chi2i = chi2(a);

      ntrial++;
      if(!(chi2i<chi2prev || drand48() <= exp(-hotfac*(chi2i-chi2prev)/2)))
	{
	  continue;
	}
      naccept++;
      for(i=1;i<=NP;++i)
	aprev[i] = a[i];
      chi2prev = chi2i;

      printf("ACCEPT %d %d %e ",naccept, ntrial, chi2i);
      for(i=1;i<=NP;++i)
	printf("%e ",a[i]);
      printf("\n");
      fflush(stdout);
    }
}

float chi2(float *a)
{
  static int iter=0;
  char aa[1000];
  FILE *fp;
  float chi2;
  if(NP<=5)
    sprintf(aa,"~/src/gaussian_process/gp_emulator2 %f %f %f %f %f > out.xxy",a[1],a[2],a[3],a[4],a[5]);
  if(NP==6)
    sprintf(aa,"~/src/gaussian_process/gp_emulator2b %f %f %f %f %f %f > out.xxy",
	    a[1],a[2],a[3],a[4],a[5],a[6]);
  if(NP==9)
    sprintf(aa,"~/src/gaussian_process/gp_emulator_hod %f %f %f %f %f %f %f %f %f> out.xxy",
	    a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);
  system(aa);
  fp = fopen("out.xxy","r");
  fscanf(fp,"%f",&chi2);
  fclose(fp);
  if(NP<=5)
    printf("ITER %d %e %f %f %f %f %f\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5]);
  if(NP==6)
    printf("ITER %d %e %f %f %f %f %f %f\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5],a[6]);
  if(NP==9)
    printf("ITER %d %e %f %f %f %f %f %f %f %f %f\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);
  fflush(stdout);
  return chi2;
}
