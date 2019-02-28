#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"

#define POWELL 0

float chi2(float *a);

void amoeba(float **p, float y[], int ndim, float ftol,
	    float (*funk)(float []), int *nfunk);
void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	    float (*func)(float []));

int NP;

int main(int argc, char **argv)
{
  int i,j, ncf=6, niter;
  float *a,**pp,*yy,FTOL=1.0E-5,chi2min,s1,dlogm,m, d[20];

  NP = ncf;

  a = vector(1,ncf);
  if(POWELL)
    pp=matrix(1,ncf,1,ncf);
  else
    pp=matrix(1,ncf+1,1,ncf);
  yy = vector(1,ncf+1);

  for(j=1;j<=ncf;++j)
    d[j] = 0.1;
  for(j=1;j<=ncf;++j)
    a[j] = 1.0;
  // if input from cmd line, take that in
  if(ncf>6)
    {
      for(j=6;j<=ncf;++j)
	a[j] = 0.1;
    }
  if(argc>1)
    {
      for(i=1;i<=ncf;++i)
	a[i] = atof(argv[i]);
    }


  for(j=1;j<=ncf;++j)
    pp[1][j]=a[j];
  yy[1]=chi2(a);
      

  if(POWELL)
    {
      for(i=1;i<=ncf;++i)
	{
	  for(j=1;j<=ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(i=1;i<=ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	  yy[i+1]=chi2(a);	  
	  for(j=1;j<=ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[ncf]-=d[ncf];
    }

  if(POWELL) 
    {
      powell(a,pp,ncf,FTOL,&niter,&chi2min,chi2);
    }
  else
    {
      amoeba(pp,yy,ncf,FTOL,chi2,&niter);
      for(i=1;i<=ncf;++i)a[i]=pp[1][i];
    }

  fprintf(stderr,"DONE!\n");
  for(i=1;i<=ncf;++i)printf("%e\n",a[i]);
  fflush(stdout);
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
    sprintf(aa,"~/src/gaussian_process/gp_emulator4 %e %e %e %e %e %e > out.xxy",
	    a[1],a[2],a[3],a[4],a[5],a[6]);
  if(NP==9)
    sprintf(aa,"~/src/gaussian_process/gp_emulator_hod %f %f %f %f %f %f %f %f %f> out.xxy",
	    a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);
  fprintf(stderr,"%s\n",aa);
  system(aa);
  fp = fopen("out.xxy","r");
  fscanf(fp,"%f",&chi2);
  fclose(fp);
  if(NP<=5)
    printf("ITER %d %e %e %e %e %e %e\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5]);
  if(NP==6)
    printf("ITER %d %e %e %e %e %e %e %e\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5],a[6]);
  if(NP==9)
    printf("ITER %d %e %e %e %e %e %e %e %e %e %e\n",iter++,chi2,
	   a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);
  fflush(stdout);
  return chi2;
}
