#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "/Users/tinker/src/lib/cutil.h"

// how do I run this? 
/* 
 * use the file input.bat, which has the COSMOMC file names and the columns to be used. argv[1]
 * use the hypercube output, ie MCMC hy[percube solution. argv[2]
 */

void jacobi(float **a, int n, float d[], float **v, int *nrot);
void gaussj(float **a, int n, float **b, int m);

int main(int argc, char **argv)
{
  int nfiles, i, np, pindx[100], npoints=0, j, k, offset, n1, nrot, n, n0, nsamp, ii;
  FILE *fp, *fp1;
  float xx[100], **x, xbar[100], **cov1, **cov, **evect, *eval, **tmp1, dx[100], da[100], *weight,
    wtot=0,dr;
  char fname[1000],aa[1000];
  float a, b, theta;

  fp = openfile(argv[1]);
  nfiles = filesize(fp);
  fprintf(stderr,"%d files in [%s]\n",nfiles,argv[1]);

  // first, lets get the number of parameters and the chain sizes
  npoints = 0;
  for(i=1;i<=nfiles;++i)
    {
      fscanf(fp,"%s",fname); // chain filename
      fscanf(fp,"%d",&np); // how many parameters? (better be same!)
      fp1 = openfile(fname);
      npoints += filesize(fp1);
      fclose(fp1);
      fgets(aa,1000,fp);
    }
  x = matrix(1,npoints,1,np);
  weight = vector(1,npoints);
  rewind(fp);
  
  offset=0;
  for(i=1;i<=nfiles;++i)
    {
      fscanf(fp,"%s",fname); // chain filename
      fscanf(fp,"%d",&np); // how many parameters?
      for(j=1;j<=np;++j) // which parameters?
	fscanf(fp,"%d",&pindx[j]);
      
      fp1 = openfile(fname);
      n1 = filesize(fp1);
      if(i==1) n0 = n1;
      for(j=1;j<=n1;++j)
	{
	  weight[j+offset] = n0*1./n1; //acount for fact that some chains have more points
	  for(k=1;k<=pindx[np];++k)
	    fscanf(fp1,"%f",&xx[k]);
	  fgets(aa,1000,fp1);
	  for(k=1;k<=np;++k)
	    x[j+offset][k] = xx[pindx[k]]; 
	  wtot += weight[j+offset];
	}
      fclose(fp1);
      offset = n1;
    }
  fprintf(stderr,"%d %f %f\n",npoints, wtot, weight[i]);


  // subtract off the means of each parameter
  for(i=1;i<=np;++i)
    xbar[i] = 0;
  for(i=1;i<=npoints;++i)
    for(j=1;j<=np;++j)
      xbar[j] += x[i][j]*weight[i]/wtot;
  for(i=1;i<=npoints;++i)
    for(j=1;j<=np;++j)
      x[i][j] -= xbar[j];

  // now get covariance matrix for parameters
  cov1 = matrix(1,np,1,np);
  cov  = matrix(1,np,1,np);
  evect = matrix(1,np,1,np);
  eval = vector(1,np);
  tmp1 = matrix(1,np,1,1);

  for(i=1;i<=npoints;++i)
    {
      for(j=1;j<=np;++j)
	{
	  for(k=1;k<=np;++k)
	    cov1[j][k]+=x[i][j]*x[i][k]*weight[i]*weight[i];
	}
    }  
  for(i=1;i<=np;++i)
    for(j=1;j<=np;++j)
      cov[i][j] = cov1[i][j] = cov1[i][j]/wtot;
  
  jacobi(cov1,np,eval,evect,&nrot);
  gaussj(evect,np,tmp1,1);

  // lets make ellipses
  i = 2;
  j = 5;
  a = (cov[i][i] + cov[j][j])/2 + sqrt((cov[i][i]-cov[j][j])*(cov[i][i]-cov[j][j])/4+cov[i][j]);
  a = sqrt(fabs(a));
  b = (cov[i][i] + cov[j][j])/2 - sqrt((cov[i][i]-cov[j][j])*(cov[i][i]-cov[j][j])/4+cov[i][j]);
  b = sqrt(fabs(b));
  theta = 0.5*atan(cov[i][j]*2/(cov[i][i]-cov[j][j]));
  printf("%e %e %e %e %e %e %e\n",a,b,theta,xbar[i],xbar[j],sqrt(cov[i][i]),sqrt(cov[j][j]));

  // print out the mean
  for(i=1;i<=np;++i)
    printf("%e ",xbar[i]);
  printf("\n");

  //print out the eigenvectors and values
  for(i=1;i<=np;++i)
    {
      printf("%e ",sqrt(eval[i]));
      for(j=1;j<=np;++j)
	printf("%e ",evect[i][j]);
      printf("\n");
    }

  // now lets step along each eigenvector individually
  for(n=1;n<=10000;++n)
    {
      
      for(i=1;i<=np;++i)
	da[i] = (2*drand48()-1)*3;

      /* get the distance in the eigenspace
       */
      dr = 0;
      for(i=1;i<=np;++i)
	dr += da[i]*da[i];
      dr = sqrt(dr);
      /*      
      for(i=1;i<=np;++i)
	da[i] = 0;
      da[n] = 1;
      */
      for(i=1;i<=np;++i)
	for(dx[i]=0,j=1;j<=np;++j)
	  dx[i] += da[j]*sqrt(eval[j])*evect[j][i];      

      for(i=1;i<=np;++i)
	dx[i] += xbar[i];
      
      for(i=1;i<=np;++i)
	printf("%e ",dx[i]);
      printf("%e\n",dr);
    }

  // now lets read in the hypercube output and transform things
  // fp = fopen("lh.out","r");
  fp = fopen(argv[2],"r");
  fprintf(stderr,"Taking hypercube input from [%s]\n",argv[2]);
  nsamp = filesize(fp);
  fp1 = fopen("LH_eigenspace.out","w");
  for(ii=1;ii<=nsamp;++ii)
    {
      for(j=1;j<=np;++j)
	fscanf(fp,"%f",&xx[j]);

      for(i=1;i<=np;++i)
	for(dx[i]=0,j=1;j<=np;++j)
	  dx[i] += 4*(2*xx[j]-1)*sqrt(eval[j])*evect[j][i];
      
      for(i=1;i<=np;++i)
	dx[i] += xbar[i];

      // get the extra point for the sum of nuetrino masses
      /*
      fscanf(fp,"%f",&xx[0]);
      xx[0] = xx[0]*1.7 + 2.6;
      */
      
      for(i=1;i<=np;++i)
	fprintf(fp1,"%e ",dx[i]);
      //fprintf(fp1,"%e",xx[0]);
      fprintf(fp1,"\n");
    }

}
