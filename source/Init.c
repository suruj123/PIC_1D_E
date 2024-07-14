#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>
#include <fftw3.h>

#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9

#define MASTER 0

 void Init(){

 char dummy[128];
  FILE *fp;
  fp = fopen("input-data","r");
  fscanf(fp, "%s %lf", dummy, &GAMMA);
  fscanf(fp, "%s %lf", dummy, &rCut);
  fscanf(fp, "%s %lf", dummy, &kappa);
  fscanf(fp, "%s %lf", dummy, &deltaT);
  fscanf(fp, "%s %d", dummy, &stepAvg);
  fscanf(fp, "%s %d", dummy, &stepEquil);
  fscanf(fp, "%s %d", dummy, &stepLimit);
  fscanf(fp, "%s %d", dummy, &stepDump);
  fscanf(fp, "%s %d", dummy, &stepTrajectory);
  fscanf(fp, "%s %c", dummy, &thermo);
  fscanf(fp, "%s %c", dummy, &BC);
  fscanf(fp, "%s %d", dummy, &limitCorrAv);
  fscanf(fp, "%s %d", dummy, &nBuffCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr2);
  fscanf(fp, "%s %d", dummy, &nValCorr);
  fscanf(fp, "%s %d", dummy, &stepCorr);
  fscanf(fp, "%s %d", dummy, &limitAcfAv);
  fscanf(fp, "%s %d", dummy, &nBuffAcf);
  fscanf(fp, "%s %d", dummy, &nValAcf);
  fscanf(fp, "%s %d", dummy, &stepAcf);
  fscanf(fp, "%s %lf", dummy, &rangeRdf);
  fscanf(fp, "%s %d", dummy, &limitRdf);
  fscanf(fp, "%s %d", dummy, &sizeHistRdf);
  fscanf(fp, "%s %d", dummy, &stepRdf);
  fscanf(fp, "%s %d", dummy, &limitVel);
  fscanf(fp, "%s %lf", dummy, &rangeVel);
  fscanf(fp, "%s %d", dummy, &sizeHistVel);
  fscanf(fp, "%s %d", dummy, &stepVel);
  fscanf(fp, "%s %c", dummy, &cellAlgo);
  fscanf(fp, "%s %lf", dummy, &lambda_c);
  fscanf(fp, "%s %d", dummy, &Nx);
  fscanf(fp, "%s %d", dummy, &Ny);
  fscanf(fp, "%s %d", dummy, &Nz);

  

  fclose(fp);

  FILE *fpSTATE;
  if((fpSTATE = fopen("../STATE","r"))==NULL)
    fprintf(fpresult,"Could not open ../STATE file\n");
  //fscanf(fpSTATE, "%s %lf", dummy, &timeNow);
  //fscanf(fpSTATE, "%s %d", dummy, &nAtom);
  //fscanf(fpSTATE, "%s %lf", dummy, &region[1]);
  //fscanf(fpSTATE, "%s %lf", dummy, &region[2]);

  //density = nAtom/(region[1]*region[2]);

  region[1] = 50.0;

  nAtom = 40000;


  cells[1] = region[1] / rCut;
  //cells[2] = region[2] / rCut;
  //cellList = (int *)malloc((nAtom + cells[1] * cells[2] + 1) * sizeof(int));
  regionH[1] = 0.5*region[1];
 // regionH[2] = 0.5*region[2];
  
  	
  rx = (double*)malloc( (nAtom + 1) * sizeof(double));
  ry = (double*)malloc( (nAtom + 1) * sizeof(double));

  vx = (double*)malloc( (nAtom + 1) * sizeof(double));
  vy = (double*)malloc( (nAtom + 1) * sizeof(double));

  ax = (double*)malloc( (nAtom + 1) * sizeof(double));
  ay = (double*)malloc( (nAtom + 1) * sizeof(double));

  fax = (double*)malloc( (nAtom + 1) * sizeof(double));
  fay = (double*)malloc( (nAtom + 1) * sizeof(double));

  dx1 = (double*)malloc( (nAtom + 1) * sizeof(double));
  p = (double*)malloc( (nAtom + 1) * sizeof(double));

  alpha = (double*)malloc( (nAtom + 1) * sizeof(double));
  Q = (double*)malloc( (nAtom + 1) * sizeof(double));
  Qold = (double*)malloc( (nAtom + 1) * sizeof(double));

  b = (double*)malloc( (nAtom + 1) * sizeof(double));
  a = (double**)malloc( (nAtom + 1) * sizeof(double));


  //Charge

  //lambda_c = 10;

  lambda_c = 1;


  //norm = 1.0/(9*lambda_c*lambda_c*lambda_c*lambda_c);
  //norm = 1.0/(9*lambda_c*lambda_c);
  norm = 1.0;


  //norm = 1.0/(9.0);
  //norm = 1.0;
  //Interpolation
  
  ngrids_X = Nx;
  ngrids_Y = Ny;

  xs = (double*)malloc( (nAtom + 1) * sizeof(double));
  ys = (double*)malloc( (nAtom + 1) * sizeof(double));

  ef_x = (double*)fftw_malloc( (nAtom + 1) * sizeof(double));
  ef_y = (double*)fftw_malloc( (nAtom + 1) * sizeof(double));
  
  xg = (double*)malloc( (Nx + 1) * sizeof(double));
  yg = (double*)malloc( (Ny + 1) * sizeof(double));

  qg = (double*)fftw_malloc( (Nx + 1) * sizeof(double));
  //rho = (double*)malloc( (Nx + 1) * sizeof(double));

  psi = (double*)fftw_malloc( (Nx + 1) * sizeof(double));
  psinew = (double*)fftw_malloc( (Nx + 1) * sizeof(double));

  ef_xg = (double*)malloc( (Nx + 1) * sizeof(double));
  ef_yg = (double*)malloc( (Nx + 1) * sizeof(double));


   rho = (double*) fftw_malloc(sizeof(double) * Nx);
   rho_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
   phi_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
   phi = (double*) fftw_malloc(sizeof(double) * Nx);



  //fftw_plan forward_plan = fftw_plan_dft_r2c_1d(Nx, rho, rho_k, FFTW_ESTIMATE);
  //fftw_plan backward_plan = fftw_plan_dft_c2r_1d(Nx, phi_k, phi, FFTW_ESTIMATE);



  del_x = region[1]/(Nx-1);
  //del_y = region[2]/(Ny-1);

  h_area = 1.0/(del_x);


  int n;


   for(n = 0; n <= nAtom; n ++){

   a[n] = (double*)malloc( (nAtom + 1) * sizeof(double));

   }



  //#pragma acc enter data copyin(rx[0:3], ry[0:3],ax[0:3], ay[0:3])
  //#pragma acc enter data copyin(region[0:3], regionH[0:3])
	
  int idx;
   
    
    for(n = 1; n <= nAtom; n ++){

    Q[n] = 1.0;

    fscanf(fpSTATE, "%lf %lf", &rx[n], &vx[n]);

    }

    fclose(fpSTATE);

  
  for(n = 1; n <= nAtom; n ++){
  //Q[n] = 1.0;
  //printf("%lf %lf %lf %lf %lf taskid=%d\n", dx[n], rx[n], ry[n], vx[n], vy[n],taskid);
  //printf("%lf %lf\n",rx[n], vx[n]);
  }
 
  //printf("%lf %lf\n", region[1], regionH[1]);
}
