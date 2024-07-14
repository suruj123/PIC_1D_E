#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>
#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9

#define MASTER 0
void Close(){

//fclose(fpresult);
 //fclose(fpxyz);
  
  //#pragma acc exit data copyout(rx[0:3], ry[0:3], ax[0:3], ay[0:3])
  //#pragma acc exit data copyout(region[0:3], regionH[0:3])
  free(rx);
  free(ry);
  free(vx);
  free(vy);
  free(ax);
  free(ay);
  free(xs);
  free(ys);
  free(xg);
  free(yg);
  free(qg);
  //free(rho);
  free(Q);
  free(ef_x);
  free(ef_y);
  free(ef_xg);
  free(ef_yg);
  //free(psi);
  free(psinew);
  free(dx1);
  

    //fftw_destroy_plan(forward_plan);
    //fftw_destroy_plan(backward_plan);
 
    fftw_free(rho);
    fftw_free(rho_k);
    fftw_free(phi_k);
    fftw_free(phi);



}
