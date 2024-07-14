//#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>


void PV_update(){

int n;


 for(n=1;n<=nAtom;n++){

   ax[n] = -Q[n]*ef_x[n];
   //ay[n] = -Q[n]*ef_y[n];


 }

   for(n=1;n<=nAtom;n++){

   vx[n]  += deltaT * ax[n];
   rx[n]  += deltaT * vx[n];

   }


   double v, vv = 0;

   for(n=1;n<=nAtom;n++){

   v = vx[n] - 0.5 * deltaT * ax[n];

   vv += v*v;


   }

   KE = 0.5*vv;

   KE = KE/nAtom;

   PE = 0;

   for(n=1;n<=nAtom;n++){

   PE += (ef_x[n])*(ef_x[n]);// + (ef_y[n])*(ef_y[n]);

   }

   PE = (0.5*PE)/nAtom;

   //printf("pe = %lf %lf\n", PE, timeNow);


   for(n=1;n<=nAtom;n++){

   rx[n] -= region[1]*rint(rx[n]/region[1]);
   //ry[n] -= region[2]*rint(ry[n]/region[2]);

   }


}
