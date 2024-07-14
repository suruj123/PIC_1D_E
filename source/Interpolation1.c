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

 void Interpolation1(){
 
 int i, j, k, n, m;

 double dx, dy, tx, ty;
 
 //del_x = region[1]/(Nx-1);
 //del_y = region[2]/(Ny-1);
 
 //h_area = 1.0/(del_x*del_y);


 printf("interpolation 1 start\n");

 for(n=1;n<=nAtom;n++){
 xs[n] = rx[n] + regionH[1];

 //printf("xs = %lf\n", xs[n]);
 //ys[n] = ry[n] + regionH[2];
 }

 for(n=1;n<=Nx;n++){
 xg[n] = (n-1)*del_x;
 }

/* for(n=1;n<=Ny;n++){
 yg[n] = (n-1)*del_y;
 }
*/ 
 
for(i=1;i<=Nx;i++){
 
   
   qg[i] = 0.0;   

   //printf("qg = %lf\n", qg[n]);

 } 

for(n=1;n<=nAtom;n++){

  i = (int)(xs[n]/del_x) + 1; 
  //j = (int)(ys[n]/del_y) + 1; 

  j = i+1;

 
  tx = xg[j] - xs[n];
  ty = xs[n] - xg[i];
 
 
  qg[i]       =  qg[i]     +  Q[n]*tx*h_area; //(i, j)  

  qg[j]       =  qg[j]     +  Q[n]*ty*h_area; //(i, j)  

  
 }
 

double sum = 0;

FILE *fp;
fp = fopen("rho_info", "w");


for(i=1;i<=Nx;i++){
 
   
   //rho[n] = qg[n]/(del_x * del_y);   
   //rho[n] = qg[n]/(region[1] * region[2]);   
   rho[i] = qg[i]*(512/40000.0);   

   sum += qg[i];
   //printf("qg = %lf\t rho = %lf %d\n", qg[i], rho[i], Nx);
   fprintf(fp, "%lf\n", rho[i]);

}

fclose(fp);

printf("sum = %lf\n", sum);
 
 
}
