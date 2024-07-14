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

 void Interpolation2(){


 int i, j, n, jp;


 for(n=1;n<=nAtom;n++){

 i = (int)(xs[n]/del_x) + 1; 
 
 j = i+1;

 //jp = (j+1) % Nx;
 
 if(j > Nx){

 j = 1;

 printf("hit\n");

 }

 //j = (int)(ys[n]/del_y) + 1; 

double  dx =  xg[j] - xs[n];

double  dy =  xs[n] - xg[i];
 
//double  tx = del_x - dx;
//double  ty = del_y - dy;


 ef_x[n] = ef_xg[i]*dx*h_area + ef_xg[j]*dy*h_area; 

/*
 ef_x[n] = ef_xg[(j-1)*Nx + i] * tx * ty * h_area + 
	   ef_xg[(j-1)*Nx + i + 1] * dx * ty * h_area + 
	   ef_xg[(j)*Nx + i] * tx * dy * h_area +
	   ef_xg[(j)*Nx + i + 1] * dx * dy * h_area;


 ef_y[n] = ef_yg[(j-1)*Nx + i] * tx * ty * h_area + 
	   ef_yg[(j-1)*Nx + i + 1] * dx * ty * h_area + 
	   ef_yg[(j)*Nx + i] * tx * dy * h_area +
	   ef_yg[(j)*Nx + i + 1] * dx * dy * h_area;

*/

 }


 printf("done Interp_2");


}
