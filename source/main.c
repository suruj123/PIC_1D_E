//#include "mpi.h"
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

/*const gsl_rng_type * T;
gsl_rng * rnd;*/

void Init();
void AllocArrays();
void Interpolation1();
void PS();
void Efield();
void Interpolation2();
void PV_update();
void Close();



int main(int argc, char **argv){


//time_t t1, t2;


if(taskid == MASTER){

sprintf(dirprefix,"../output/");
prefix = strcat(dirprefix, argv[1]);


sprintf(result, "%s.result", prefix);
fpresult = fopen (result, "w");


sprintf (xyz, "%s.xyz", prefix);
fpxyz = fopen (xyz, "w");








 

for(stepCount=1;stepCount<=stepLimit;stepCount++){
 
  


timeNow=stepCount*deltaT;




Interpolation1();     //particle to mesh

PS();                 //Poison solver

Efield();             //Electric field determination

Interpolation2();     //Mesh to particle

PV_update();          //Prticle push



int n;
// Print the trajectories

if(stepCount % stepTrajectory == 0){
fprintf(fpxyz, "%d\n", nAtom);
fprintf(fpxyz, "timeNow %lf region[1] %lf \n", timeNow, region[1]);

for(n = 1 ; n <= nAtom ; n ++){
fprintf(fpxyz, "%d\t %lf\t %lf\t %lf\n", n, rx[n], vx[n], Q[n]);
fprintf(fpxyz, "\n");
}}


if(stepCount % stepDump == 0){
  char DUMP[256];
  FILE *fpDUMP;
  sprintf (DUMP, "%s.STATE", prefix);
  fpDUMP = fopen (DUMP, "w");

  fprintf(fpDUMP, "timeNow %lf\n", timeNow);
  fprintf(fpDUMP, "nAtom %d\n", nAtom);
  fprintf(fpDUMP, "region[1] %lf\n", region[1]);

  for(n = 1 ; n <= nAtom ; n ++)
  fprintf(fpDUMP, "%d\t %lf\t %lf\t %lf\n", n, rx[n], vx[n], Q[n]);
  fclose(fpDUMP);
}



//Calculate total energy
  aPE += PE;

  aKE += KE;

  if(stepCount % stepAvg == 0){

   aPE = aPE/stepAvg;

   aKE = aKE/stepAvg;


   double TE = aPE + aKE;

  fprintf(fpresult,"%lf\t %lf\t %lf\t %lf\n", timeNow, aKE, aPE, TE);
  fflush(fpresult);

  aKE = 0;
  aPE = 0;

  }


  
  fclose(fpxyz);
  fclose(fpresult);

 


Close();

return 0;
}
