#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>



void PS(){

printf("PS start %d\n", Nx);

double sorcoeff = 1.00;
int  t_end;

int  n;

double locfirst;

 for(n=1;n<=Nx;n++){

  psi[n] = 0.0;

  psinew[n] = 0.0;

  //printf("Hi from for loop 2\n");

 }


int iter = 0;
double error = 1.0;
double tol = 1.0e-6;
while(error>tol){
//while(error>tol && iter < t_end)
double sumerr = 0.0;



  for(n=1;n<=Nx;n++){

  if(n==1){
  locfirst = psi[n+1] + psi[Nx];
  }else if(n==(Nx)){
  locfirst = psi[1] + psi[n-1];
  }else{
  locfirst = psi[n+1] + psi[n-1];
  }




  psinew[n] = (1 - sorcoeff)*psi[n] + sorcoeff*0.50*(locfirst  - Sqr(del_x)*(rho[n] - 1));
  sumerr += Sqr(psinew[n] - psi[n]);

  //printf("Hi from for loop 1\n");

  //printf("%lf %lf %d\n", psi[nx], rho[nx], iter);


  }





 //for(n=1;n<=Nx;n++){
 for(n=1;n<=Nx;n++){

  psi[n] = psinew[n];

  //printf("Hi from for loop 2\n");

 }


error = sqrt(sumerr/(Nx));
printf("%0.16lf\t%d\t %lf\n", error, iter, timeNow);
iter++;
} 

printf("PS done\n");


}
