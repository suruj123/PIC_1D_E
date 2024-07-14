#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"mpi.h"
#include"global.h"
#include<time.h>
#include<string.h>
#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9

#define MASTER 0

 void Efield(){


  int n, i, j;	 

  double locfirst;

  double h_x = region[1]/(Nx-1);

   printf("Start E_field");
 

  for(n=1;n<=Nx;n++){


	    

      if(n==1){
      locfirst = psi[n+1] - psi[Nx];
      }else if(n==(Nx)){
      locfirst = psi[1] - psi[n-1];
      }else{
      locfirst = psi[n+1] - psi[n-1];
      }



      ef_xg[n] = -locfirst/(2*h_x);
  
      printf("%lf \t%d %lf\n", ef_xg[n], n, psi[n]);
  
    }


 printf("done E_field");

}
