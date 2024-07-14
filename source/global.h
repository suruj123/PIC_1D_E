/*-----------------------------------------------------------------------------

    MPMD-v2.0 : MULTI POTENTIAL MOLECULAR DYNAMICS-version 2.0 
    A parallel classical molecular dynamics code
    Copyright (C) 2018  Harish Charan, charan.harish@gmail.com

    This program is free software but a proper permission must be taken: you can 
    redistribute it and/or modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
    more details.

    You should have received a copy of the GNU General Public License along 
    with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    See the README file in the top-level MPMD-v2.0 directory.

-----------------------------------------------------------------------------*/



typedef double real;

#include <fftw3.h>

#define	NDIM 3
#define Sqr(x) ((x) * (x))
#define SignR(x, y) (((y) >= 0) ? (x) : (- (x)))

double	*rx, *ry, *rz, *vx, *vy, *vz, *ax, *ay, *az, *dx1 ,*p,*fax, *fay, *faz, *alpha, *Qold, **a, *b;

double *Q;

real    *A, *Ax, *Ay;
real    *TValSum;

real	region[NDIM+1], regionH[NDIM+1], GAMMA, temperature, 
	deltaT, density, kinEnergy, potEnergy, potEnergy1, pressure, rCut, kappa,
	sKinEnergy, sPressure, sTotEnergy, ssKinEnergy, ssPressure, 
	ssTotEnergy, timeNow, totEnergy, totEnergy1, uSum, virSum, svirSum,
	sPotEnergy, ssPotEnergy, potEnergy,row, sTotEnergy1, ssTotEnergy1,
	vSum, vvSum, keConfig, sKEConfig, ssKEConfig, sPotEnergy1, ssPotEnergy1, PE, KE,
        aKE, aPE;


int     npoints, taskid, k, numtasks, first, actual,remnant;
int     source, start, buffer[2];
int     dest, offset, npoint;

char    BC, thermo, cellAlgo;

int	*atomType;


int	initUcell[NDIM+1], moreCycles, stepAvg, 
	stepCount, stepEquil, stepLimit, stepTrajectory, stepDump;

int	*cellList, cells[NDIM+1],nAtom;

int 	rank, size, master;
real	*fax, *fay, fuSum, fvirSum;

// Spacetime correlations for k = k_x
real	**spacetimeFluct,**complxFuncOrg, *complxFunc,*spacetimeFluctAv;
int	*indexFluct, countFluctAv;
int	nFunCorr;

// Spacetime Correlations for k = k_y
real	**spacetimeFluct2,**complxFuncOrg2, *complxFunc2,*spacetimeFluctAv2;
int	*indexFluct2, countFluctAv2;
int      limitCorrAv, nBuffCorr, nFunCorr2, nValCorr, stepCorr;

// Viscosity
real  	rfAtom, frfAtom;
real  	 **viscAcf, *viscAcfOrg, *viscAcfAv, viscAcfInt;
int     nValAcf, nBuffAcf, stepAcf, countAcfAv, limitAcfAv, *indexAcf;

// Radial distribution function
real 	**histRdf,*HistRdf, rangeRdf, rBin;
int	countRdf, limitRdf, sizeHistRdf, stepRdf;

// Diffusion by Einstein's defination 
real	**rDiffuseTrue, **rDiffuseOrg, **xDiffuse, **yDiffuse, **rDiffuse, *xDiffuseAv,*yDiffuseAv,*rDiffuseAv,*MSDAv;
int	*indexDiffuse, countDiffuseAv ;

// Velocity auto correlation function
real	**VeloAcfOrg, **VeloAcf, *VeloAcfAv, diffuseVeloAcfInt;
int	*indexVeloAcf, countVeloAcfAv ; 

// Velocity distribution
real *histVel, rangeVel;
int countVel, limitVel, sizeHistVel, stepVel;

//Charge
double lambda_c, norm;

//Interpolation
double del_x, del_y, del_z, h_volume, h_area;

double  *xs, *ys, *zs, *xg, *yg, *zg;//,

double  *qg, *dqg_dx, *dqg_dy, *dqg_dz;//, *rho;

fftw_complex *rho_k, *phi_k;

double  *rho, *phi;

double *ef_x, *ef_y;

double *dq_dx, *dq_dy, *dq_dz;

int Nx, Ny, Nz;

//Poisson Solver

double *psi, *psinew;

int ngrids_X, ngrids_Y;


//E_field

double *ef_xg, *ef_yg;


// Output files prefixes
char    *prefix;
char    dirprefix[250];

char    result[250];
FILE    *fpresult;

char 	xyz[256];
FILE 	*fpxyz;

char 	charge[256];
FILE 	*fpcharge;

//char 	interpolation[256];
//FILE 	*fpinterpolation;

char    rdf[256];
FILE    *fprdf;
