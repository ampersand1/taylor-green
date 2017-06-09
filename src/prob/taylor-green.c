#include "copyright.h"
/*============================================================================*/
/*! \file taylor-green.c
 *  \brief Problem generator for taylor-green vortex problem.
 *
 * REFERENCE: Paradigmatic flow for small-scale magnetohydrodynamics
 * E. Lee et. al				      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include <string.h>


#ifndef MHD
#error : The problem generator orszag_tang.c only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,is,ie,js,je,ks,ke,nx1,nx2,nx3;
  Real L, p, p0, b0,v0,d0, ***ax, ***ay, ***az, x1,x2,x3,gamma;

  /*  grid setting  */
  is = pGrid->is;
  ie = pGrid->ie;

  js = pGrid->js;
  je = pGrid->je;

  ks = pGrid->ks;
  ke = pGrid->ke;

  /*Based on orszag-tang.c file*/
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

    //This is a way of declaring certain variables.
    //When I put this in the "illegal instruction 4" error went away.
    if ((ax = (Real***)calloc_3d_array(nx3,nx2, nx1, sizeof(Real))) == NULL) {
        ath_error("[taylor-green]: Error allocating memory for vector pot\n");
    }

    if ((ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
        ath_error("[taylor-green]: Error allocating memory for vector pot\n");
    }

    if ((az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
        ath_error("[taylor-green]: Error allocating memory for vector pot\n");
    }

    if (pGrid->Nx[2] == 1) {
       ath_error("[taylor-green]: This problem can only be run in 3D\n");
     }

  /* Variables Set in Input file */

  b0              = par_getd("problem","b0");
  v0              = par_getd("problem","v0");
  d0              = par_getd("problem","d0");
  p0              = par_getd("problem","p0");
  L               = par_getd("problem","L");

 /*internal energy density -> pressure in this case */
/*********Functional Form of Pressure*************************/
//p= p0+(d0*(v0*v0)/16.0)*(cos(2*x1/L)+cos(2*x2/L))*(cos(2*x3/L)+2);
//p = (1.0/gamma - ((v0*v0)/16.0)*(cos(2*x1/L)+cos(2*x2/L))*(cos(2*x3/L)+2)) ;
/******************************************************/

  /* Initialize vector potential */

for(k=ks; k<=ke+1; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      x1 -= 0.5*pGrid->dx1;
      x2 -= 0.5*pGrid->dx2;
      x3 -= 0.5*pGrid->dx3;

      /*Magnetic potential currently set to insulating*/
      ax[k][j][i] = -b0*sin(x1)*cos(x2)*cos(x3) ;
      ay[k][j][i] =  b0*cos(x1)*sin(x2)*cos(x3) ;
      az[k][j][i] =  0 ;

    }
  }
}
  /******************************************************/
  /* Initialize density, momentum, face-centered fields */
for  (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      pGrid->U[k][j][i].d = d0;
      pGrid->U[k][j][i].M1 = d0*v0*sin(x1)*cos(x2)*cos(x3);
      pGrid->U[k][j][i].M2 = -d0*v0*cos(x1)*sin(x2)*cos(x3);
      pGrid->U[k][j][i].M3 = 0.0;
      //pGrid->U[k][j][i].E = 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
              //+ SQR(pGrid->U[k][j][i].B3c)) + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
              //+ SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
      //pGrid->U[k][j][i].E = (1.0/gamma - ((v0*v0)/16.0)*(cos(2*x1/L)+cos(2*x2/L))*(cos(2*x3/L)+2))/(gamma-1.0) ;

      /*Curl of A*/
      pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 - (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
      pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 - (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
      pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 - (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
    }
  }
}

/******************************************************/
/******************************************************/

/* boundary conditions on interface B */
for  (k=ks; k<=ke; k++) {
 for (j=js; j<=je; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 - (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
  }
 }
}

for (k=ks; k<=ke; k++){
 for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie; i++)   {
    pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 - (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
  }
 }
}
for (k=ks; k<=ke+1; k++){
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++)  {
    pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 - (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
  }
 }
}

  /******************************************************/
/* initialize total energy and cell-centered B */
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
    pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
    pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);

//p/Gamma_1

#ifndef ISOTHERMAL
    pGrid->U[k][j][i].E = p0/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
             + SQR(pGrid->U[k][j][i].B3c))
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
              + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif

  }
 }
}
  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/
/*
static Real current1(const GridS *pG, const int i, const int j, const int k)
 {  return ((pG->B3i[k][j][i]-pG->B3i[k][j-1][i])/pG->dx2 -
            (pG->B2i[k][j][i]-pG->B2i[k-1][j][i])/pG->dx3);
 }


 static Real current2(const GridS *pG, const int i, const int j, const int k)
 {  return ((pG->B1i[k][j][i]-pG->B1i[k][j-1][i])/pG->dx2 -
            (pG->B2i[k][j][i]-pG->B2i[k][j][i-1])/pG->dx1);
 }


 static Real current3(const GridS *pG, const int i, const int j, const int k)
 {  return ((pG->B2i[k][j][i]-pG->B2i[k][j][i-1])/pG->dx1 -
            (pG->B1i[k][j][i]-pG->B1i[k][j-1][i])/pG->dx2);
 }
*/
void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;

}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
//  if(strcmp(expr,"J1")==0) return current1;
//  if(strcmp(expr,"J2")==0) return current2;
//  if(strcmp(expr,"J3")==0) return current3;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
