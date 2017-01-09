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

#ifndef MHD
#error : The problem generator orszag_tang.c only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is,ie,js,je,ks,ke,i,j,k,nx1,nx2,nx3;
  Real b0,v0,d0, **ax, **ay, **az;
  
  /*  is js ks seems like start for the 3 directions i,j,k   */
  is = pGrid->is;
  ie = pGrid->ie;

  js = pGrid->js;
  je = pGrid->je;

  ks = pGrid->ks;
  ke = pGrid->ke;

  /*Not exactly sure what this is for but it was in the orszag-tang.c file*/
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  /* Temporarily set B0 and v0 */
  b0 = 1.0;
  v0 = 1.0;
  d0 = 1.0;

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
      ax[k][j] = -b0*sin(x1)*cos(x2)*cos(x3) ; 
      ay[k][i] =  b0*cos(x1)*sin(x2)*cos(x3) ;
      az[j][i] =  0 

    }
  }
}

  /******************************************************/
  /* Initialize density, momentum, face-centered fields */
for  (k=ks; k<=ke; j++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      pGrid->U[k][j][i].d = d0;
      pGrid->U[k][j][i].M1 = d0*v0*sin(x1)*cos(x2)*cos(x3);
      pGrid->U[k][j][i].M2 = -d0*v0*cos(x1)sin(x2)*cos(x3);
      pGrid->U[k][j][i].M3 = 0.0;

      /*Curl of A*/
      pGrid->B1i[k][j][i] = (az[j+1][i] - az[j][i])/pGrid->dx2 - (ay[k+1][i] - ay[k][i])/pGrid->dx3;
      pGrid->B2i[k][j][i] = (ax[k+1][j] - ax[k][j])/pGrid->dx3 - (az[j][i+1] - az[j][i])/pGrid->dx1;
      pGrid->B3i[k][j][i] = (ay[k][i+1] - ay[k][i])/pGrid->dx1 - (ax[k][j+1] - ax[k][j])/pGrid->dx2;
    }
  }
}

/******************************************************/
/* NOTE TO SELF!!!!
Look through paper for what boundary condition for B is supposed to be.  
Also intialize the total energy and cell-centered B
Consider using Orszag-tang.c values of b0, v0. 

Also what is p0?  
What should the total energy be?
*/
  /******************************************************/

/* boundary conditions on interface B */
for  (k=ks; k<=ke; j++) {  
 for (j=js; j<=je; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[ks][j][i] = (az[j+1][i] - az[j][i])/pGrid->dx2;
  }
 }
  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->B2i[ks][j][i] =-(az[j][i+1] - az[j][i])/pGrid->dx1;
  }
 }
}

  /******************************************************/
/* initialize total energy and cell-centered B */

  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->U[ks][j][i].B1c = 0.5*(pGrid->B1i[ks][j][i]+pGrid->B1i[ks][j][i+1]);
    pGrid->U[ks][j][i].B2c = 0.5*(pGrid->B2i[ks][j][i]+pGrid->B2i[ks][j+1][i]);
    pGrid->U[ks][j][i].B3c = 0.0;
#ifndef ISOTHERMAL
    pGrid->U[ks][j][i].E = p0/Gamma_1
        + 0.5*(SQR(pGrid->U[ks][j][i].B1c) + SQR(pGrid->U[ks][j][i].B2c)
             + SQR(pGrid->U[ks][j][i].B3c))
        + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)
              + SQR(pGrid->U[ks][j][i].M3))/pGrid->U[ks][j][i].d;
#endif
  }}




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
