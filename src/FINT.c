//-----------------------------------------------------------------------------
//Author:        C. Letertre
//               B. Schorr
//               Sven Schumann (schumans@kph.uni-mainz.de)
//Description:   TMultiDim: Multi-dimensional interpolation class
//               KERNLIB
//-----------------------------------------------------------------------------
// FINT.c
//
// Multidimensional linear interpolation
// Taken from the KERNLIB function E104, converted from Fortran to C,
// and packaged in ROOT class TMultiDim
//-----------------------------------------------------------------------------

#include <stdio.h>

//-----------------------------------------------------------------------------

//Select minimum & maximum values from 2 given integers
inline int min(int const A, int const B){ return A < B ? A : B; }
inline int max(int const A, int const B){ return A > B ? A : B; }

//-----------------------------------------------------------------------------

//Fortran-to-C ported multidimensional interpolation, ripped and cleaned up
//from CERNLib Fortran function FINT (E104)
double FINT(int* NARG, double* ARG, int* NENT, double* ENT, double* TABLE)
{
  double WEIGHT[1024];
  double X, H, ETA;
  double tmp, fint;
  int INDEX[1024];
  int LMIN, LMAX;
  int ISTEP, ISHIFT;
  int KNOTS;
  int NDIM;
  int LOCA, LOCB, LOCC;
  int N, K, I;

  //Adjust for different array start indices in Fortran (1) and C (0)
  --TABLE;
  --ENT;
  --NENT;
  --ARG;

  fint = 0.0;
  if(((*NARG) < 1) || ((*NARG) > 10)) goto l300;

  LMAX = 0;
  ISTEP = 1;
  KNOTS = 1;
  INDEX[0] = 1;
  WEIGHT[0] = 1.0;
  for(N=1; N<(*NARG)+1; N++)
  {
    X = ARG[N];
    NDIM = NENT[N];
    LOCA = LMAX;
    LMIN = LMAX + 1;
    LMAX+=NDIM;
    if(NDIM > 2) goto l010;
    if(NDIM==1)  goto l100;
    H = X - ENT[LMIN];
    if(H==0.0) goto l090;
    ISHIFT = ISTEP;
    if((X-ENT[LMIN+1])==0.0) goto l021;
    ISHIFT = 0;
    ETA = H / (ENT[LMIN+1] - ENT[LMIN]);
    goto l030;
l010:
    LOCB = LMAX + 1;
l011:
    LOCC = (LOCA+LOCB) / 2;
    tmp = X - ENT[LOCC];
    if(tmp < 0.0) goto l012;
    if(tmp > 0.0) goto l013;
    goto l020;
l012:
    LOCB = LOCC;
    goto l014;
l013:
    LOCA = LOCC;
l014:
    if((LOCB-LOCA) > 1) goto l011;
    LOCA =  min(max(LOCA, LMIN), LMAX-1);
    ISHIFT = (LOCA - LMIN) * ISTEP;
    ETA = (X - ENT[LOCA]) / (ENT[LOCA+1] - ENT[LOCA]);
    goto l030;
l020:
    ISHIFT = (LOCC - LMIN) * ISTEP;
l021:
    for(K=1; K<KNOTS+1; K++)
    {
      INDEX[K-1]+=ISHIFT;
    }
    goto l090;
l030:
    for(K=1; K<KNOTS+1; K++)
    {
      INDEX[K-1]+=ISHIFT;
      INDEX[K-1+KNOTS] = INDEX[K-1] + ISTEP;
      WEIGHT[K-1+KNOTS] = WEIGHT[K-1] * ETA;
      WEIGHT[K-1]-=WEIGHT[K-1+KNOTS];
    }
    KNOTS*=2;
l090:
    ISTEP*=NDIM;
l100:
    ; //NOP statement
  }
  for(K=1; K<KNOTS+1; K++)
  {
    I = INDEX[K-1];
    fint+=(WEIGHT[K-1] * TABLE[I]);
  }
  return fint;

l300:
  printf("       FUNCTION FINT ... NARG =%6d NOT WITHIN RANGE\n", (*NARG));
  return fint;
}

//-----------------------------------------------------------------------------
