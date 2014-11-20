//-----------------------------------------------------------------------------
//Author:        C. Letertre
//               B. Schorr
//               Sven Schumann (schumans@kph.uni-mainz.de)
//Description:   TMultiDim: Multi-dimensional interpolation class
//               KERNLIB
//-----------------------------------------------------------------------------
// TMultiDim.cc
//
// Multidimensional linear interpolation
// Taken from the KERNLIB function E104, converted from Fortran to C,
// and packaged in ROOT class TMultiDim
//-----------------------------------------------------------------------------

#include "TMultiDim.h"

ClassImp(TMultiDim)

//-----------------------------------------------------------------------------

TMultiDim::TMultiDim() : TObject()
{
  //Initialise array pointers
  fValue = NULL;
  fCoord = NULL;
  fNBins = NULL;
  fBin   = NULL;
  fPoint = NULL;
  //Default c'tor starts with no dimensions at all
  fNDim = 0;
}

//-----------------------------------------------------------------------------

TMultiDim::TMultiDim(UInt_t ndim) : TObject()
{
  //Initialise array pointers
  fValue = NULL;
  fCoord = NULL;
  fNBins = NULL;
  fBin   = NULL;
  fPoint = NULL;
  //Apply given number of dimensions
  SetNDim(ndim);
}

//-----------------------------------------------------------------------------

TMultiDim::~TMultiDim()
{
  //Clean up arrays
  delete[] fValue;
  delete[] fCoord;
  delete[] fNBins;
  delete[] fBin;
  delete[] fPoint;
}

//-----------------------------------------------------------------------------

void TMultiDim::AllocMem()
{
  //Determine total number of data values (N = m_n * m_n-1 * ... * m_1 * m_0)
  fNValue = fNBins[0];
  for(UInt_t i=1; i<fNDim; i++)
    fNValue*=fNBins[i];
  //Clear existing data array (if existing) and allocate new
  if(fValue) delete[] fValue;
  fValue = new Double_t[fNValue];

  //Determine total number of coordinate values (N = p_n + p_n-1 + ... + p_1 + p_0)
  fNCoord = fNBins[0];
  for(UInt_t i=1; i<fNDim; i++)
    fNCoord+=fNBins[i];
  //Clear existing data array (if existing) and allocate new
  if(fCoord) delete[] fCoord;
  fCoord = new Double_t[fNCoord];
}

//-----------------------------------------------------------------------------

UInt_t TMultiDim::GetValueIdx()
{
  UInt_t Idx; //Final array index
  UInt_t Bin; //Bin position in currently processed dimension

  //Consider an n-dimensional C array with maximum ranges m_n and indices i_n,
  //i.e. defined as
  //  Arr[m_n][m_n-1]...[m_0]
  //and accessed as
  //  Arr[i_n][i_n-1]...[i_0]
  //with i_n within [0, m_n - 1]
  //This array can be flattened to a 1-dimensional C array
  //  Arr[m_n * m_n-1 * ... * m_0]
  //For a given set of i_n indices, the 1-dimensional index j is then given by
  //  j = i_n * (m_n-1 * m_n-2 * ... * m_0) + i_n-1 * (m_n-2 * m_n-3 * ... * m_0)
  //    + ... + i_1 * m_0 + i_0

  //The following loops calculate the equivalent 1-dim index j for n-dim. indices i_n within ranges m_n
  //with the following translations:
  //  Idx = j
  //  fBin[n] = i_n
  //  fNBins[n] = m_n

  //The position i_0 in the lowest dimension n = 0 can be taken directly.
  Idx = fBin[0];
  //Then, process all dimensions n > 0
  for(UInt_t t=fNDim-1; t>0; t--)
  {
    Bin = fBin[t];      //Select the position i_n within the current dimension, ...
    for(UInt_t s=t; s>0; s--)
      Bin*=fNBins[s-1]; //...multiply with the sizes of all lower dimensions, ...
    Idx+=Bin;           //...and add to flattened index.
  }

  return Idx;
}

//-----------------------------------------------------------------------------

UInt_t TMultiDim::GetCoordIdx(UInt_t dim, UInt_t bin)
{
  UInt_t Idx; //Final array index

  Idx = 0;
  //Skip all bins in previous dimensions...
  for(UInt_t t=0; t<dim; t++)
    Idx+=fNBins[t];
  //...and select chosen bin
  Idx+=bin;

  return Idx;
}

//-----------------------------------------------------------------------------

void TMultiDim::SetNDim(UInt_t ndim)
{
  if((ndim>=1) && (ndim<=10)) //Check for valid dimension ranges
  {
    fNDim = ndim;

    //Clear (if existing) and allocate array (=vector) for number of bins in each dimension
    if(fNBins) delete[] fNBins;
    fNBins = new UInt_t[fNDim];

    //Clear (if existing) and allocate array (=vector) for current bin in each dimension
    if(fBin) delete[] fBin;
    fBin = new UInt_t[fNDim];

    //Clear (if existing) and allocate array (=vector) for coordinate position in each dimension
    if(fPoint) delete[] fPoint;
    fPoint = new Double_t[fNDim];
  }
}

//-----------------------------------------------------------------------------

void TMultiDim::SetCoord(UInt_t dim, UInt_t bin, Double_t coord)
{
  //Check for valid dimension and bin ranges
  if(dim>=fNDim) return;
  if(bin>=fNBins[dim]) return;

  fCoord[GetCoordIdx(dim, bin)] = coord;
}

//-----------------------------------------------------------------------------

void TMultiDim::SetCoord(UInt_t dim, Double_t* coord)
{
  //Check for valid dimension range
  if(dim>=fNDim) return;

  //Store all coordinate values for chosen dimension
  for(UInt_t t=0; t<fNBins[dim]; t++)
    fCoord[GetCoordIdx(dim, t)] = coord[t];
}

//-----------------------------------------------------------------------------

Double_t TMultiDim::GetCoord(UInt_t dim, UInt_t bin)
{
  //Check for valid dimension and bin ranges
  if(dim>=fNDim) return 0.0;
  if(bin>=fNBins[dim]) return 0.0;

  return fCoord[GetCoordIdx(dim, bin)];
}

//-----------------------------------------------------------------------------
