//-----------------------------------------------------------------------------
//Author:        C. Letertre
//               B. Schorr
//               Sven Schumann (schumans@kph.uni-mainz.de)
//Description:   TMultiDim: Multi-dimensional interpolation class
//               KERNLIB
//-----------------------------------------------------------------------------
// TMultiDim.h
//
// Multidimensional linear interpolation
// Taken from the KERNLIB function E104, converted from Fortran to C,
// and packaged in ROOT class TMultiDim
//-----------------------------------------------------------------------------

#ifndef __TMultiDim__
#define __TMultiDim__

#include "TObject.h"

//---------------------------------------------------------------------------

extern "C"
{
  //Multidim. linear interpolation ripped and ported from kernlib E104
  Double_t FINT(UInt_t&, Double_t*, UInt_t*, Double_t*, Double_t*);
}
//---------------------------------------------------------------------------

class TMultiDim : public TObject
{
 private:
  //Data members
  UInt_t    fNDim;   //Number of dimensions of parameter space
  UInt_t    fNValue; //Total number of data points in n-dim. parameter space
  UInt_t    fNCoord; //Total number of coordinate positions of n-dim. parameter space
  UInt_t*   fNBins;  //[fNDim]:   Array (=vector) for number of bins in each dimension
  UInt_t*   fBin;    //[fNDim]:   Array (=vector) for current bin positions in each dimension
  Double_t* fValue;  //[fNValue]: (Flattened) array of data values for each point in n-dim. parameter space
  Double_t* fCoord;  //[fNCoord]: Array for coordinate positions of n-dim. parameter space
  Double_t* fPoint;  //[fNDim]:   Array (=vector) for point in n-dim. parameter space for which to interpolate
  //Internal methods
  UInt_t    GetValueIdx();                       //Calculates flattened 1-dim. array index from n-dim. bin position vector
  UInt_t    GetCoordIdx(UInt_t dim, UInt_t bin); //Calculates index in coordinate array for given dimension and bin position
  void      AllocMem();                          //Reserves memory for current settings of dimensions and binning 
 public:
  //C'tor/d'tor
  TMultiDim();            //Default c'tor
  TMultiDim(UInt_t ndim); //C'tor with dimension setting
  virtual ~TMultiDim();   //Default d'tor
  //Set/get number of dimensions
  virtual UInt_t   GetNDim() { return fNDim; }
  virtual void     SetNDim(UInt_t ndim);
  //Set/get number of bins in each dimension. Exists with n-dim. bin vector and overloaded with explicit bin positions for 1 to 10 dimensions
  virtual void     SetNBins(UInt_t* n);
  virtual void     SetNBins(UInt_t n0);
  virtual void     SetNBins(UInt_t n0, UInt_t n1);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8);
  virtual void     SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9);
  virtual UInt_t   GetNBins(UInt_t dim) { if(fNBins && (dim < fNDim)) return fNBins[dim]; else return 0;}
  //Set/get data value at given bin. Exists with n-dim. bin vector and overloaded with explicit bin positions for 1 to 10 dimensions
  virtual void     SetValue(UInt_t* n, Double_t value);
  virtual void     SetValue(UInt_t n0, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, Double_t value);
  virtual void     SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9, Double_t value);
  virtual Double_t GetValue(UInt_t* n);
  virtual Double_t GetValue(UInt_t n0);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8);
  virtual Double_t GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9);
  //Set/get coordinate positions for given dimension and bin position
  virtual void     SetCoord(UInt_t dim, UInt_t bin, Double_t coord);
  virtual void     SetCoord(UInt_t dim, Double_t* coord);
  virtual Double_t GetCoord(UInt_t dim, UInt_t bin);
  //Perform interpolation at given point. Exists with n-dim. point vector and overloaded with explicit point coordinates for 1 to 10 dimensions
  virtual Double_t Interpolate(Double_t* x);
  virtual Double_t Interpolate(Double_t x0);
  virtual Double_t Interpolate(Double_t x0, Double_t x1);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7, Double_t x8);
  virtual Double_t Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7, Double_t x8, Double_t x9);

  ClassDef(TMultiDim, 1)
};

//---------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t* n)
{
  for(UInt_t t=0; t<fNDim; t++)
  {
    fBin[t] = n[t]; if(n[t]>=fNBins[t]) return 0.0; //Cancel if invalid bin position is requested
  }

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0)
{
  if(fNDim!=1) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1)
{
  if(fNDim!=2) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2)
{
  if(fNDim!=3) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n1>=fNBins[1]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3)
{
  if(fNDim!=4) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n3>=fNBins[3]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4)
{
  if(fNDim!=5) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n3>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n4>=fNBins[4]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5)
{
  if(fNDim!=6) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n3>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n4>=fNBins[4]) return 0.0;
  fBin[5] = n5; if(n5>=fNBins[5]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6)
{
  if(fNDim!=7) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n3>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n4>=fNBins[4]) return 0.0;
  fBin[5] = n5; if(n5>=fNBins[5]) return 0.0;
  fBin[6] = n6; if(n6>=fNBins[6]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7)
{
  if(fNDim!=8) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n2>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n3>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n4>=fNBins[4]) return 0.0;
  fBin[5] = n5; if(n5>=fNBins[5]) return 0.0;
  fBin[6] = n6; if(n6>=fNBins[6]) return 0.0;
  fBin[7] = n7; if(n7>=fNBins[7]) return 0.0;

  return fValue[GetValueIdx()]; //Check if called function is compatible with number of dimensions
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8)
{
  if(fNDim!=9) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n1>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n1>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n1>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n1>=fNBins[4]) return 0.0;
  fBin[5] = n5; if(n1>=fNBins[5]) return 0.0;
  fBin[6] = n6; if(n1>=fNBins[6]) return 0.0;
  fBin[7] = n7; if(n1>=fNBins[7]) return 0.0;
  fBin[8] = n8; if(n1>=fNBins[8]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::GetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9)
{
  if(fNDim!=10) return 0.0; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n1>=fNBins[0]) return 0.0; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return 0.0;
  fBin[2] = n2; if(n1>=fNBins[2]) return 0.0;
  fBin[3] = n3; if(n1>=fNBins[3]) return 0.0;
  fBin[4] = n4; if(n1>=fNBins[4]) return 0.0;
  fBin[5] = n5; if(n1>=fNBins[5]) return 0.0;
  fBin[6] = n6; if(n1>=fNBins[6]) return 0.0;
  fBin[7] = n7; if(n1>=fNBins[7]) return 0.0;
  fBin[8] = n8; if(n1>=fNBins[8]) return 0.0;
  fBin[9] = n9; if(n1>=fNBins[9]) return 0.0;

  return fValue[GetValueIdx()];
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t* n, Double_t value)
{
  for(UInt_t t=0; t<fNDim; t++) //Check if called function is compatible with number of dimensions
  {
    fBin[t] = n[t]; if(n[t]>=fNBins[t]) return; //Cancel if invalid bin position is requested
  }

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------


inline void TMultiDim::SetValue(UInt_t n0, Double_t value)
{
  if(fNDim!=1) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, Double_t value)
{
  if(fNDim!=2) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, Double_t value)
{
  if(fNDim!=3) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, Double_t value)
{
  if(fNDim!=4) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, Double_t value)
{
  if(fNDim!=5) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, Double_t value)
{
  if(fNDim!=6) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;
  fBin[5] = n5; if(n5>=fNBins[5]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, Double_t value)
{
  if(fNDim!=7) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;
  fBin[5] = n5; if(n5>=fNBins[5]) return;
  fBin[6] = n6; if(n6>=fNBins[6]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, Double_t value)
{
  if(fNDim!=8) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;
  fBin[5] = n5; if(n5>=fNBins[5]) return;
  fBin[6] = n6; if(n6>=fNBins[6]) return;
  fBin[7] = n7; if(n7>=fNBins[7]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, Double_t value)
{
  if(fNDim!=9) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;
  fBin[5] = n5; if(n5>=fNBins[5]) return;
  fBin[6] = n6; if(n6>=fNBins[6]) return;
  fBin[7] = n7; if(n7>=fNBins[7]) return;
  fBin[8] = n8; if(n8>=fNBins[8]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetValue(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9, Double_t value)
{
  if(fNDim!=10) return; //Check if called function is compatible with number of dimensions

  fBin[0] = n0; if(n0>=fNBins[0]) return; //Cancel if invalid bin position is requested
  fBin[1] = n1; if(n1>=fNBins[1]) return;
  fBin[2] = n2; if(n2>=fNBins[2]) return;
  fBin[3] = n3; if(n3>=fNBins[3]) return;
  fBin[4] = n4; if(n4>=fNBins[4]) return;
  fBin[5] = n5; if(n5>=fNBins[5]) return;
  fBin[6] = n6; if(n6>=fNBins[6]) return;
  fBin[7] = n7; if(n7>=fNBins[7]) return;
  fBin[8] = n8; if(n8>=fNBins[8]) return;
  fBin[9] = n9; if(n9>=fNBins[9]) return;

  fValue[GetValueIdx()] = value;
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t* n)
{
  for(UInt_t t=0; t<fNDim; t++)
    fNBins[t] = n[t];

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0)
{
  if(fNDim!=1) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1)
{
  if(fNDim!=2) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2)
{
  if(fNDim!=3) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3)
{
  if(fNDim!=4) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4)
{
  if(fNDim!=5) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5)
{
  if(fNDim!=6) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;
  fNBins[5] = n5;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6)
{
  if(fNDim!=7) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;
  fNBins[5] = n5;
  fNBins[6] = n6;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7)
{
  if(fNDim!=8) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;
  fNBins[5] = n5;
  fNBins[6] = n6;
  fNBins[7] = n7;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8)
{
  if(fNDim!=9) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;
  fNBins[5] = n5;
  fNBins[6] = n6;
  fNBins[7] = n7;
  fNBins[8] = n8;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline void TMultiDim::SetNBins(UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6, UInt_t n7, UInt_t n8, UInt_t n9)
{
  if(fNDim!=10) return; //Check if called function is compatible with number of dimensions

  fNBins[0] = n0;
  fNBins[1] = n1;
  fNBins[2] = n2;
  fNBins[3] = n3;
  fNBins[4] = n4;
  fNBins[5] = n5;
  fNBins[6] = n6;
  fNBins[7] = n7;
  fNBins[8] = n8;
  fNBins[9] = n9;

  AllocMem();
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0)
{
  if(fNDim!=1) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1)
{
  if(fNDim!=2) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2)
{
  if(fNDim!=3) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3)
{
  if(fNDim!=4) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4)
{
  if(fNDim!=5) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5)
{
  if(fNDim!=6) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;
  fPoint[5] = x5;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6)
{
  if(fNDim!=7) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;
  fPoint[5] = x5;
  fPoint[6] = x6;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7)
{
  if(fNDim!=8) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;
  fPoint[5] = x5;
  fPoint[6] = x6;
  fPoint[7] = x7;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7, Double_t x8)
{
  if(fNDim!=9) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;
  fPoint[5] = x5;
  fPoint[6] = x6;
  fPoint[7] = x7;
  fPoint[8] = x8;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, Double_t x7, Double_t x8, Double_t x9)
{
  if(fNDim!=10) return 0.0; //Check if called function is compatible with number of dimensions

  fPoint[0] = x0;
  fPoint[1] = x1;
  fPoint[2] = x2;
  fPoint[3] = x3;
  fPoint[4] = x4;
  fPoint[5] = x5;
  fPoint[6] = x6;
  fPoint[7] = x7;
  fPoint[8] = x8;
  fPoint[9] = x9;

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

inline Double_t TMultiDim::Interpolate(Double_t* x)
{
  for(UInt_t t=0; t<fNDim; t++)
    fPoint[t] = x[t];

  return FINT(fNDim, fPoint, fNBins, fCoord, fValue);
}

//-----------------------------------------------------------------------------

#endif
