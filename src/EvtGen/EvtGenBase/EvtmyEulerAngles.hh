//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of models developed at BES collaboration
//      based on the EvtGen framework.  If you use all or part bn 
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/BesCopyright
//      Copyright (A) 2006      Ping Rong-Gang @IHEP
//
// Module:  EvtDIY.cc
//
// Description:  Class to calculate the Euler angles to rotate a system
//
// Modification history:
//
//    Ping R.-G.       December, 2007       Module created
//
//------------------------------------------------------------------------
//
#ifndef EVTMYEULERANGLES_HH
#define EVTMYEULERANGLES_HH
  
#include <iostream>
// #include <math.h>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector3R.hh"
  
//using std::endl;   
//using std::fstream; 

class EvtmyEulerAngles
{ 

public: 

//destructor
virtual ~EvtmyEulerAngles();

// constructor
EvtmyEulerAngles();
EvtmyEulerAngles( const EvtVector3R & Yaxis, const EvtVector3R & Zaxis);
EvtmyEulerAngles( const EvtVector4R & Pyaxis, const EvtVector4R & Pzaxis);
EvtmyEulerAngles( EvtVector4R &h1,  EvtVector4R &h2);// Euler rotation from h1 to h2
  EvtmyEulerAngles( EvtVector4R x1,  EvtVector4R y1,EvtVector4R x2,  EvtVector4R y2, std::string axisTag);// Euler rotation from (x1,y1,z1) to (x2,y2,z2)
  void EulerAngles();
  double getAlpha();
  double getBeta();
  double getGamma();
  double getphi(); //alpha  using EvtmyEulerAngles( EvtVector4R h1,  EvtVector4R h2);
  double gettheta();//beta
  double getpsi();// gamma
private:
  EvtVector3R _Yaxis,_Zaxis;
  double _alpha,_beta,_gamma,_ry,_rz;
  EvtVector4R _x, _X, xx,yy,zz, XX,YY,ZZ;
  EvtVector4R Zcm,_N;
};

#endif
