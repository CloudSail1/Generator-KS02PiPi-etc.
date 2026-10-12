//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of models developed at BES collaboration
//      based on the EvtGen framework.  If you use all or part
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
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <fstream>   
#include <stdio.h>
#include <cstdlib>
#include "EvtGenBase/EvtmyEulerAngles.hh"
#include "EvtGenBase/EvtReport.hh"
using namespace std;  //::endl;

EvtmyEulerAngles::~EvtmyEulerAngles(){}
EvtmyEulerAngles::EvtmyEulerAngles(const EvtVector3R & Yaxis, const EvtVector3R & Zaxis){
  _Yaxis=Yaxis;
  _Zaxis=Zaxis;
   EulerAngles();
}

EvtmyEulerAngles::EvtmyEulerAngles( EvtVector4R &h1, EvtVector4R & h2){
  EvtVector4R z0(0,0,0,1);
  Zcm=z0;
  _x=h1;
  _X=h2;
  yy=Zcm.cross(_x);
  zz=h1;
  xx=yy.cross(zz);

  YY=Zcm.cross(_X);
  ZZ=_X;
  XX=YY.cross(ZZ);
  _N=zz.cross(ZZ);
  
}

EvtmyEulerAngles::EvtmyEulerAngles( EvtVector4R x1, EvtVector4R  y1, EvtVector4R x2, EvtVector4R  y2, std::string axisTag){
  
  EvtVector4R z0(0,0,0,1);
  if(axisTag=="xy"){
    xx = x1;
    yy = y1;
    zz=xx.cross(yy);
    
    XX = x2;
    YY = y2;
    ZZ=XX.cross(YY);
  }else if(axisTag=="yz"){
    yy=x1;
    zz=y1;
    xx=yy.cross(zz);

    YY=x2;
    ZZ=y2;
    XX=YY.cross(ZZ);  
  }else{
    std::cout<<"EvtmyEulerAngles::EvtmyEulerAngles bad axisTag option"<<std::endl;
    abort();
  }

  _N=zz.cross(ZZ);
  //std::cout<<" _N "<<zz<<ZZ<<_N<<yy<<std::endl;
  
}

double EvtmyEulerAngles::getphi(){
  double cp=yy.dot(_N)/yy.d3mag()/_N.d3mag();
  if(fabs(cp-1)<0.00001) return 0;
  double acoscp = acos(cp);
  double xyz=zz.dot(yy.cross(_N));
  if(xyz>=0){return acoscp;}else{ return 2*3.1415926-acoscp;}

}

double EvtmyEulerAngles::gettheta(){
  double cp=zz.dot(ZZ)/zz.d3mag()/ZZ.d3mag();
  return acos(cp);
}


double EvtmyEulerAngles::getpsi(){
  double cp=YY.dot(_N)/YY.d3mag()/_N.d3mag();
  if(fabs(cp-1)<0.00001) return 0;
  double xyz=ZZ.dot(_N.cross(YY));
  if(xyz>=0){return acos(cp);}else{ return 2*3.1415926-acos(cp);}
}

EvtmyEulerAngles::EvtmyEulerAngles( const EvtVector4R & Pyaxis, const EvtVector4R & Pzaxis){
  for (int i=1;i<4;i++){
    _Yaxis.set(i-1,Pyaxis.get(i));
    _Zaxis.set(i-1,Pzaxis.get(i));
  }
  EulerAngles();
}

EvtmyEulerAngles::EvtmyEulerAngles(){
}

double EvtmyEulerAngles::getAlpha(){
  return _alpha;
}

double EvtmyEulerAngles::getBeta(){
  return _beta;
}

double EvtmyEulerAngles::getGamma(){
  return _gamma;
}

void EvtmyEulerAngles::EulerAngles(){
  // to calculate Euler angles with y-convention, i.e. R=R(Z, alpha).R(Y,beta).R(Z,gamma)
  double pi=3.1415926;
    _ry=_Yaxis.d3mag();
    _rz=_Zaxis.d3mag();

    if(_ry==0 ||_rz==0) {
      report(ERROR,"") << "Euler angle calculation specified by zero modules of the axis!"<<endl;
      report(ERROR,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();
    }
   double tolerance=1e-10;
   bool Y1is0=fabs(_Yaxis.get(0))<tolerance;
   bool Y2is0=fabs(_Yaxis.get(1))<tolerance;
   bool Y3is0=fabs(_Yaxis.get(2))<tolerance;
   bool Z1is0=fabs(_Zaxis.get(0))<tolerance;
   bool Z2is0=fabs(_Zaxis.get(1))<tolerance;
   bool Z3is0=fabs(_Zaxis.get(2))<tolerance;

   if(Y1is0 && Y3is0 && Z1is0 && Z2is0 ){
     _alpha=0; _beta=0; _gamma=0;  
   return;
   }

   if( Z1is0 && Z2is0 && !Y2is0){
      _alpha=0; _beta=0;
      _gamma=acos(_Yaxis.get(0)/_ry);
	if(_Yaxis.get(1)<0) _gamma=2*pi - _gamma;
      return;
   }

  // For general case to calculate Euler angles  
  // to calculate beta  

   if(Z1is0 && Z2is0) {
     _beta=0;
   } else{ _beta =acos(_Zaxis.get(2)/_rz);}

   //to calculate alpha

   if(_beta==0){
     _alpha=0.0;
   }else {
     double cosalpha=_Zaxis.get(0)/_rz/sin(_beta);
     if(fabs(cosalpha)>1.0) cosalpha=cosalpha/fabs(cosalpha);
     _alpha=acos(cosalpha);
     if(_Zaxis.get(1)<0.0) _alpha=2*pi - _alpha;
   }

   //to calculate gamma, alpha=0 and beta=0 case has been calculated, so only alpha !=0 and beta !=0 case left

   double singamma=_Yaxis.get(2)/_ry/sin(_beta);
   double cosgamma=(-_Yaxis.get(0)/_ry-cos(_alpha)*cos(_beta)*singamma)/sin(_alpha);
  if(fabs(singamma)>1.0) singamma=singamma/fabs(singamma);
  _gamma=asin(singamma);
  if(singamma>0 && cosgamma<0 ) _gamma=pi - _gamma;  // _gamma>0
  if(singamma<0 && cosgamma<0 ) _gamma=pi - _gamma;  //_gamma<0
  if(singamma<0 && cosgamma>0 ) _gamma=2*pi + _gamma; //_gamma<0
 
  
}
