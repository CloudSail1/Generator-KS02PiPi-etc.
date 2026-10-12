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
// Module: EvtDIY.cc
//
// Description: Routine for users to create model to use the model DIY 
//
// Modification history:
//
//    Ping R.-G.       December, 2006       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <math.h>
#include "RooComplex.h"
#include <complex>
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtKine.hh"

#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtHelSys.hh"
#include "EvtGenModels/EvtEtap2gpipiB.hh"//cov form, not normalized, get momentum in CMS
#include <string>
#include "TMath.h"
using namespace std;  



// ******** eta' decay to pi+ pi- gamma with angular distribution
class angular_boxrho2{
public:
  angular_boxrho2(EvtVector4R pd1, EvtVector4R pd2,EvtVector4R pd3, bool type_switch){
    _pd[0]=pd1;
    _pd[1]=pd2;
    _pd[2]=pd3;
    ang_on=type_switch;
  }
  double lamb(double a, double b, double c);

  double amps(int xx);

private:
  EvtVector4R _pd[3];
  bool ang_on;
};

  double angular_boxrho2::amps(int flag){
      EvtVector4R prho=_pd[0]+_pd[1];
      EvtVector4R mprho(prho.get(0),-1*prho.get(1),-1*prho.get(2),-1*prho.get(3));
      //cout<<"pos:"<<prho.get(0)<<", "<<prho.get(1)<<", "<<prho.get(2)<<", "<<prho.get(3)<<endl;
      //cout<<"neg:"<<mprho.get(0)<<", "<<mprho.get(1)<<", "<<mprho.get(2)<<", "<<mprho.get(3)<<endl;
      double m2=(_pd[0]+_pd[1]).mass2(); 
      double m=(_pd[0]+_pd[1]).mass(); 
      //*****************angluar part******************
      EvtVector4R pi_rhocms = boostTo(_pd[0],mprho);
      EvtVector4R gam_rhocms = boostTo(_pd[2],mprho);
      EvtVector4R etap_rhocms = boostTo(prho+_pd[2],mprho);
      double coshel_pi = pi_rhocms.dot(etap_rhocms)/pi_rhocms.d3mag()/etap_rhocms.d3mag();
      //double coshel_pi2 = pi_rhocms.dot(gam_rhocms)/pi_rhocms.d3mag()/gam_rhocms.d3mag();
    //  cout<<(prho+_pd[2]).get(1)<<", "<<(prho+_pd[2]).get(2)<<", "<<(prho+_pd[2]).get(3)<<", "<<(prho+_pd[2]).get(0)<<endl;
      //cout<<"hel1:"<<coshel_pi<<endl;
      //cout<<"hel2:"<<coshel_pi2<<endl;
      double ang_part = 1-coshel_pi*coshel_pi;
      //*****************angluar part******************
      //                                             //VALUE      // ERROR 
      //the values of parameters in model=rho0+omega+rho(1405)+box are obtained from fitting to my data
      double Mrho  =7.76565e-01;       
      double Grho  =1.50839e-01;      
      double delta =1.60865e-03;       
      double argdelta= 6.71184e-02; 
      double delta2=   4.15936e-01 ; 
      double argdelta2= 2.00844e-02;
      double Eetap=   -2.13798e+01; 
      
      double Momega = 0.78265;
      double Gomega = 0.00849;
      double Mrhop = 1.465;
      double Grhop = 0.400;

      double PI = 3.14159265;
      double mpi=0.13957018;
      //double kgamma = (0.95778*0.95778-m2)/(2*0.95778);   
      //double qpim = sqrt(0.25*m2 - pow(0.13957018,2));
      double kgamma=0;
      if(m<0.95778) kgamma = (0.95778*0.95778-m*m)/(2*0.95778);
      double qpim=0;
      if(m>2*mpi) qpim = sqrt(0.25*m*m - pow(mpi,2));
      
      double qpimrho = sqrt(0.25*Mrho*Mrho - pow(mpi,2));
      double qpimrhop = sqrt(0.25*Mrhop*Mrhop - pow(mpi,2));

      double hm = (2./PI)*(qpim/m)*log((m+2*qpim)/(2*mpi));

      double hmrho  = (2./PI)*(qpimrho/Mrho)*log((Mrho+2*qpimrho)/(2.*mpi));
      double hmrhop = (2./PI)*(qpimrhop/Mrhop)*log((Mrhop+2*qpimrhop)/(2.*mpi));

      double d=(3./PI)*(mpi*mpi/pow(qpimrho,2))*log((Mrho+2*qpimrho)/(2*mpi))+(Mrho/(2*PI*qpimrho))-((mpi*mpi*Mrho)/(PI*pow(qpimrho,3)));
      double drhop=(3./PI)*(mpi*mpi/pow(qpimrhop,2))*log((Mrhop+2*qpimrhop)/(2*mpi))+(Mrhop/(2*PI*qpimrhop))-((mpi*mpi*Mrhop)/(PI*pow(qpimrhop,3)));

      double dhds = hmrho*(pow(8*pow(qpimrho,2),-1)-pow(2*Mrho*Mrho,-1))+pow(2*PI*Mrho*Mrho,-1);
      double dhdsrhop = hmrhop*(pow(8*pow(qpimrhop,2),-1)-pow(2*Mrhop*Mrhop,-1))+pow(2*PI*Mrhop*Mrhop,-1);

      double fm2 = Grho*(Mrho*Mrho/pow(qpimrho,3))*(qpim*qpim*(hm-hmrho)+(Mrho*Mrho-m*m)*qpimrho*qpimrho*(dhds));
      double fm2rhop = Grhop*(Mrhop*Mrhop/pow(qpimrhop,3))*(qpim*qpim*(hm-hmrhop)+(Mrhop*Mrhop-m*m)*qpimrhop*qpimrhop*(dhdsrhop));

      double Grhom = Grho*pow((qpim/qpimrho),3)*(Mrho/m);
      double Grhopm = Grhop*pow((qpim/qpimrhop),3)*(Mrhop/m);

      double coefficient = (1./(48*PI*PI*PI))*pow(kgamma,2)*pow(qpim,2);
      double coe1= pow((2*sqrt(48*PI*pow(Mrho,-4))),2);
      
      RooComplex denomGrhom(Mrho*Mrho-m*m+fm2,-1*Mrho*Grhom);
      RooComplex real(Mrho*Mrho*(1+d*Grho/Mrho),0);
      RooComplex BWrho= real/denomGrhom;
     
      RooComplex denomBWomega(Momega*Momega-m*m,-Momega*Gomega);
      RooComplex real1(Momega*Momega,0);
      RooComplex BWomega = real1/denomBWomega;
      RooComplex cdelta(delta*cos(argdelta),delta*sin(argdelta));
      RooComplex part2=BWrho*cdelta*BWomega*(m*m)/(Momega*Momega);

      RooComplex denomGrhopm(Mrhop*Mrhop-m*m+fm2rhop,-1*Mrhop*Grhopm);
      RooComplex realpm(Mrhop*Mrhop*(1+drhop*Grhop/Mrhop),0);
      RooComplex BWrhop = realpm/denomGrhopm;
      RooComplex cdelta2(delta2*cos(argdelta2),delta2*sin(argdelta2));
      RooComplex denom(delta2*cos(argdelta2)+1.0,delta2*sin(argdelta2));
      RooComplex part3=cdelta2*BWrhop;
     
      RooComplex eof(Eetap,0);
      RooComplex amp=(( BWrho* (cdelta*BWomega*(m*m)/(Momega*Momega)+1.0) + cdelta2*BWrhop)/(cdelta2+1)*sqrt(coe1)+ eof);
      double total;
      if(flag==0)total=coefficient* amp.abs2();

      //DIY for omega in case of rho0-omega-rho'-box
      else if(flag==2) total=coefficient*coe1* (part2/(cdelta2+1)).abs2();
    
      //DIY for rho' in case of rho0-omega-rho'-box
       else if(flag==3) total=coefficient*coe1* (part3/(cdelta2+1)).abs2();
      
       //DIY for box in case of rho0-omega-rho'-box
       else if(flag==4) total=coefficient* (eof.abs2());
      
      double amp2;
      if(ang_on) amp2 = total*ang_part;
      else amp2 = total;
      return amp2;
}
//  ************** end class angular_boxrho2


//////////////////***** DIYshero  Model
double EvtEtap2gpipiB::AmplitudeSquare(){
  EvtVector4R dp1=GetDaugMomCM(0);
  EvtVector4R dp2=GetDaugMomCM(1);
  EvtVector4R dp3=GetDaugMomCM(2);
  //cout<<"CMS:"<<(dp1+dp2+dp3).get(0)<<","<<(dp1+dp2+dp3).get(1)<<","<<(dp1+dp2+dp3).get(2)<<","<<(dp1+dp2+dp3).get(3)<<endl;
  //cout<<(dp1.get(0)+dp2.get(0)+dp3.get(0))<<endl;
  //cout<<(dp1.get(1)+dp2.get(1)+dp3.get(1))<<endl;
  //cout<<(dp1.get(2)+dp2.get(2)+dp3.get(2))<<endl;
  //cout<<(dp1.get(3)+dp2.get(3)+dp3.get(3))<<endl;
  angular_boxrho2 Rho2pipi(dp1,dp2,dp3,true);
  //angular_boxrho2 Rho2pipi(dp1,dp2,dp3,false);
  int xbd=getArg(0);
  //cout<<"xbd= "<<xbd<<" amps= "<<Rho2pipi.amps(xbd) <<endl;
  return Rho2pipi.amps(xbd);
}

double EvtEtap2gpipiB::AmplitudeSquare_p(){
  EvtVector4R dp1=GetDaugMomCM(0);
  EvtVector4R dp2=GetDaugMomCM(1);
  EvtVector4R dp3=GetDaugMomCM(2);
  angular_boxrho2 tmp(dp1,dp2,dp3,false);
  int xbd=getArg(0);
  return tmp.amps(xbd);
}
