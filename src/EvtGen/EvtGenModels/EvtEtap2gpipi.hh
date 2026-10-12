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
// Module:  EvtEtap2gpipi.cc
//
//  for etaprime -> gamma pipi without box contributions
//  see PHYSICAL REVIEW LETTERS 120, 242003 (2018) and refs. therein
//
// Description: Model provided by user, see the mannual
//
// Modification history:
//
//    Ping R.-G.       December, 2006       Module created
//
//------------------------------------------------------------------------

#ifndef EvtETAP2GPIPI_HH 
#define EvtETAP2GPIPI_HH 

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtVector4R.hh"

class EvtParticle;

class EvtEtap2gpipi:public  EvtDecayIncoherent  {

public:
  
  EvtEtap2gpipi() {}
  virtual ~EvtEtap2gpipi();

  void getName(std::string& name);

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *p); 

  double AmplitudeSquare();
  double AmplitudeSquare_p();
  
  EvtVector4R GetDaugMomLab(int i){return _p4Lab[i];}
  EvtVector4R GetDaugMomCM(int i){return _p4CM[i];}
  EvtVector4R GetDaugMomHel(int i){return _p4Hel[i];}


  int GetNdaug(){return _nd;}   

private:

 static int nrun;
 static double max_amps;

 EvtVector4R _p4Lab[10],_p4CM[10],_p4Hel[10];
 int _nd;

};


#endif

