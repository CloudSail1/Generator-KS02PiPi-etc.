//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtKS0IB.hh
//
// Description: 
//      Generator for decay K_S0 -> gamma pi pi.
//      See Study of the Photon Spectrum in the Decay K0(s) –> pi+ pi- gamma
//      DOI: 10.1016/0370-2693(76)90542-6
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//    cloudsail   June 4, 2026            Fork from EvtPi0Dalitz.hh
//
//------------------------------------------------------------------------

#ifndef EVTKS0IB_HH
#define EVTKS0IB_HH

#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtDecayBase.hh"

class EvtParticle;

class EvtKS0IB: public  EvtDecayProb  {

public:

  EvtKS0IB() {}
  virtual ~EvtKS0IB();

  void getName(std::string& name);
  EvtDecayBase* clone();

  void init();
  // comment out to auto-tuning probmax
  // void initProbMax();

  void decay(EvtParticle *p);
};

#endif
