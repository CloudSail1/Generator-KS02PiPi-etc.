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
// Module: EvtGen/EvtKS0IB.cc
//
// Description: 
//      Generator for decay K_S0 -> gamma pi pi.
//      See Study of the Photon Spectrum in the Decay K0(s) –> pi+ pi- gamma
//      DOI: 10.1016/0370-2693(76)90542-6
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//    cloudsail   June 4, 2026            Fork from EvtPi0Dalitz.cc
//
//------------------------------------------------------------------------

#include "EvtGenModels/EvtKS0IB.hh"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cmath>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtConst.hh"
using std::fstream;

EvtVector3R toVector3R(EvtVector4R p4) {
  return EvtVector3R(p4.get(1),p4.get(2),p4.get(3));
}

EvtVector3R toUnit(EvtVector3R p3) {
  double p2 = p3.get(0)*p3.get(0)
    + p3.get(1)*p3.get(1)
    + p3.get(2)*p3.get(2);
  double p = sqrt(p2);
  return EvtVector3R(p3.get(0)/p, p3.get(1)/p, p3.get(2)/p);
}

EvtKS0IB::~EvtKS0IB() {}

void EvtKS0IB::getName(std::string& model_name){

    model_name="KS0_IB";

}

EvtDecayBase* EvtKS0IB::clone(){

  return new EvtKS0IB;

}

// comment out to auto-tuning probmax
// void EvtKS0IB::initProbMax(){

//   setProbMax(1);

// }


void EvtKS0IB::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

}


void EvtKS0IB::decay( EvtParticle *p){

  EvtVector4R pip, pim, gamma, ks0;
  double weight = p->initializePhaseSpace(getNDaug(),getDaugs(),0.00000002,0,1);
  setWeight(weight);

  pip   = p->getDaug(0)->getP4();
  pim   = p->getDaug(1)->getP4();
  gamma = p->getDaug(2)->getP4();
  ks0   = p->getP4();

  /*
  Refer to paper (see head), PDF defined by:
  PDF(E_gamma, cos(theta))
  =(2*alpha/pi)
    *BranchingFraction(KS_0->pi pi)
    *(1/E_gamma)
    *(1-2*E_gamma/m_KS_0)
    *beta^3*sin(theta)^2
    /beta_0*(1-beta^2*cos(theta)^2)^2

  beta
  =sqrt(
    1-4*m_pi^2
    /m_K^2-2*m_K*m_gamma
  )
  theta is the angle beteew gamma and one pion in dipion frame
  E_gamma is the energy of gamma in KS_0 frame

  just caculate cos(theta) and E_gamma
  */

  // caculate cos(theta)
  EvtVector4R dipion_p4(pip+pim);
  EvtVector4R dipion_gamma = boostTo(gamma, dipion_p4);
  EvtVector4R dipion_pip   = boostTo(pip,   dipion_p4);
  double      cos_theta    = toUnit(toVector3R(dipion_pip))
    .dot(toUnit(toVector3R(dipion_gamma)));

  // caculate E_gamma
  double e_gamma = gamma.get(0); // already in KS_0 rest frame

  // caculate beta
  double beta_0 = sqrt(
    1
    -4*pip.mass2()/ks0.mass2()
  );
  double beta   = sqrt(
    1
    -4*pip.mass2()
    /(ks0.mass2()
      -2*ks0.mass()*e_gamma
    )
  );

  double bf = 0.692; // average branching fraction for KS_0 -> pi pi in PDG 2026
  double pd = // probability density
    (2*(1.0/137)/EvtConst::pi)
    *bf
    *(1/e_gamma)
    *(1-2*e_gamma/ks0.mass())
    *pow(beta,3)
    *(1-cos_theta*cos_theta)
    /(beta_0
      *pow(
        1-beta*beta*cos_theta*cos_theta
        ,2
      )
    );

  report(NOTICE,"EvtGen") << "[EvtKS0IB.cc] e_gamma   : " << e_gamma       << endl;
  report(NOTICE,"EvtGen") << "[EvtKS0IB.cc] cos_theta : " << cos_theta     << endl;
  report(NOTICE,"EvtGen") << "[EvtKS0IB.cc] weight    : " << weight        << endl;
  report(NOTICE,"EvtGen") << "[EvtKS0IB.cc] prob      : " << pd            << endl;
  
  setProb(pd);
  return;
}