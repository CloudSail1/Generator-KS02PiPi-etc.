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
// Module:  EvtEtap2gpipiB.cc
//
// Description: Model provided by user, see the mannual     
//
// Modification history:
//
//    Ping R.-G.       December, 2006       Module created
//    qinlq       Sep., 2011       Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtPhotonParticle.hh" 
#include "EvtGenBase/EvtNeutrinoParticle.hh"
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtRaritaSchwingerParticle.hh"
#include "EvtGenBase/EvtHighSpinParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtEtap2gpipiB.hh"
#include "EvtGenBase/EvtRandom.hh"
#include <string>
using namespace std;  //::endl;

EvtEtap2gpipiB::~EvtEtap2gpipiB() {}

void EvtEtap2gpipiB::getName(std::string& model_name){

  model_name="Etap2gpipiB";     

}

EvtDecayBase* EvtEtap2gpipiB::clone(){

  return new EvtEtap2gpipiB;

}


void EvtEtap2gpipiB::init(){

// check that there are 1 arguments:angular distribution parameter 
  checkNArg(1);
  checkNDaug(3);
  EvtSpinType::spintype parenttype = EvtPDL::getSpinType(getParentId());
}
void EvtEtap2gpipiB::initProbMax(){
  noProbMax();
}


int EvtEtap2gpipiB::nrun=1;
double EvtEtap2gpipiB::max_amps=0.0;

void EvtEtap2gpipiB::decay( EvtParticle *p ){

double amps,SamAmps,rd1;


// calculated the max amplitude square in 20000 events, is it enough larger?

   if(nrun==1){
     int ir,nd;
     for(ir=0;ir<=2000000;ir++){
    // for(ir=0;ir<=20000;ir++){
  loop0:
       p->initializePhaseSpace(getNDaug(),getDaugs());
       int nd=p->getNDaug(),i;
       _nd=nd;
       for(i=0;i<=nd-1;i++){
          _p4Lab[i]=p->getDaug(i)->getP4Lab();
          _p4CM[i]=p->getDaug(i)->getP4();
       }
       amps=AmplitudeSquare();
       if(amps<0) goto loop0;
       if(amps>max_amps) max_amps=amps*1.01;

       //cout<<"amps= "<<amps<<endl;
       //cout<<"max_amps= "<<max_amps<<endl;

       if(nrun<2){//check decay card
        EvtId pid0= p->getId();//etaprime
        EvtId pid1= p->getDaug(0)->getId();//pi+
        EvtId pid2= p->getDaug(1)->getId();//pi-
        EvtId pid3= p->getDaug(2)->getId();//gamma
        if(pid0!= EvtPDL::getId("eta'")) { std::cout<<"Etap2gpipiB need specify the parent as etaprime"<<std::endl;abort();}
        if(pid3!= EvtPDL::getId("gamma")) { std::cout<<"Etap2gpipiB need specify the daughters  as pi+ pi- gamma"<<std::endl;abort();}
       }
       nrun++;
     }
  }

//cout<<"max_amps= "<<max_amps<<endl; 


if(max_amps==0.0) {report(ERROR,"EvtGen")<<"The decay amplitude square should be positive number"<<endl;abort();}

//cout<<"max_amp="<<max_amps<<endl;

loop:
  p->initializePhaseSpace(getNDaug(),getDaugs());
  int i;
  for(i=0;i<=p->getNDaug()-1;i++){
          _p4Lab[i]=p->getDaug(i)->getP4Lab(); 
          _p4CM[i]=p->getDaug(i)->getP4();
      }
//  Put phase space results into the daughters.
   amps=AmplitudeSquare();
   if(amps < 0) goto loop;
   SamAmps=amps/max_amps;
   rd1=EvtRandom::Flat(0.0, 1.0);
  if(rd1>=SamAmps) goto loop;
  return ;
}


