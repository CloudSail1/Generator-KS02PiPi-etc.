#ifndef EVTDTOKPIENU_HH
#define EVTDTOKPIENU_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtDToKpienu:public  EvtDecayProb  {

public:

   EvtDToKpienu() {}
   virtual ~EvtDToKpienu();

   void getName(std::string& name);
   EvtDecayBase* clone();

   void init();
   void initProbMax();
   void decay(EvtParticle *p);

private:
   void KinVGen(EvtVector4R vp4_K, EvtVector4R vp4_Pi, EvtVector4R vp4_Lep, EvtVector4R vp4_Nu, int charm, double& m2, double& q2, double& cosV, double& cosL, double& chi);
   double calPDF(double m2, double q2, double cosV, double cosL, double chi);
   void ResonanceP(double m, double q, double mV, double mA, double V_0, double A1_0, double A2_0, double m0, double width0, double rBW, double& amplitude, double& delta, EvtComplex& F11, EvtComplex& F21, EvtComplex& F31);
   void NRS(double m, double q, double rS, double rS1, double a_delta, double b_delta, double mA, double m0, double width0, double& amplitude, double& delta, EvtComplex& F10);
   void ResonanceD(double m, double q, double mV, double mA, double TV_0, double T1_0, double T2_0, double m0, double width0, double rBW, double& amplitude, double& delta, EvtComplex& F12, EvtComplex& F22, EvtComplex& F32);
   double getPStar(double m, double m1, double m2);
   double getF1(double m, double m0, double m_c1, double m_c2, double rBW);
   double getF2(double m, double m0, double m_c1, double m_c2, double rBW);
   double getWidth0(double m, double m0, double m_c1, double m_c2, double width0);
   double getWidth1(double m, double m0, double m_c1, double m_c2, double width0, double rBW);
   double getWidth2(double m, double m0, double m_c1, double m_c2, double width0, double rBW);
   EvtComplex getCoef(double rho, double phi);

   int nAmps;
   int type[5];

   double rS;
   double rS1;
   double a_delta;
   double b_delta;
   double m0_1430_S;
   double width0_1430_S;

   double mV;
   double mA;
   double V_0;
   double A1_0;
   double A2_0;
   double m0;
   double width0;
   double rBW;
   double rho;
   double phi;
   double m0_1410;
   double width0_1410;
   double rho_1410;
   double phi_1410;
   double TV_0;
   double T1_0;
   double T2_0;
   double m0_1430;
   double width0_1430;
   double rho_1430;
   double phi_1430;

   double mD;
   double mPi;
   double mK;
   double Pi;
   double root2;
   double root2d3;
   double root1d2;
   double root3d2;
};

#endif
