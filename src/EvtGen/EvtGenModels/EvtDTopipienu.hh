#ifndef EVTDTOPIPIENU_HH
#define EVTDTOPIPIENU_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtDTopipienu:public  EvtDecayProb  {

public:

   EvtDTopipienu() {}
   virtual ~EvtDTopipienu();

   void getName(std::string& name);
   EvtDecayBase* clone();

   void init();
   void initProbMax();
   void decay(EvtParticle *p);

private:
   void KinVGen(EvtVector4R vp4_K, EvtVector4R vp4_Pi, EvtVector4R vp4_Lep, EvtVector4R vp4_Nu, int charm, double& m2, double& q2, double& cosV, double& cosL, double& chi);
   double calPDF(double m2, double q2, double cosV, double cosL, double chi);
   void ResonanceGS(double m, double q, double massD, double massPi1, double massPi2, EvtComplex& F11, EvtComplex& F21, EvtComplex& F31);
   void ResonancePGScbw(double m, double q, EvtComplex& F11, EvtComplex& F21, EvtComplex& F31);
   void ResonanceSBugg(double m, double q, EvtComplex& F10);
   double getPStar(double m, double m1, double m2);
   double getF1(double m, double m0, double m_c1, double m_c2, double rBW);
   double getF2(double m, double m0, double m_c1, double m_c2, double rBW);
   double getWidth0(double m, double m0, double m_c1, double m_c2, double width0);
   double getWidth1(double m, double m0, double m_c1, double m_c2, double width0, double rBW);
   double getWidth2(double m, double m0, double m_c1, double m_c2, double width0, double rBW);
   EvtComplex getCoef(double rho, double phi);

   inline double getGx(double m0, double p0, double m_c1, double m_c2);
   inline double getFx(double mr2, double sx, double p, double p0, double m_c1, double m_c2);
   inline double getHx(double sx, double p, double m_c1, double m_c2);
   inline double getdh(double mr2, double p0, double m_c1, double m_c2);
   inline double getG1(double sx, double mr);
   inline double getZ(double sx, double Mr2);
   inline double getRho(double sx, double mX);
   inline double getWidthrho(double m, double m0, double width0, double p, double p0);
   inline EvtComplex getrho(double sx, double mX);

   int first;
   int last;
   double ProbMax;
   int type[5];

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
   double BF;

   double m0_omega;
   double width0_omega;
   double rho_omega;
   double phi_omega;
   double BF_omega;

   double m0_S;
   double rho_S;
   double phi_S;
      
   double Dp_mD;
   double Dp_mPi1;
   double Dp_mPi2;
   double D0_mD;
   double D0_mPi1;
   double D0_mPi2;

   double Pi;
   double root2;
   double root2d3;
   double root1d2;
   double root3d2;
   double mKa;
   double mPi;
   double mEt;

};

#endif
