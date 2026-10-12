#ifndef EVTD0TOKPIPIPI_HH
#define EVTD0TOKPIPIPI_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtD0ToKpipipi:public EvtDecayProb  {

public:

   EvtD0ToKpipipi() {}
   virtual ~EvtD0ToKpipipi();
   void getName(std::string& name);
   EvtDecayBase* clone();
   void init();
   void initProbMax();
   void decay(EvtParticle *p);

private:
   EvtComplex KPiSFormfactor(double sa, double sb, double sc, double r);
   double calPDF(double Km[], double Pip1[], double Pip2[], double Pim[]);
   EvtComplex D2VV(double P1[], double P2[], double P3[], double P4[], int g[]);
   EvtComplex D2AP_A2VP(double P1[], double P2[], double P3[], double P4[], int g[], int flag);
   //flag = 1, V = K*; flag = 2, V = rho
   EvtComplex D2AP_A2SP(double P1[], double P2[], double P3[], double P4[], int flag);
   //flag = 1, S = K*; flag = 2, S = rho
   EvtComplex D2PP_P2VP(double P1[], double P2[], double P3[], double P4[], int flag);
   //flag = 1, V = K*; flag = 2, V = rho 
   EvtComplex D2VP_V2VP(double P1[], double P2[], double P3[], double P4[], int flag);
   //flag = 1, (K*Pi)V; flag = 2, (rhoK)V
   EvtComplex D2VS(double P1[], double P2[], double P3[], double P4[], int g, int flag);
   //flag = 1, V = K*; flag = 2, V = rho
   EvtComplex D2TS(double P1[], double P2[], double P3[], double P4[], int flag);
   //flag = 1, T = K*; flag = 2, T = rho
   EvtComplex PHSP(double Km[], double Pip[]);

   double calDalEva(double P1[], double P2[], double P3[]);
   EvtComplex Spin_factor(double P1[], double P2[], double P3[], int spin);
   EvtComplex getProp(double s[], int flag);
   EvtComplex rhofactor(double sx, double sdau);

   double dot(double *a1, double *a2)const;
   double Qabcs(double sa, double sb, double sc)const;
   double barrier(double l, double sa, double sb, double sc, double r)const;
   void calt1(double daug1[], double daug2[], double t1[]) const;
   void calt2(double daug1[], double daug2[], double t2[][4]) const;

   EvtComplex propogator(double mass, double width, double sx) const;
   double wid(double mass, double sa, double sb, double sc, double r, int l) const;
   double h(double m, double q) const;
   double dh(double mass, double q0) const;
   double f(double mass, double sx, double q0, double q) const;
   double d(double mass, double q0) const;
   double Flatte_rhoab(double sa, double sb, double sc) const;
   double rhoab_L(double sa, double sb, double sc, int L) const;
   EvtComplex propagatorRBW(double mass, double width, double sa, double sb, double sc, double r, int l)const;
   EvtComplex propagatorGS(double mass, double width, double sa, double sb, double sc, double r, int l)const;
   EvtComplex propagatorFlatte(double mass, double width, double sx, double *sb, double *sc)const;
   EvtComplex propagatorK1270(double mass, double sx, double *sb, double *sc, double *g)const;
   EvtComplex propagator980(double mass, double sx, double *sb, double *sc, double *g)const;
   EvtComplex propagatorf600(double mass, double width, double sx)const;
   EvtComplex propagatorsigma500(double mass, double width, double sx, double *g)const;

   double G[4][4],E[4][4][4][4];
   double mass[4];
   double width[4];
   double rho[24];
   double phi[24];

   double mD;
   double rRes;
   double rD;
   double metap;
   double mkstr;
   double mk0;
   double mass_Kaon;
   double mass_Pion;
   double mass_Pi0;
   double math_pi;

   double pi;
   double mpi;
   double g1;
   double g2;
};

#endif
