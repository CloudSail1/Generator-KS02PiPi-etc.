#ifndef EVTDTOKSPIPIPI_HH
#define EVTDTOKSPIPIPI_HH

#include "EvtGenBase/EvtDecayProb.hh"
class EvtParticle;

class EvtDToKSpipipi:public EvtDecayProb  {

public:

   EvtDToKSpipipi() {}
   virtual ~EvtDToKSpipipi();

   void getName(std::string& name);
   EvtDecayBase* clone();

   void init();
   void initProbMax();
   void decay(EvtParticle *p);

private:
   EvtComplex KPiSFormfactor(double sa, double sb, double sc, double r);
   double calPDF(double Km[], double Pip1[], double Pip2[], double Pim[]);
   EvtComplex D2AP_A2VP(double P1[], double P2[], double P3[], double P4[], int L);
   EvtComplex D2AP_A2SP(double P1[], double P2[], double P3[], double P4[]);
   EvtComplex D2PP_P2VP(double P1[], double P2[], double P3[], double P4[]);
   EvtComplex D2VP_V2VP(double P1[], double P2[], double P3[], double P4[]);
   EvtComplex PHSP(double Km[], double Pip[]);
   EvtComplex getprop(double daug1[], double daug2[], double mass, double width, int flag, int L);
   EvtComplex rhoab(double sa, double sb, double sc);
   EvtComplex rho4Pi(double sa);

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
   double mrho, Grho, mKstr, GKstr, msigma, Gsigma;
   double ma1, Ga1, mK1270, GK1270, mK1400, GK1400, mK1460, GK1460, mK1650, GK1650;
   double rho_omega, phi_omega;

   double mD; 
   double rD; 
   double metap;
   double mkstr;
   double mk0;
   double mass_Kaon;
   double mass_Pion;
   double math_pi;

   double pi;
   double mpi;
   double g1;
   double g2;

};

#endif
