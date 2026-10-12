//--------------------------------------------------------------------------
// Environment:
//      This software is part of models developed at BES collaboration
//      based on the EvtGen framework.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/BesCopyright
//      Copyright (A) 2006      Ping Rong-Gang @IHEP
//
// Module: EvtDToKSpipipi.cc
//         the necessary file: EvtDToKSpipipi.hh
//
// Description: D -> pi pi e+ nu, 
//              PHYSICAL REVIEW Letter 122, 062001 (2019)
//
// Modification history:
//   Liaoyuan Dong    Jan.16, 2020       Module created
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtDTopipienu.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include <stdlib.h>

EvtDTopipienu::~EvtDTopipienu() {}

void EvtDTopipienu::getName(std::string& model_name){
   model_name="DTopipienu";
}

EvtDecayBase* EvtDTopipienu::clone(){
   return new EvtDTopipienu;
}

void EvtDTopipienu::init(){
   static EvtId DM=EvtPDL::getId("D-");
   static EvtId DP=EvtPDL::getId("D+");
   static EvtId D0=EvtPDL::getId("D0");
   static EvtId D0B=EvtPDL::getId("anti-D0");
        
   checkNArg(0);
   checkNDaug(4);
   checkSpinParent(EvtSpinType::SCALAR);
   checkSpinDaughter(0,EvtSpinType::SCALAR);
   checkSpinDaughter(1,EvtSpinType::SCALAR);

   std::cout << "EvtDTopipienu ==> Initialization !" << std::endl;

   EvtId parnum=getParentId();
   if ( parnum == D0 || parnum == D0B ) {
      first = 0;
      last  = 1;
      ProbMax = 405000;
   } else if (parnum == DP || parnum == DM ) {
      first = 1;
      last  = 4;
      ProbMax = 591000;      
   }

   type[0] = 0;
   type[1] = 1;
   type[2] = 2;
   type[3] = 3;

   mV = 2.01;
   mA = 2.42;
   V_0 = 1.6948;
   A1_0 = 1;
   A2_0 = 0.84489;

   m0 = 0.77526;
   width0 = 0.14910;
   rBW = 3.0;
   rho = 1.0;
   phi = 0.0;
   BF  = 1.0;

   m0_omega = 0.78265;
   width0_omega = 0.00849;
   rho_omega = 0.12902;
   phi_omega = 2.9285;
   BF_omega  = 1.0;

   m0_S  =    0.953;
   rho_S =  135.27;
   phi_S =    3.4044;
      
   Dp_mD   =  1.86962;
   Dp_mPi1 =  0.13957;
   Dp_mPi2 =  0.13957;
   D0_mD   =  1.86486;
   D0_mPi1 =  0.13957;
   D0_mPi2 =  0.1349766;

   Pi = atan2(0.0,-1.0);
   root2 = sqrt(2.);
   root2d3 = sqrt(2./3);
   root1d2 = sqrt(0.5);
   root3d2 = sqrt(1.5);

   mKa = 0.493677;
   mPi = 0.13957;
   mEt = 0.547853;
}

void EvtDTopipienu::initProbMax() {
   cout << "EvtDTopipienu: setProbMax = " << ProbMax << endl;
   setProbMax(ProbMax);
}

void EvtDTopipienu::decay( EvtParticle *p){
/*
   double maxprob = 0.0;
   for(int ir=0;ir<=60000000;ir++){
      p->initializePhaseSpace(getNDaug(),getDaugs());
      EvtVector4R _pi1 = p->getDaug(0)->getP4();
      EvtVector4R _pi2 = p->getDaug(1)->getP4();
      EvtVector4R _e   = p->getDaug(2)->getP4();
      EvtVector4R _nu  = p->getDaug(3)->getP4();

      int pid = EvtPDL::getStdHep(p->getDaug(0)->getId());
      int charm;
      if(pid == -211) charm = 1;
      else charm = -1;
      double m2, q2, cosV, cosL, chi;
      KinVGen(_pi1, _pi2, _e, _nu, charm, m2, q2, cosV, cosL, chi);
      double _prob = calPDF(m2, q2, cosV, cosL, chi);
      if(_prob>maxprob) {
          maxprob=_prob;
          std::cout << "Max PDF = " << ir << " charm= " << charm << " prob= " << _prob << std::endl;
      }
   } 
   std::cout << "Max!!!!!!!!!!! " << maxprob<< std::endl;
*/
   p->initializePhaseSpace(getNDaug(),getDaugs());
   EvtVector4R pi1 = p->getDaug(0)->getP4();
   EvtVector4R pi2 = p->getDaug(1)->getP4();
   EvtVector4R e   = p->getDaug(2)->getP4();
   EvtVector4R nu  = p->getDaug(3)->getP4();

   int pid = EvtPDL::getStdHep(p->getDaug(0)->getId());
   int charm;
   if(pid == -211) charm = 1;
   else charm = -1;
   double m2, q2, cosV, cosL, chi;
   KinVGen(pi1, pi2, e, nu, charm, m2, q2, cosV, cosL, chi);
   double prob = calPDF(m2, q2, cosV, cosL, chi);
   setProb(prob);
   return;
}

void EvtDTopipienu::KinVGen(EvtVector4R vp4_K, EvtVector4R vp4_Pi, EvtVector4R vp4_Lep, EvtVector4R vp4_Nu, int charm, double& m2, double& q2, double& cosV, double& cosL, double& chi)
{
    EvtVector4R vp4_KPi = vp4_K + vp4_Pi;
    EvtVector4R vp4_W = vp4_Lep + vp4_Nu;
    m2 = vp4_KPi.mass2();
    q2 = vp4_W.mass2();

    EvtVector4R boost;
    boost.set(vp4_KPi.get(0), -vp4_KPi.get(1), -vp4_KPi.get(2), -vp4_KPi.get(3));
    EvtVector4R vp4_Kp = boostTo(vp4_K, boost);
    cosV = vp4_Kp.dot(vp4_KPi)/(vp4_Kp.d3mag()*vp4_KPi.d3mag());

    boost.set(vp4_W.get(0), -vp4_W.get(1), -vp4_W.get(2), -vp4_W.get(3));
    EvtVector4R vp4_Lepp = boostTo(vp4_Lep, boost);
    cosL = vp4_Lepp.dot(vp4_W)/(vp4_Lepp.d3mag()*vp4_W.d3mag());

    EvtVector4R V = vp4_KPi/vp4_KPi.d3mag();
    EvtVector4R C = vp4_Kp.cross(V);
    C/=C.d3mag();
    EvtVector4R D = vp4_Lepp.cross(V);
    D/=D.d3mag();
    double sinx = C.cross(V).dot(D);
    double cosx = C.dot(D);
    chi = sinx>0? acos(cosx): -acos(cosx);
    if(charm==-1) chi=-chi;
}

double EvtDTopipienu::calPDF(double m2, double q2, double cosV, double cosL, double chi) {
    double m = sqrt(m2);
    double q = sqrt(q2); 

    //begin to calculate form factor
    EvtComplex F10(0.0,0.0);
    EvtComplex F11(0.0,0.0);
    EvtComplex F21(0.0,0.0);
    EvtComplex F31(0.0,0.0);
    EvtComplex F12(0.0,0.0);
    EvtComplex F22(0.0,0.0);
    EvtComplex F32(0.0,0.0);

    EvtComplex f10(0.0,0.0);
    EvtComplex f11(0.0,0.0);
    EvtComplex f21(0.0,0.0);
    EvtComplex f31(0.0,0.0);    

    EvtComplex coef(0.0,0.0);

    for(int index=first; index<last; index++)
    {
        switch(type[index])
        {
            case 0: //calculate form factor of P wave (rho)
                {
                    ResonanceGS(m, q, D0_mD, D0_mPi1, D0_mPi2, f11, f21, f31);
                    coef = getCoef(rho, phi);
                    F11 = F11+ coef*f11;
                    F21 = F21+ coef*f21;
                    F31 = F31+ coef*f31;
                    break;
                }
            case 1: //calculate form factor of P wave (GS)
                {
                    ResonanceGS(m, q, Dp_mD, Dp_mPi1, Dp_mPi2, f11, f21, f31);
                    coef = getCoef(rho, phi);
                    F11 = F11+ coef*f11;
                    F21 = F21+ coef*f21;
                    F31 = F31+ coef*f31;
                    break;
                }
            case 2: //calculate form factor of omega
                {
                    ResonancePGScbw(m, q, f11, f21, f31);
                    coef = getCoef(rho_omega, phi_omega);
                    F11 = F11+ coef*f11;
                    F21 = F21+ coef*f21;
                    F31 = F31+ coef*f31;
                    break;
                }
            case 3: //calculate form factor of S wave (BW corrected by Bugg)
                {
                    ResonanceSBugg(m, q, f10);
                    coef = getCoef(rho_S, phi_S);
                    F10 = F10 + coef*f10;
                    break;
                }
            default:
                {
                    std::cout<<"No such form factor type!!!"<<std::endl;
                    break;
                }
        }
    }

    //begin to calculate pdf value
    double I,I1,I2,I3,I4,I5,I6,I7,I8,I9;

    double cosV2 = cosV*cosV;
    double sinV = sqrt(1.0-cosV2);
    double sinV2 = sinV*sinV;

    EvtComplex F1 = F10 + F11*cosV + F12*(1.5*cosV2-0.5);
    EvtComplex F2 = F21*root1d2 + F22*cosV*root3d2;
    EvtComplex F3 = F31*root1d2 + F32*cosV*root3d2;

    I1 =  0.25*( abs2(F1) + 1.5*sinV2*(abs2(F2) + abs2(F3) ) );
    I2 = -0.25*( abs2(F1) - 0.5*sinV2*(abs2(F2) + abs2(F3) ) );
    I3 = -0.25*( abs2(F2) - abs2(F3) )*sinV2;
    I4 = real( conj(F1)*F2 )*sinV*0.5;
    I5 = real( conj(F1)*F3 )*sinV;
    I6 = real( conj(F2)*F3 )*sinV2;
    I7 = imag( conj(F2)*F1 )*sinV;
    I8 = imag( conj(F3)*F1 )*sinV*0.5;
    I9 = imag( conj(F3)*F2 )*sinV2*(-0.5);
                                                                            
    double coschi = cos(chi);
    double sinchi = sin(chi);
    double sin2chi = 2.0*sinchi*coschi;
    double cos2chi = 1.0 - 2.0*sinchi*sinchi;

    double sinL = sqrt(1.-cosL*cosL);
    double sinL2 = sinL*sinL;
    double sin2L = 2.0*sinL*cosL;
    double cos2L = 1.0 - 2.0*sinL2;

    I = I1 + I2*cos2L + I3*sinL2*cos2chi + I4*sin2L*coschi + I5*sinL*coschi + I6*cosL + I7*sinL*sinchi + I8*sin2L*sinchi + I9*sinL2*sin2chi;
    return I;
}

void EvtDTopipienu::ResonanceGS(double m, double q, double massD, double massPi1, double massPi2, EvtComplex& F11, EvtComplex& F21, EvtComplex& F31)
{
    double pKPi = getPStar(massD, m, q);
    double mD2 = massD*massD;
    double m2 = m*m;
    double m02 = m0*m0;
    double q2 = q*q;
    double mV2 = mV*mV;
    double mA2 = mA*mA;
    double summDm = massD+m;
    double V  = V_0 /(1.0-q2/(mV2));
    double A1 = A1_0/(1.0-q2/(mA2));
    double A2 = A2_0/(1.0-q2/(mA2));
    double A = summDm*A1;
    double B = 2.0*massD*pKPi/summDm*V;

    //construct the helicity form factor
    double H0 = 0.5/(m*q)*((mD2-m2-q2)*summDm*A1-4.0*(mD2*pKPi*pKPi)/summDm*A2);
    double Hp = A-B;
    double Hm = A+B;

    //calculate alpha
    //double B_Kstar = 2./3.0; //B_Kstar = Br(Kstar(892)->k pi)
    double pStar  = getPStar(m, massPi1, massPi2);
    double pStar0 = getPStar(m0, massPi1, massPi2);
    double alpha = sqrt(3.0*Pi*BF/(pStar0*width0));

    //construct amplitudes of (non)resonance
    double F = getF1(m, m0, massPi1, massPi2, rBW);
    EvtComplex C(m02*(1.0+width0*getGx(m0, pStar0, massPi1, massPi2)/m0)*F, 0.0);
    double AA = m02-m2+width0*getFx(m02, m2, pStar, pStar0, massPi1, massPi2);
    double BB = -m0*getWidthrho(m, m0, width0, pStar, pStar0);
    EvtComplex tmp(AA,BB);
    EvtComplex amp = C/tmp;

    double alpham2 = alpha*2.0;
    F11 = amp*alpham2*q*H0*root2;
    F21 = amp*alpham2*q*(Hp+Hm);
    F31 = amp*alpham2*q*(Hp-Hm);
}

void EvtDTopipienu::ResonancePGScbw(double m, double q, EvtComplex& F11, EvtComplex& F21, EvtComplex& F31)
{
    double pKPi = getPStar(Dp_mD, m, q);
    double mD2 = Dp_mD*Dp_mD;
    double m2 = m*m;
    double m02 = m0_omega*m0_omega;
    double mR2 = m0*m0;
    double q2 = q*q;
    double mV2 = mV*mV;
    double mA2 = mA*mA;
    double summDm = Dp_mD+m;
    double V  = V_0 /(1.0-q2/(mV2));
    double A1 = A1_0/(1.0-q2/(mA2));
    double A2 = A2_0/(1.0-q2/(mA2));
    double A = summDm*A1;
    double B = 2.0*Dp_mD*pKPi/summDm*V;
    //construct the helicity form factor
    double H0 = 0.5/(m*q)*((mD2-m2-q2)*summDm*A1-4.0*(mD2*pKPi*pKPi)/summDm*A2);
    double Hp = A-B;
    double Hm = A+B;

    //calculate alpha
    //double B_Kstar = 2./3.0; //B_Kstar = Br(Kstar(892)->k pi)
    double pStar  = getPStar(m, Dp_mPi1, Dp_mPi2);
    double pStar0 = getPStar(m0_omega, Dp_mPi1, Dp_mPi2);
    double alpha = sqrt(3.0*Pi*BF_omega/(pStar0*width0_omega));
 
    //construct amplitudes of (non)resonance
    double F = getF1(m, m0_omega, Dp_mPi1, Dp_mPi2, rBW);
    EvtComplex amp1(0.0, 0.0);
    EvtComplex amp2(0.0, 0.0);
    // CBW
    EvtComplex C1(m0_omega*width0_omega*F,0.0);
    double AA1 = m02-m2;
    double BB1 = -m0_omega*width0_omega;
    EvtComplex tmp1(AA1,BB1);
    amp1 = C1/tmp1;
    // GS
    pStar0 = getPStar(m0, Dp_mPi1, Dp_mPi2);
    EvtComplex C2(mR2*(1.0+width0*getGx(m0, pStar0, Dp_mPi1, Dp_mPi2)/m0), 0.0);
    double AA2 = mR2-m2+width0*getFx(mR2, m2, pStar, pStar0, Dp_mPi1, Dp_mPi2);
    double BB2 = -m0*getWidthrho(m, m0, width0, pStar, pStar0);
    EvtComplex tmp2(AA2,BB2);
    amp2 = C2/tmp2;
    EvtComplex amp = amp1*amp2;

    double alpham2 = alpha*2.0;
    F11 = amp*alpham2*q*H0*root2;
    F21 = amp*alpham2*q*(Hp+Hm);
    F31 = amp*alpham2*q*(Hp-Hm);
}

void EvtDTopipienu::ResonanceSBugg(double m, double q, EvtComplex& F10)
{
    double pKPi = getPStar(Dp_mD, m, q);
    double m2  = m*m;
    double q2  = q*q;
    double mA2 = mA*mA;
   
    double sA    = 0.41*mPi*mPi;
    double mr    = m0_S;         //0.953;
    double mr2   = m0_S*m0_S;    //0.908209;// 0.953*0.953;
    double alpha = 1.3;
    double g4pi  = 0.011;

    EvtComplex ciR(1.0, 0.0);
    EvtComplex ciM(0.0, 1.0);
    EvtComplex Gamma1(0.0, 0.0);
    EvtComplex Gamma2(0.0, 0.0);
    EvtComplex Gamma3(0.0, 0.0);
    EvtComplex Gamma4(0.0, 0.0);

    Gamma1 = getrho(m2,mPi)*getG1(m2,mr)*(m2-sA)/(mr2-sA);
    Gamma2 = getrho(m2,mKa)*0.6*getG1(m2,mr)*(m2/mr2)*exp(-alpha*fabs(m2-4.0*mKa*mKa));
    Gamma3 = getrho(m2,mEt)*0.2*getG1(m2,mr)*(m2/mr2)*exp(-alpha*fabs(m2-4.0*mEt*mEt));
    if (m > 4*mPi) Gamma4 = ciR*mr*g4pi*(1.0+exp(7.082-2.845*mr2))/(1.0+exp(7.082-2.845*m2));

    double AA = mr2-m2-getG1(m2,mr)*(m2-sA)*getZ(m2, mr2)/(mr2-sA);
    EvtComplex amp = ciR/(ciR*AA-ciM*(Gamma1+Gamma2+Gamma3+Gamma4));
    F10 = amp*pKPi*Dp_mD/(1.0-q2/mA2);
}

double EvtDTopipienu::getPStar(double m, double m1, double m2)
{
    double s = m*m;
    double s1 = m1*m1;
    double s2 = m2*m2;
    double x = s + s1 - s2;
    double t = 0.25*x*x/s - s1;
    double p;
    if (t>0.0) {
        p = sqrt(t);
    } else {
        std::cout << " Hello, pstar is less than 0.0" << std::endl;
        p = 0.04;
    }
    return p;
}

double EvtDTopipienu::getF1(double m, double m0, double m_c1, double m_c2, double rBW)
{
    double pStar  = getPStar(m, m_c1, m_c2);
    double pStar0 = getPStar(m0, m_c1, m_c2);
    double rBW2 = rBW*rBW;
    double pStar2  = pStar*pStar;
    double pStar02 = pStar0*pStar0;
    double B  = 1./sqrt(1.+rBW2*pStar2);
    double B0 = 1./sqrt(1.+rBW2*pStar02);
    double F = pStar/pStar0*B/B0;
    return F; 
}

double EvtDTopipienu::getF2(double m, double m0, double m_c1, double m_c2, double rBW)
{
    double pStar  = getPStar(m, m_c1, m_c2);
    double pStar0 = getPStar(m0, m_c1, m_c2);
    double rBW2 = rBW*rBW;
    double pStar2  = pStar*pStar;
    double pStar02 = pStar0*pStar0;
    double B  = 1./sqrt((rBW2*pStar2-3.)*(rBW2*pStar2-3.)+9.*rBW2*pStar2);
    double B0 = 1./sqrt((rBW2*pStar02-3.)*(rBW2*pStar02-3.)+9.*rBW2*pStar02);
    double F = pStar2/pStar02*B/B0;
    return F;
}

double EvtDTopipienu::getWidth0(double m, double m0, double m_c1, double m_c2, double width0)
{
    double pStar = getPStar(m, m_c1, m_c2);
    double pStar0 = getPStar(m0, m_c1, m_c2);
    double width = width0*pStar/pStar0*m0/m;
    return width;
}

double EvtDTopipienu::getWidth1(double m, double m0, double m_c1, double m_c2, double width0, double rBW)
{
    double pStar = getPStar(m, m_c1, m_c2);
    double pStar0 = getPStar(m0, m_c1, m_c2);
    double F = getF1(m, m0, m_c1, m_c2, rBW);
    double width = width0*pStar/pStar0*m0/m*F*F;
    return width;
}

double EvtDTopipienu::getWidth2(double m, double m0, double m_c1, double m_c2, double width0, double rBW)
{
    double pStar = getPStar(m, m_c1, m_c2);
    double pStar0 = getPStar(m0, m_c1, m_c2);
    double F = getF2(m, m0, m_c1, m_c2, rBW);
    double width = width0*pStar/pStar0*m0/m*F*F;
    return width;
}

EvtComplex EvtDTopipienu::getCoef(double rho, double phi)
{
    EvtComplex coef ( rho*cos(phi), rho*sin(phi) );
    return coef;
}

inline double EvtDTopipienu::getGx(double m0, double p0, double m_c1, double m_c2)
{
    double Gg = 0;
    double MPI = 0.5*(m_c1+m_c2);
    Gg = 3*MPI*MPI*log((m0+2*p0)/(2*MPI))/(Pi*p0*p0) + m0/(2*Pi*p0) - MPI*MPI*m0/(Pi*p0*p0*p0);
    return Gg;
}
	
inline double EvtDTopipienu::getFx(double mr2, double sx, double p, double p0, double m_c1, double m_c2)
{
    double Fx = 0;
    Fx = mr2/(pow(p0,3))*(p*p*(getHx(sx,p,m_c1,m_c2)-getHx(mr2,p0,m_c1,m_c2))+(mr2-sx)*p0*p0*getdh(mr2,p0,m_c1,m_c2));
    return Fx;
}

inline double EvtDTopipienu::getHx(double sx, double p, double m_c1, double m_c2)
{
    double m  = sqrt(sx); 
    double Hx = 0;
    Hx = 2*p*log((m+2*p)/(m_c1+m_c2))/(Pi*m);
    return Hx;
}

inline double EvtDTopipienu::getdh(double mr2, double p0, double m_c1, double m_c2)
{
  double mass = sqrt(mr2);
  double dh   = getHx(mass,p0,m_c1,m_c2)*(1.0/(8*p0*p0)-1.0/(2*mass*mass))+1.0/(2*Pi*mass*mass);
  return dh;
}

inline double EvtDTopipienu::getG1(double m2, double Mr)
{
    double b1  = 1.302;
    double b2  = 0.340;
    double A   = 2.426;
    double Mr2 = Mr*Mr;// 0.953*0.953;
    double gg1 = Mr*(b1+b2*m2)*exp(-(m2-Mr2)/A);
    return gg1;
}

inline double EvtDTopipienu::getZ(double m2, double Mr2)
{
    double zz = (getRho(m2,mPi)*log((1.0-getRho(m2,mPi))/(1.0+getRho(m2,mPi)))
                -getRho(Mr2,mPi)*log((1.0-getRho(Mr2,mPi))/(1.0+getRho(Mr2,mPi))))/Pi;

    return zz;
}

inline double EvtDTopipienu::getRho(double m2, double mX)
{
    double rho = 0.0;
    if ((1.0-4.0*mX*mX/m2)>0)
         rho = sqrt(1.0-4.0*mX*mX/m2);
    else rho = 0.0; 
    return rho;
}

inline EvtComplex EvtDTopipienu::getrho(double m2, double mX)
{
    EvtComplex rho(0.0, 0.0);
    EvtComplex ciR(1.0, 0.0);
    EvtComplex ciM(0.0, 1.0);
    if ((1.0-4.0*mX*mX/m2)>0)
         rho = ciR*sqrt(1.0-4.0*mX*mX/m2);
    else rho = ciM*sqrt(4.0*mX*mX/m2-1.0);
    return rho;
}
inline double EvtDTopipienu::getWidthrho(double m, double m0, double width0, double p, double p0 )
{
    double widthRho = 0.0;
    widthRho = width0*pow(p/p0, 3)*(m0/m);
    return widthRho;
}
