//---- read out the momentum for each particle, then using them to calculate the amplitude squared, for example:
/*********** using the following decay cards ****************
Decay J/psi 
1 Lambda0 anti-Lambda0 PHSP;
Enddecay
 
Decay Lambda0
1 p+ pi- PHSP;
Enddecay
 
Decay anti-Lambda0
1 anti-p- pi+ PHSP;
Enddecay

End
*////////
//----------------------------------
  EvtVector4R vpho,lambda,lambdabar,proton,antiproton;
  vpho      = part->getP4();              //Jpsi momentum  
  lambda   = part->getDaug(0)->getP4();  //lambda momentum 
  lambdabar= part->getDaug(1)->getP4();  //anti-Lambda  momentum
  proton    = part->getDaug(0)->getDaug(0)->getP4();  //proton momentum
  antiproton= part->getDaug(1)->getDaug(0)->getP4();  //anti-proton momentum

  EvtHelSys angles0(vpho,lambda);
  double theta = angles0.getHelAng(1);
  double phi   = angles0.getHelAng(2);
 
  EvtHelSys angles1(lambda,proton);
  double theta1 = angles1.getHelAng(1);
  double phi1   = angles1.getHelAng(2);

  EvtHelSys angles2(lambdabar,antiproton);
  double theta2 = angles2.getHelAng(1);
  double phi2   = angles2.getHelAng(2);
  
  m_vangs.clear();
  double  Delta=0.74;
  double  alphaL =0.75;
  double  alphaLb=-0.75;
  double  alpha   =0.461;
  m_vangs.push_back(theta);
  m_vangs.push_back(phi);
  m_vangs.push_back(theta1);
  m_vangs.push_back(phi1);
  m_vangs.push_back(theta2);	
  m_vangs.push_back(phi2);
   amps=
     (2*(2 + alpha + alpha*cos(2*theta) - sqrt(1 - pow(alpha,2))*alphaLb*sin(Delta)*sin(phi2)*sin(2*theta)*sin(theta2)) +
      alphaL*(2*sqrt(1 - pow(alpha,2))*sin(Delta)*sin(phi1)*sin(2*theta)*sin(theta1) +
               alphaLb*(2*sin(theta1)*(sqrt(1 - pow(alpha,2))*cos(Delta)*cos(phi1)*cos(theta2)*sin(2*theta) +
                                        2*(cos(phi1)*cos(phi2) + alpha*sin(phi1)*sin(phi2))*pow(sin(theta),2)*sin(theta2)) -
                         2*cos(theta1)*((1 + 2*alpha + cos(2*theta))*cos(theta2) +
                                        sqrt(1 - pow(alpha,2))*cos(Delta)*cos(phi2)*sin(2*theta)*sin(theta2)))))/2.;

  
