/*********** using the following decay cards ****************
Decay J/psi
1 gamma eta_c  PHSP;
Enddecay
 
Decay eta_c
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
    EvtVector4R pjsi,petac,pLam,pLamb,pp,pm;
    pjsi = part->getP4();
    petac= part->getDaug(1)->getP4();  //etac momentum 
    pLam= part->getDaug(1)->getDaug(0)->getP4();  //Lambda momentum
    pLamb= part->getDaug(1)->getDaug(1)->getP4();  //Lambda momentum

    pp = part->getDaug(1)->getDaug(0)->getDaug(0)->getP4();  //p+ from Lambda 
    pm = part->getDaug(1)->getDaug(1)->getDaug(0)->getP4();  //anti-p- from anti-Lambda 

    EvtHelSys angles0(pjsi,petac);
    double theta0 = angles0.getHelAng(1);
    double phi0   = angles0.getHelAng(2);
 
   EvtVector4R yL =pLam.cross(pp);  //standard helicity system
   EvtVector4R yLb=pLamb.cross(pm);  //standard helicity system

   EvtVector4R yeL =petac.cross(pLam);
   EvtVector4R yeLb =petac.cross(pLamb);

   
    EvtmyEulerAngles angles2(yeL,pLam,yL,pp,"yz");
    double theta2 = angles2.gettheta();
    double phi2   = angles2.getphi();

    EvtmyEulerAngles angles3(yeLb,pLamb,yLb,pm,"yz");
    double theta3 = angles3.gettheta(); 
    double phi3   = angles3.getphi();

//std::cout<<angles2.getpsi()<<" "<<angles3.getpsi()<<std::endl;

    double aL=0.75;
    double aLb=-aL;

    amps=(3+cos(2*theta0)) * (1+ aL*aLb*cos(theta2)*cos(theta3) - aL*aLb* sin(theta2)*sin(theta3) * cos(phi2+phi3) );

