#include <math.h>


double eps,stop,step,hmin,t0;
double Rg,Na,sigma,mu,rol,n,Pinf,cl,K,ap,alfae,Tinf,Rv,alfaM,epsilon,sigmap,e0,L,cpl,kl,alfa,P0,Patm;
double Amp,f,w,fi,Amp2,frec2,w2,fi2;
double a_h2o,a_Ar,b_h2o,b_Ar,a_Arh2o,b_Arh2o;
double Cvh2o,CvAr,CvO2,CvN2;
double a1,a2,a3,a4;
double Mh2o,MAr,MO2,MN2,MH2,MH,MO,MOH,MN,MNH,MNH2,MNH3,MN2H,MNO,MNO2,MN2O,MHNO;
double k1H2,k1H,k1O,k1O2,k1OH,k1h2o,k1N,k1NH,k1NH2,k1NH3,k1N2H,k1NO,k1NO2,k1N2O,k1HNO,k1N2,k1Ar;
double k2H2,k2H,k2O,k2O2,k2OH,k2h2o,k2N,k2NH,k2NH2,k2NH3,k2N2H,k2NO,k2NO2,k2N2O,k2HNO,k2N2,k2Ar;
int N,NumEsp;
double Vi,Rh2o,Ri,R0i;
double Tmin, Tmax,T0;
double mpmin, mpmax;
double Cmin,Cmax;
double ngas0,nh2o0,nh20,nh0,no0,no20,noh0,nn0,nnh0,nnh20,nnh30,nn2h0,nno0,nno20,nn2o0,nhno0,nn20,nar0,nh2o20,nho20;
double nt0Ar,RAr;
double Dif,DN2,DO2,DAr,C0N2,C0O2,C0Ar;
int tcal;
int rapidez, yndisoc,yndif;
double beta, pvap0;
int modelo;
int NCh,NCh2,NCh3,Nint,Npoints,Nvar,Nsalida;
double gama;
int numciclos;
double Ab1,Ab2,Ab3,Ab4,Ab5,Ab6,Cb1,Cb2,Cb3,Cb4,Cb5,Cb6,betab1,betab2,betab3,betab4,betab5,betab6;
int Ndatos,NdatosN;
double tmaxexp,tminexp;
double Rmax;
int nn;
double nu;
double Intensidad,duracion;
long double hPlanck,clight;
int RayTay,iespectro;
long double lambdapuntos,TempCorte,CTR,DTR,EION,E2,Z;
double t1,t2,t5,t6,t15,t16,t33,t34,t39;
double Ppert,Rpert,Vpert;



void parametros(){

   //Universal Constants
   Rg=8.31451135945; //constante de los gases [J/gmol K]
   Na=6.0221367E23; //Numero de Avogadro [#/mol]
   K= 1.38065805307E-23; //constante de Boltzmann []

   //constantes
   n=5.0; //constante [1] Yasui
   ap=0.827; //constante [1] Yasui
   sigmap=0.4e-18; //[] Yasui
   e0=1.0e-3; //[] no usadoen ningun lado


   /*Liquid properties(agua) */
   rol=998.00399; //densidad del liquido [kg/m3]
   mu=0.000959;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4181.0;//calor especifico del liquido [J/kg K]
   kl=0.606;//conductividad termica del agua [W/m K] NO USADO
   sigma=0.0727;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
   nu=mu/rol;
   T0=293.15; //Initial temperature [K]
   Tinf=293.15; //Liquid temperature far from bubble [K]
   P0=1.013e5;
   Pinf=0.925E5;  //presion en el infinito [Pa]
   Patm=1.013e5;
   pvap0=133.322*exp(18.3036-3816.44/(T0-46.13));

/*
   //Liquid properties (85% H2SO4)
   rol=1778.6; //densidad del liquido [kg/m3]
   mu=0.015;  //viscosidad [Pa s]
   cl=1473.0; //velocidad del sonido del liquido [m/s]
   cpl=1829.282;//calor especifico del liquido [J/kg K]
   kl=0.35;//conductividad termica del agua [W/m K]
   sigma=0.0712;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
   nu=mu/rol;
   T0=293.15; //Initial temperature [K]
   Tinf=293.15; //Liquid temperature far from bubble [K]
   Pinf=0.92E5;  //presion en el infinito [Pa]
   P0=1.013E5;  //presion de una atmosfera [Pa]
   Patm=1.013E5;
   pvap0=100.0;
*/

  /*
   //Liquid properties (98% H2SO4)
   rol=1831.0; //densidad del liquido [kg/m3]
   mu=0.025;  //viscosidad [Pa s]
   cl=1470.0; //velocidad del sonido del liquido [m/s]
   cpl=1382.0;//calor especifico del liquido [J/kg K]
   kl=0.314;//conductividad termica del agua [W/m K]
   sigma=0.0712;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
   nu=mu/rol;
   T0=293.15; //Initial temperature [K]
   Tinf=293.15; //Liquid temperature far from bubble [K]
   Pinf=0.2E5;  //presion en el infinito [Pa]
   P0=1.013E5;  //presion de una atmosfera [Pa]
   pvap0=0.15;
   Patm=1.013e5;
 */


   //SL Parameters
   Intensidad=0.; //Potencia
   duracion=2.0e-8; //Duracion del pulso
   Amp=0.0e5;  //amplitud de la onda acustica [Pa]
   f=28300.0; //frecuencia de la onda acustica
   w=2.0*3.14159*f;  //frecuencia angular de la onda acustica [1/s]
   fi=3.14159; //fase de la onda [rad]
   Amp2=0.0*20./250.0*Amp;
   frec2=11.0*f;
   w2=2.0*3.14159*frec2;
   fi2=3.14159;

//
//D2O
   //Vapor Properties
   a_h2o=0.5583; //constante de Van der Waals del agua [m6 Pa/mol]
   b_h2o=3.0877e-5; //constante de Van der Waals del agua [m3/mol]
   Mh2o=0.020; //masa molar del gas [kg/gmol]
   Cvh2o=3.0*Rg; //calor especifico a volumen constante del vapor de agua
   k1h2o=-0.02028; //constantes de la conductividad del vapor [W/m K]
   k2h2o=1.15006e-4;
   Rv=Rg/Mh2o; //constante de los gases para el vapor de agua [J/kg K]
   beta=1.0;
//
//
/*H2O */
/*   a_h2o=0.5536; //constante de Van der Waals del agua [m6 Pa/mol]
   b_h2o=3.049e-5; //constante de Van der Waals del agua [m3/mol]
   Mh2o=0.018; //masa molar del gas [kg/gmol]
   Cvh2o=3.0*Rg; //calor especifico a volumen constante del vapor de agua
   k1h2o=-0.01332; //constantes de la conductividad del vapor [W/m K]
   k2h2o=1.00655e-4;
   Rv=Rg/Mh2o; //constante de los gases para el vapor de agua [J/kg K]
   beta=1.0;
*/

   //Gas Properties - Argon
   a_Ar=0.1345; //constante de Van der Waals del argon []
   b_Ar=3.219e-5;//constante de Van der Waals del argon []
   MAr=0.03942; //masa molar del argon [kg/gmol]
   CvAr=3.0/2.0*Rg; //heat capacity at constant volume gas[J/kg K]
   k1Ar=0.009; //constantes de la conductividad del gas [W/m K]
   k2Ar=3.2e-5;
   gama=Rg/(3.0*Rg)+1.0;

   //Gas Properties - Helio
/*   a_Ar=0.00346; //constante de Van der Waals del argon []
   b_Ar=2.38e-5;//constante de Van der Waals del argon []
   MAr=0.004002602; //masa molar del argon [kg/gmol]
   CvAr=3.0/2.0*Rg; //heat capacity at constant volume gas[J/kg K]
   k1Ar=0.072438133; //constantes de la conductividad del gas [W/m K]
   k2Ar=0.000271733;
   gama=Rg/(3.0*Rg)+1.0;
*/
   //Gas Properties - Oxigeno
   //a_Ar=0.136; //constante de Van der Waals del argon []
   //b_Ar=3.183e-5;//constante de Van der Waals del argon []
   MO2=0.032; //masa molar del argon [kg/gmol]
   CvO2=3.0/2.0*Rg; //heat capacity at constant volume gas[J/kg K]
   k1O2=0.009; //constantes de la conductividad del gas [W/m K]
   k2O2=3.2e-5;

   //Gas Properties - Nitrogeno
   //a_Ar=0.139; //constante de Van der Waals del argon []
   //b_Ar=3.913e-5;//constante de Van der Waals del argon []
   MN2=0.028; //masa molar del argon [kg/gmol]
   CvN2=5.0/2.0*Rg; //heat capacity at constant volume gas[J/kg K]
   k1N2=669e-5; //constantes de la conductividad del gas [W/m K]
   k2N2=6.24357e-5;

   //Liquid-Vapor properties
//   L=2258987.94;// H2O calor latente de vaporizacion [J/kg]
   L=2066700.0;// D2O calor latente de vaporizacion [J/kg]
   alfae=1.0; //coeficiente de acomodacion termico [1]
   alfaM=0.1; //[1]//coeficiente de acomodacion


   //Temporal integration control
   eps=1.0e-7;//error
   step=1.0e-11; //Time step [s]
   hmin=1.0e-16; //el paso minimo que debe tener en cuenta
   numciclos=1; //numero de ciclos que calcula NO USADO
   epsilon=1.0e-6; //Convergence criterion para las cantidades que se calculan con el metodo de biseccion[1]
   N=100; // cantidad de pasos utilizados para integrar la funcion error (ERF())

   //Initial conditions (t=0)
   Ri=425.e-6; //radio maximo que adquiere
   Vi=0.0; //velocidad inicial [m/s]
   R0i=24.080e-6; //radio de gases no condensables iniciales (O2 H2)


   //Computed Values
   a_Arh2o=sqrt(a_Ar*a_h2o); //VW gas-vapor
   b_Arh2o=pow(0.5*(pow(b_Ar,1.0/3.0)+pow(b_h2o,1.0/3.0)),3.0); //VW gas-vapor

   //Limits
   Tmin=200.0; //Minimum temperature para el calculo de biseccion [K]
   Tmax=20000.0; //Maximum temperature para el calculo de biseccion [K]
   mpmin=-5.5e4; //Minimum mp para el calculo de biseccion [kg/m2 s]
   mpmax=5.5e4; //Maximum mp para el calculo de biseccion [kg/m2 s]
   Cmin=0.0;
   Cmax=4.0;
   Rmax=1.0e-6;


   //coeficientes de difusion
   Dif=2.39e-5; //dentro de la burbuja-vapor (h20) en gas (ar y demas) NO USADO
   DN2=2.6e-9;  //no usado
   DO2=2.4e-9;  //no usado
   DAr=2.5e-9;  //Ar en H2O a 25 grados C


   //concentraciones de saturacion a p0 y t0
   C0N2=0.0; // NO USADO
   C0O2=0.0; // NO USADO
   //C0H2=0.0; // NO USADO
   C0Ar=0.061; //BUSCAR CONCENTRACION DE SATURACION DE ARGON EN MEZCLA DE SULFURICO NO USADO


   //constantes de disociacion para el oxigeno y el nitrogeno
   Ab1=1.8125e12; //no usado
   Ab2=1.0695e13; //no usado
   Ab3=3.3917e12; //no usado
   Ab4=1.9024e11; //no usado
   Ab5=0.0e11;  //no usado
   Ab6=4.8027e11; //no usado
   Cb1=59419.3; //no usado
   Cb2=59419.3; //no usado
   Cb3=59419.3; //no usado
   Cb4=89000.0; //no usado
   Cb5=89000.0; //no usado
   Cb6=89000.0; //no usado
   betab1=-1.0; //no usado
   betab2=-1.0; //no usado
   betab3=-1.0; //no usado
   betab4=-0.5; //no usado
   betab5=-0.5; //no usado
   betab6=-0.5; //no usado


   //switches del modelo

   tcal=1; //switch utilizado para poner o sacar transferencia de calor.
           // Si es Cero no hay, otro numero y hay TC
   rapidez=1; //switch utilizado para tomar a Gama como una funcion de mp, o como constante igual
	      //a uno. Si es 1 es rapido, otro numero y es lento
   modelo=2;  //if modelo = 1 -> Yasui, modelo = 2 -> Toegel
	      //if modelo = 3 -> Hao-Prosperetti, modelo = 4 -> Hibrido entre H-P y Toegel
   yndisoc=1;
   yndif=1;


   //Numeros
   NCh=-1; //Cantidad de Polinomios de Chebyshev. Para usar los modelos de Yasui y Toegel
           //se debe poner NCh=-1.
   NCh2=NCh+1;//cantidad de terminos totales (+1 porque esta el termino cero)
   NCh3=NCh+6+1; //cantidad de variables a calcular (las usuales + los polinomios de Chebyshev)
   NumEsp=19; //Numero de especies totales
   Nvar=NCh3+NumEsp+2;


   //Inestabilidades
   nn=2; //orden de inestabilidad de SI
   Ppert=1.0; //Amplitud de deformacion (Chinos)
   Rpert=0.0; //no usada
   Vpert=0.0; //no usada

   //Masas Molares D2O
   MH2=0.004;
   MH=0.002;
   MO=0.016;
   MOH=0.018;
   MN=0.014;
   MNH=0.015;
   MNH2=0.016;
   MNH3=0.017;
   MN2H=0.029;
   MNO=0.03;
   MNO2=0.046;
   MN2O=0.044;
   MHNO=0.031;

//
   /*Masas Molares H20*/
/*
   MH2=0.002;
   MH=0.001;
   MO=0.016;
   MOH=0.017;
   MN=0.014;
   MNH=0.015;
   MNH2=0.016;
   MNH3=0.017;
   MN2H=0.029;
   MNO=0.03;
   MNO2=0.046;
   MN2O=0.044;
   MHNO=0.031;
*/

   //conductividades
   k1H2=k1O2; //ordenadas al origen
   k1H=k1O2;
   k1O=k1O2;
   k1OH=k1O2;
   k1N=k1O2;
   k1NH=k1O2;
   k1NH2=k1O2;
   k1NH3=k1O2;
   k1N2H=k1O2;
   k1NO=k1O2;
   k1NO2=k1O2;
   k1N2O=k1O2;
   k1HNO=k1O2;
   k1O2=k1O2;

   k2H2=k2O2; //pendientes (respecto de la temperatura: k1i+k2i*T)
   k2H=k2O2;
   k2O=k2O2;
   k2OH=k2O2;
   k2N=k2O2;
   k2NH=k2O2;
   k2NH2=k2O2;
   k2NH3=k2O2;
   k2N2H=k2O2;
   k2NO=k2O2;
   k2NO2=k2O2;
   k2N2O=k2O2;
   k2HNO=k2O2;

   //Valores relacionados con el espectro
   iespectro=6;
   TempCorte=1.;
   lambdapuntos=60.0;
   Z=1;
   //aAr=0.106e-9; %effective atomic diameter Ar (covalente)
   //aH=0.031e-9; %effective atomic diameter H (covalente)

   //DTR=-0.6e-20*pow((0.031e-9/0.106e-9),2);
   //CTR=1.6e-20/1.6e-19*pow((0.031e-9/0.106e-9),2);

   DTR=-0.6e-20;
   CTR=1.6e-20/1.6e-19;
   hPlanck=6.626e-34;
   clight=2.998e8;
   EION=23.8*1.6e-19;
   E2=19.75*1.6e-19;
//   EION=13.59843449*1.6e-19;
//   E2=10.19880606*1.6e-19;

   //Parametros de disociacion
   t1=1.0;
   t2=1.0;
   t5=1.0;
   t6=1.0;
   t15=-1.0;
   t16=1.0;
   t33=-1.0;
   t34=-1.0;
   t39=1.0;

   /*Valores para comparacion con experimentos
   tminexp=-9.39E-6;
   tmaxexp=4.06E-5;
   Ndatos=5000;
   NdatosN=7000;
   */



   //Parametros del agua a distintas temperaturas

/*
   Tinf=275;

   rol=1000.0; //densidad del liquido [kg/m3]
   mu=0.001652;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4211.0;//calor especifico del liquido [J/kg K]
   kl=0.574;//conductividad termica del agua [W/m K]
   sigma=0.0753;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=280;

   rol=1000.0; //densidad del liquido [kg/m3]
   mu=0.001422;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4198.0;//calor especifico del liquido [J/kg K]
   kl=0.582;//conductividad termica del agua [W/m K]
   sigma=0.0748;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=285;

   rol=1000.0; //densidad del liquido [kg/m3]
   mu=0.001225;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4189.0;//calor especifico del liquido [J/kg K]
   kl=0.590;//conductividad termica del agua [W/m K]
   sigma=0.0743;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=290;

   rol=999.0009991; //densidad del liquido [kg/m3]
   mu=0.001080;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4184.0;//calor especifico del liquido [J/kg K]
   kl=0.598;//conductividad termica del agua [W/m K]
   sigma=0.0737;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=295;

   rol=998.00399; //densidad del liquido [kg/m3]
   mu=0.000959;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4181.0;//calor especifico del liquido [J/kg K]
   kl=0.606;//conductividad termica del agua [W/m K]
   sigma=0.0727;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=300;

   rol=997.00897; //densidad del liquido [kg/m3]
   mu=0.000855;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4179.0;//calor especifico del liquido [J/kg K]
   kl=0.613;//conductividad termica del agua [W/m K]
   sigma=0.0717;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=305; //32 grados

   rol=995.0248; //densidad del liquido [kg/m3]
   mu=0.000769;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4178.0;//calor especifico del liquido [J/kg K]
   kl=0.620;//conductividad termica del agua [W/m K]
   sigma=0.0709;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=310;

   rol=993.048659; //densidad del liquido [kg/m3]
   mu=0.000695;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4178.0;//calor especifico del liquido [J/kg K]
   kl=0.628;//conductividad termica del agua [W/m K]
   sigma=0.0700;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
*/
/*
   Tinf=315;

   rol=991.08; //densidad del liquido [kg/m3]
   mu=0.000631;  //viscosidad [Pa s]
   cl=1480.0; //velocidad del sonido del liquido [m/s]
   cpl=4179.0;//calor especifico del liquido [J/kg K]
   kl=0.634;//conductividad termica del agua [W/m K]
   sigma=0.0692;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;

    pvap0=133.322*exp(18.3036-3816.44/(Tinf-46.13));
*/
//
   Tinf=4.8+273.15; //D2O

   rol=1105.6; //densidad del liquido [kg/m3]
   mu=0.0020713;  //viscosidad [Pa s]
   cl=1324.9; //velocidad del sonido del liquido [m/s]
   cpl=4211.3;//calor especifico del liquido [J/kg K]
   kl=0.56492;//conductividad termica del agua [W/m K]
   sigma=0.074903;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
   pvap0=669.81;
//
//
   Tinf=35+273.15; //D2O

   rol=1101.8; //densidad del liquido [kg/m3]
   mu=0.000869;  //viscosidad [Pa s]
   cl=1415.1; //velocidad del sonido del liquido [m/s]
   cpl=4236.50;//calor especifico del liquido [J/kg K]
   kl=0.60562;//conductividad termica del agua [W/m K]
   sigma=0.07031;  //tension superficial [J/m2]
   alfa=kl/cpl/rol;
   pvap0=4949;
//
//
   alfa=kl/cpl/rol;
   nu=mu/rol;
   T0=293.15; //Initial temperature [K]
   P0=1.013e5;
   Pinf=0.925E5;  //presion en el infinito [Pa]
   Patm=1.013e5;
//
}
