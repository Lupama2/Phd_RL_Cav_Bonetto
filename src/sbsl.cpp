//---------------------------------------------------------------------------

#pragma hdrstop

//---------------------------------------------------------------------------

#pragma argsused


#include <vcl.h>
#include<math.h>
#include <iostream.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream.h>
#include <iomanip.h>


//***************INICIO DE LLAMADO DE FUNCIONES*********************************
void derivada(double t,double *y,double TB,double mp,double k,double *colpoints,
              double **Q, double **R, double *dydt);


void rk4(double *y,double *dydx,double x,double h, double TB, double mp,double k,
         double *colpoints, double **Q, double **R,double *yout);



void rkqc(double *y,double *dydx,double *x,double htry,double eps,double *yscal,
          double *hdid,double *hnext,double TB, double mp, double k,double *colpoints,
          double **Q,double **R, double *yout);

void parametros();

double min(double arg1, double arg2);
double fact (int n);
double erf(double omega,int N);
double calculamp(double mpmin,double mpmax,double epsilon, double pvap, double pv,
                 double T,double TLi,double TB,int N);
double calculaTBT(double Tmin,double Tmax,double epsilon,double TLi,double k,
                  double lth,double np,double T,double m);
double calculaTBY(double Tmin,double Tmax,double epsilon,double TLi,double k,
                  double landa,double np,double T,double m);
double calculaTBHP(double Tmin,double Tmax,double epsilon,double TLi,double qpp,
                   double np,double m);
double calculaCinfCsat2(double Cmin,double Cmax,double epsilon,double f1,
                        double f2,double f3,double f4,double f5,double f6,
                        double f7,double f8, double f9,double nt,double porcgn,
                        double porcn2,double porco2);
double PolCheb(double y, int NCh);
double DerPCheb(double y, int NCh);
double DerSCheb(double y, int NCh);

void Conductividad1(double *k1);
void Conductividad2(double *k2);
void MasaMolar(double *Masa);

long double espectro(long double RADIUS,long double TEMPERATURE, long double NAR,
                     long double Z,long double LAMBDA,long double EION,
                     long double CTR, long double DTR, long double E2);


//*************FIN DE LLAMADO DE FUNCIONES**************************************



//*************INICIO DE LA PARTE PRINCIPAL*************************************

void main(){

   char *nombre1="res1.dat";
   char *nombre2="res2.dat";
   char *nombre3="res3.dat";
   char *nombre4="res4.dat";
   char *nombre5="res5.dat";
   char *nombre6="res6.dat";
   char *nombre7="res7.dat";


   ofstream out1(nombre1);
   if(!out1){
      cout<<"no pude abrir el archivo de salida 1"<<endl;
      exit(-1);
   }

   ofstream out2(nombre2);
   if(!out2){
      cout<<"no pude abrir el archivo de salida 2"<<endl;
      exit(-1);
   }

   ofstream out3(nombre3);
   if(!out3){
      cout<<"no pude abrir el archivo de salida 3"<<endl;
      exit(-1);
   }

   ofstream out4(nombre4);
   if(!out4){
      cout<<"no pude abrir el archivo de salida 4"<<endl;
      exit(-1);
   }

   ofstream out5(nombre5);
   if(!out5){
      cout<<"no pude abrir el archivo de salida 5"<<endl;
      exit(-1);
   }

   ofstream out6(nombre6);
   if(!out6){
      cout<<"no pude abrir el archivo de salida 6"<<endl;
      exit(-1);
   }

   ofstream out7(nombre7);
   if(!out7){
      cout<<"no pude abrir el archivo de salida 7"<<endl;
      exit(-1);
   }


   double y3old;
   y3old=0.0;

   int i,j;
   int kk;
   int cont;
   double sum;


   parametros();


   extern double eps,stop,step;
   extern double Ri,R0i,Vi,T0,nt0h2o,nt0Ar,RAr,P0;
   extern double ngas0,nh2o0,nh20,nh0,no0,no20,noh0,nn0,nnh0,nnh20,nnh30,nn2h0,nno0,nno20,nn2o0,nhno0,nn20,nar0,nh2o20,nho20;
   extern double ap,alfaM,alfae,e0;
   extern double sigmap,Na,Rg,Rv,n,kl,cpl,rol,Tinf,L,sigma,mu,nu,K,alfa;
   extern double Mh2o,MAr,MO2,MN2,MH2,MH,MO,MOH,MN,MNH,MNH2,MNH3,MN2H,MNO,MNO2,MN2O,MHNO;
   extern double w,fi,Amp,cl,f,w2,fi2,Amp2,frec2;
   extern double gama, Rh2o, avw,Pinf;
   extern double Cvh2o,CvAr,CvO2,CvN2;
   extern double Tc,Pc;
   extern double Tmin,Tmax;
   extern double epsilon;
   extern double mpmin,mpmax;
   extern double a_h2o,a_Ar,b_h2o,b_Ar,a_Arh2o,b_Arh2o;
   extern double a1,a2,a3,a4;
   extern double k1H2,k1H,k1O,k1O2,k1OH,k1h2o,k1N,k1NH,k1NH2,k1NH3,k1N2H,k1NO,k1NO2,k1N2O,k1HNO,k1N2,k1Ar;
   extern double k2H2,k2H,k2O,k2O2,k2OH,k2h2o,k2N,k2NH,k2NH2,k2NH3,k2N2H,k2NO,k2NO2,k2N2O,k2HNO,k2N2,k2Ar;
   extern double Dif,DN2,DO2,DAr;
   extern int modelo;
   extern int rapidez;
   extern int NCh,NCh2,NCh3,Nvar,Nsalida;
   extern int N,Npoints,Ndatos,NumEsp;
   extern int numciclos,nn;
   extern double beta, pvap0;
   extern double Ab1,Ab2,Ab3,Ab4,Ab5,Ab6,Cb1,Cb2,Cb3,Cb4,Cb5,Cb6,betab1,betab2,betab3,betab4,betab5,betab6;
   extern double C0N2,C0O2,C0Ar;
   extern double Rmax;
   extern double tminexp,tmaxexp;
   extern double Cmin,Cmax;
   extern double hP,qe,ep0,landa2,Eion,me,ctr,dtr;
   extern double Intensidad,duracion;
   extern double nh2o0,ngas0;
   extern long double hPlanck,clight;
   extern int iespectro;
   extern long double lambdapuntos;
   extern long double TempCorte;
   extern int RayTay;
   extern int Npulso,ConPulso;
   extern long double CTR,DTR,E2,Z,EION;


   double *y=new double[Nvar];
   double *yout=new double[Nvar];
   double *dydx=new double[Nvar];
   double *yscal=new double[Nvar];
   double *Pl=new double[10];

   double *Tiempo=new double[Nsalida];
   double *Ann=new double[Nsalida];
   double *Bnn=new double[Nsalida];
   double *Radio=new double[Nsalida];


   int NumEspm1=NumEsp+1;

   double *Masa=new double[NumEspm1];
   double *k1=new double[NumEspm1];
   double *k2=new double[NumEspm1];

   double x,htry,hdid,hnext,xold;
   double a,b;


   double *aux=new double[10];
   double grad;
   double DPTCh;

   double TB;
   double ps,ps2,pg,pvap,pv,pgg;
   double mp;
   double v,V,np,nAr,ntot,ngas;
   double k;
   double landa,M,m,qpp,rog;

   double Xi,lth;
   int icont,nsal;
   double Tmaxx,Vmax,Amax,pmax;
   int Nvar2=NCh3-1;

   extern int Nsalida;
   double Rprueba,Rgas;
   double tfin;


   double *colpoints=new double[NCh2];


   double **Ch;
   Ch = new double*[NCh2];
   for (j=0; j<=NCh; j++)
      Ch[j] = new double[NCh2];

   double **vv;
   vv = new double*[NCh2];
   for (j=0; j<=NCh; j++)
      vv[j] = new double[NCh2];

   double **Q;
   Q = new double*[NCh2];
   for (j=0; j<=NCh; j++)
      Q[j] = new double[NCh2];


   double **R;
   R = new double*[NCh2];
   for (j=0; j<=NCh; j++)
      R[j] = new double[NCh2];

   double pi=3.141592653589;

   long double *electra=new long double[lambdapuntos+1];

   long double RADIUS,TEMPERATURE,NAR;
   long double LAMBDA;
   long double plt;
   long double lambdainicial;
   long double lambdafinal;
   long double deltalambda;
   long double fluencia,fotones,potencia;
   int iLAMBDA;



   //********************INICIALIZACION*****************************************

   //definicion de los puntos de colocacion, donde i=0 corresponde al centro de
   //la burbuja, y i=NCh es y=1, la interfase

   for(i=0;i<=NCh;i++){
      colpoints[i]=cos(pi*double(NCh-i)/2.0/double(NCh));
   }

   //matriz que define el sistema de ecuaciones para los coeficientes que
   //acompa�an a los polinomios de Chebyshev

   for(j=0;j<=NCh;j++){
      for(i=0;i<=NCh;i++)
         Ch[i][j]=PolCheb(colpoints[i],2*j);
   }

   //***************************************************************************
   //************************Descomposicion QR**********************************
   //***************************************************************************


   //inicializacion de la matriz Q

   for(i=0;i<=NCh;i++){
      for(j=0;j<=NCh;j++)
         Q[i][j]=0.0;
   }

   //inicializacion de la matriz R

   for(i=0;i<=NCh;i++){
      for(j=0;j<=NCh;j++)
         R[i][j]=0.0;
   }

   if(NCh2>1){

      cont=0;
      sum=0.0;
      while(cont<=NCh){
         sum=sum+Ch[cont][0]*Ch[cont][0];
         cont=cont+1;
      }

      for(kk=0;kk<=NCh;kk++)
         Q[kk][0]=Ch[kk][0]/sqrt(sum);

      R[0][0]=sqrt(sum);



      for(j=1;j<=NCh;j++){

         for(kk=0;kk<=NCh;kk++)
            vv[kk][j]=Ch[kk][j];

         for(i=0;i<=(j-1);i++){

            cont=0;
            while(cont<=NCh){
               R[i][j]=R[i][j]+Q[cont][i]*Ch[cont][j];
               cont=cont+1;
            }

            for(kk=0;kk<=NCh;kk++)
               vv[kk][j]=vv[kk][j]-R[i][j]*Q[kk][i];

         }//fi del for i

         cont=0;
         while(cont<=NCh){
            R[j][j]=R[j][j]+vv[cont][j]*vv[cont][j];
            cont=cont+1;
            }

         R[j][j]=sqrt(R[j][j]);

         for(kk=0;kk<=NCh;kk++)
            Q[kk][j]=vv[kk][j]/R[j][j];

      }//fin del for j
   }


   Rmax=1.0e-6;
   Tmaxx=Tinf;
   Vmax=0.0;
   Amax=0.0;

   //***************************************************************************
   //**********BLOQUE DONDE SE CAMBIAN LOS PARAMETROS***************************
   //***************************************************************************
   //FABIAN

   //alfaM=0.1;
   //Amp=0.0e5;
   //Ri=400.0e-6; //radio maximo que adquier
   //RAr=30.0e-6; //radio de gases no condensables iniciales (O2 H2)
   //Rgas=0.0e-6; nitrogeno
   //Vi=0.;
   //Tinf=293.0;
   //T0=293.0;
   //Intensidad=0.0; //siempre igual a cero en este nodo


   kk=0;

   x=0.0; //Inicializacion del tiempo
   y[0]=0.0;
   y[1]=Ri; //Inicializacion del radio
   y[2]=Vi; //Inicializacion de la velocidad
   y[3]=Tinf; // Inicializacion de la temperatura
   nh2o0=Na*pvap0*(4.0/3.0*pi*y[1]*y[1]*y[1])/Rg/(Tinf); //beta*Na*pvap0*(4.0/3.0*3.14159*pow(Rprueba,3.0))/Rg/T0;
   y[4]=nh2o0;//cantidad de particulas iniciales de vapor de agua
   y[5]=Tinf; // Inicializacion de la temperatura de la pared de la burbuja en la cara externa
   //RAr=pow(pow(Ri,3.0)-pow(Rgas,3.0),1./3.0);
   nAr=Na*P0*4.0/3.0*3.14159*pow(R0i,3.0)/Rg/T0;


    if(NCh3>6){
      for(i=6;i<NCh3;i++)
         y[i]=0.0;
   }

//   nh2o0=pvap0*(4.0/3.0*pi*pow(Ri,3.0))*Na/Rg/(Tinf);

   nh20=y[Nvar2+1]=2.0/3.0*nAr;
//   nh20=y[Nvar2+1]=0.*nAr;
   nh0=y[Nvar2+2]=0.0*nAr;
   no0=y[Nvar2+3]=0.0/3.0*nAr;
   no20=y[Nvar2+4]=1.0/3.0*nAr;
//   no20=y[Nvar2+4]=0.0/3.0*nAr;
   noh0=y[Nvar2+5]=0.0;
   nh2o0=y[Nvar2+6]=nh2o0;
//   nh2o0=y[Nvar2+6]=0.;
   nn0=y[Nvar2+7]=0.0*nAr;
   nnh0=y[Nvar2+8]=0.0;
   nnh20=y[Nvar2+9]=0.0;
   nnh30=y[Nvar2+10]=0.0;
   nn2h0=y[Nvar2+11]=0.0;
   nno0=y[Nvar2+12]=0.0*nAr;
   nno20=y[Nvar2+13]=0.0;
   nn2o0=y[Nvar2+14]=0.0;
   nhno0=y[Nvar2+15]=0.0;
   nn20=y[Nvar2+16]=0.0*nAr;
//   nar0=y[Nvar2+17]=0.0*nAr;
   nar0=y[Nvar2+17]=0.0*nAr;
   nh2o20=y[Nvar2+18]=0.0*nAr;
   nho20=y[Nvar2+19]=0.0*nAr;

   //ngas0=Na*P0*4.0/3.0*3.14159*pow(Rgas,3.0)/Rg/T0; //nitrogeno

   y[NCh3+NumEsp]=Ri;
   y[NCh3+NumEsp+1]=0.0;

   //***************************************************************************
   //***************************FIN DEL BLOQUE**********************************
   //***************************************************************************


   for(i=0;i<Nvar;i++)
      yscal[i]=y[1];


   for(i=0;i<10;i++)
      Pl[i]=0.0;

   htry=step; //Inicializacion del paso optativo
   hdid=step; // Inicializacion del paso anterior
   hnext=step; // Inicializacion del paso siguiente

   fluencia=0.0;
   fotones=0.0;


   for(i=0;i<Nsalida;i++){
      Radio[i]=0.0;
      Ann[i]=0.0;
      Bnn[i]=0.0;
   }

   tfin=700.0e-6; //cuanto tiempo calcula

   out6<<"tiempo"<<" "<<"radio"<<" "<<"velocidad"<<" "<<"aceleracion"<<" "<<"temperatura"<<" "<<
   "pg"<<" "<<"nH2"<<" "<<"nH"<<" "<<"nO"<<" "<<"nO2"<<" "<<"nOH"<<" "<<"nH2O"<<" "<<"nN"<<" "<<"nNH"
   <<" "<<"nNH2"<<" "<<"nNH3"<<" "<<"nN2H"<<" "<<"nNO"<<" "<<"nNO2"<<" "<<"nN2O"<<" "<<"nHNO"<<" "<<"nN2"
   <<" "<<"nAr"<<"nH2O2"<<" "<<"nHO2"
   <<" "<<endl;

   out7<<"Rmax"<<" "<<"Ri"<<" "<<"tiempo"<<" "<<"radio"<<" "<<"velocidad"<<" "<<"aceleracion"<<" "<<"temperatura"<<" "<<
   "pg"<<" "<<"nH2"<<" "<<"nH"<<" "<<"nO"<<" "<<"nO2"<<" "<<"nOH"<<" "<<"nH2O"<<" "<<"nN"<<" "<<"nNH"
   <<" "<<"nNH2"<<" "<<"nNH3"<<" "<<"nN2H"<<" "<<"nNO"<<" "<<"nNO2"<<" "<<"nN2O"<<" "<<"nHNO"<<" "<<"nN2"
   <<" "<<"nAr"<<"nH2O2"<<" "<<"nHO2"
   <<" "<<endl;
   

   nsal=0;

   while ((x<=tfin)&& (!kbhit())){


      xold=x;

      ps=Amp*sin(w*(x+y[1]/cl)-fi);

      ps2=Amp2*sin(w2*(x+y[1]/cl)-fi2);


      //************************************************************************
      //*********************************calcula el mp**************************
      //************************************************************************

      V=4.0/3.0*pi*pow(y[1],3.0); //volumen total

      MasaMolar(Masa);

      M=0.0;

      icont=1;

      while(icont<=NumEsp){
         M=M+(Masa[icont]*y[Nvar2+icont])/Na;
         icont=icont+1;
      }


      if(iespectro!=0){

         RADIUS=y[1]*1.0e6;
         TEMPERATURE=y[3]/10000.0;
//         NAR=y[Nvar2+17];
         NAR=y[Nvar2+1]+y[Nvar2+4];
         lambdainicial=0.2;
         lambdafinal=0.8;
         deltalambda=(lambdafinal-lambdainicial)/lambdapuntos;

         if (TEMPERATURE>TempCorte){

             potencia=0.0;

             for (iLAMBDA=0;iLAMBDA<=int(lambdapuntos);iLAMBDA++){

                LAMBDA=lambdainicial+iLAMBDA*deltalambda;
                plt=espectro(RADIUS,TEMPERATURE,NAR,Z,LAMBDA,EION,CTR,DTR,E2);

                if(plt<0)
                   plt=0.0;

                fluencia=fluencia+plt*(deltalambda*1e-6)*hdid;
                potencia=potencia+plt*(deltalambda*1e-6);
                fotones=fotones+plt*(deltalambda*1e-6)*hdid/(hPlanck*clight)*(LAMBDA*1.0e-6);
                electra[iLAMBDA]=electra[iLAMBDA]+plt*hdid;

//                if(iespectro==6)
//                   out4<<setprecision(12)<<xold<<" "<<(LAMBDA*1.0e-6)<<" "<<plt<<endl;

             }
             if(iespectro==6)
                out2<<setprecision(12)<<" "<<xold<<" "<<potencia<<endl;
         }
      }

      //************************************************************************
      //************************************************************************


      rog=M/V;

      ntot=0.0;
      icont=1;

      while(icont<=NumEsp){
         ntot=ntot+y[Nvar2+icont];
         icont=icont+1;
      }

      ngas=ntot-y[Nvar2+6];

      m=M/ntot; //masa media por particula

      landa=V/(sqrt(2.0)*sigmap*ntot);  //camino libre medio

      v=V/ntot*Na; //volumen molar

      np=ntot/V; //densidad de particulas

      k=0.0;

      Conductividad1(k1);
      Conductividad2(k2);

      icont=1;

      while(icont<=NumEsp){
         k=k+(k1[icont]+k2[icont]*y[3])*y[Nvar2+icont]/ntot;
         icont=icont+1;
      }

      np=ntot/V;

      Xi=k/(8.0/2.0*nh2o0*K+5.0/2.0*nAr*K);

      if(y[2] == 0)
         lth=y[1]/pi;
      else{
         if((y[1]/pi)<(sqrt(y[1]*Xi/fabs(y[2]))))
            lth=y[1]/pi;
         else
            lth=sqrt(y[1]*Xi/fabs(y[2]));
      }


      if(modelo == 1)
         TB= calculaTBY(Tmin,Tmax,epsilon,y[5],k,landa,np,y[3],m);

      if(modelo == 2)
         TB= calculaTBT(Tmin,Tmax,epsilon,y[5],k,lth,np,y[3],m);


      if(modelo == 3){
         grad=0.0;
         j=NCh;
         while(j>=0){
            DPTCh=DerPCheb(1,2*j);
            grad=grad+y[6+j]/y[1]*DPTCh;
            j=j-1;
         }

         TB= calculaTBHP(Tmin,Tmax,epsilon,y[5],grad,np,m);
      }


      if(modelo == 4){
         grad=0.0;
         j=NCh;
         while(j>=2){
            DPTCh=DerPCheb(1,2*j);
            grad=grad+y[6+j]/y[1]*DPTCh;
            j=j-1;
         }

         grad=grad+k*(Tinf-y[3])/lth;

         TB= calculaTBHP(Tmin,Tmax,epsilon,y[5],grad,np,m);

      }

      ngas=ntot-y[Nvar2+6];

      a = a_Ar*(ngas/ntot)*(ngas/ntot)+2.0*a_Arh2o*(ngas/ntot)*(y[Nvar2+6]/ntot)+
          a_h2o*(y[Nvar2+6]/ntot)*(y[Nvar2+6]/ntot);
          //const a de VW promediada con la cantidad de particulas

      b = b_Ar*(ngas/ntot)*(ngas/ntot)+2.0*b_Arh2o*(ngas/ntot)*(y[Nvar2+6]/ntot)+
          b_h2o*(y[Nvar2+6]/ntot)*(y[Nvar2+6]/ntot);
          //const b de VW promediada con la cantidad de particulas

      pg=Rg*y[3]/(v-b)-a/(v*v); //presion del gas

      pgg=ngas/ntot*pg;

      pvap=pvap0;

      pv=y[Nvar2+6]/ntot*pg; //presion parcial del vapor de agua


      if(rapidez==1)
         mp=alfaM/sqrt(2.0*3.14159*Rv)*(pvap/sqrt(y[5])-pv/sqrt(TB));

      else
         mp=calculamp(mpmin,mpmax,epsilon,pvap,pv,y[3],y[5],TB,N);


      //************************************************************************
      //************************************************************************
      //************************************************************************

      if((nsal%20) == 0){
         out6<<setprecision(12)<<xold<<" "<<y[1]<<" "<<y[2]<<" "<<dydx[2]<<" "<<y[3]<<" "<<pg<<" "<<
         y[Nvar2+1]<<" "<<y[Nvar2+2]<<" "<<
         y[Nvar2+3]<<" "<<y[Nvar2+4]<<" "<<y[Nvar2+5]<<" "<<y[Nvar2+6]<<" "<<
         y[Nvar2+7]<<" "<<y[Nvar2+8]<<" "<<
         y[Nvar2+9]<<" "<<y[Nvar2+10]<<" "<<y[Nvar2+11]<<" "<<y[Nvar2+12]<<" "<<
         y[Nvar2+13]<<" "<<y[Nvar2+14]<<" "<<
         y[Nvar2+15]<<" "<<y[Nvar2+16]<<" "<<y[Nvar2+17]<<" "<<
         y[Nvar2+18]<<" "<<y[Nvar2+19]<<" "<<endl;
      }

//       out7<<"Rmax"<<" "<<"Ri"<<" "<<"tiempo"<<" "<<"radio"<<" "<<"velocidad"<<" "<<"aceleracion"<<" "<<"temperatura"<<" "<<
//       "pg"<<" "<<"nH2"<<" "<<"nH"<<" "<<"nO"<<" "<<"nO2"<<" "<<"nOH"<<" "<<"nH2O"<<" "<<"nH2O2"<<" "<<
//       "nHO2"<<" "<<"R0"<<endl;

/*
        if(y[3] < y3old){
         out7<<setprecision(12)<<Ri<<" "<<R0i<<" "<<xold<<" "<<y[1]<<" "<<y[2]<<" "<<dydx[2]<<" "<<y[3]<<" "<<pg<<" "<<
         y[Nvar2+1]<<" "<<y[Nvar2+2]<<" "<<
         y[Nvar2+3]<<" "<<y[Nvar2+4]<<" "<<y[Nvar2+5]<<" "<<y[Nvar2+6]<<" "<<
         y[Nvar2+7]<<" "<<y[Nvar2+8]<<" "<<
         y[Nvar2+9]<<" "<<y[Nvar2+10]<<" "<<y[Nvar2+11]<<" "<<y[Nvar2+12]<<" "<<
         y[Nvar2+13]<<" "<<y[Nvar2+14]<<" "<<
         y[Nvar2+15]<<" "<<y[Nvar2+16]<<" "<<y[Nvar2+17]<<" "<<
         y[Nvar2+18]<<" "<<y[Nvar2+19]<<" "<<endl;
         y3old=0.0;
        }
        else
          y3old=y[3];
*/

      if(modelo == 1)
         qpp=k*(TB-y[3])/n/landa+mp*L+mp/Na*Cvh2o*(TB-y[5]);

      if(modelo == 2)
         qpp=k*(Tinf-y[3])/lth+mp*L+mp/Na*Cvh2o*(TB-y[5]);

      if(modelo == 3){
         j=NCh;
         qpp=0.0;

         while(j>=0){
            DPTCh=DerPCheb(1,2*j);
            qpp=qpp+y[6+j]/y[1]*DPTCh;
            j=j-1;
         }

         qpp=qpp+mp*L+mp/Na*Cvh2o*(TB-y[5]);

      }


      if(modelo == 4){
         j=NCh;
         qpp=0.0;

         while(j>=2){
            DPTCh=DerPCheb(1,2*j);
            qpp=qpp+y[6+j]/y[1]*DPTCh;
            j=j-1;
         }

         qpp=qpp+k*(Tinf-y[3])/lth+mp*L+mp/Na*Cvh2o*(TB-y[5]);

      }



      if((nsal%1000) == 0){
         cout<<xold<<" "<<y[1]<<" "<<y[3]<<" "<<y[Nvar2+17]<<endl;
      }

      nsal=nsal+1;


      derivada(x,y,TB,mp,k,colpoints,Q,R,dydx);

      rkqc(y,dydx,&x,htry,eps,yscal,&hdid,&hnext,TB,mp,k,colpoints,Q,R,yout);

      htry=hnext;

      for(i=1;i<Nvar;i++)
         y[i]=yout[i];

      for(i=1;i<6;i++)
         yscal[i]=y[i]+hdid*yscal[i];

      for(i=6;i<Nvar;i++)
         yscal[i]=(fabs(y[i])+hdid*fabs(yscal[i]));

      for(i=0;i<Nvar;i++){
         if(y[i]==0)
            yscal[i]=1000000000.0;
      }



   } //fin del while***************************


   delete y;
   delete yout;
   delete yscal;
   delete dydx;
   delete vv;
   delete Q;
   delete R;
   delete colpoints;
   delete aux;
   delete Pl;
   delete Ch;
   delete Tiempo;
   delete Ann;
   delete Bnn;
   delete Radio;
   delete Masa;
   delete k1;
   delete k2;

   if(iespectro==6)
      out1<<setprecision(12)<<fluencia<<endl;

   if(iespectro==6){
      for(i=0;i<=int(lambdapuntos);i++)
         out3<<setprecision(12)<<((lambdainicial+i*(lambdafinal-lambdainicial)/lambdapuntos)*1.0e-6)<<" "<<
         electra[i]<<endl;
   }

   if(iespectro==6)
      out5<<fotones<<endl;


   cout<<"termine"<<endl;
   getch();
   out1.close();
   out2.close();
   out3.close();
   out4.close();
   out5.close();
   out6.close();
   out7.close();



} //fin del main


