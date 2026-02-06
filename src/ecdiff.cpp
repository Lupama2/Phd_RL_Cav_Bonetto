#include<math.h>
#include<iostream.h>

  
double erf(double Omega, int N);
double PolCheb(double y, int NCh);
double DerPCheb(double y, int NCh);
double DerSCheb(double y, int NCh);
void MasaMolar(double *Masa);
void Conductividad1(double *k1);
void Conductividad2(double *k2);
void reacciones(double *y,double VolExcl,double ntot,double V,int Nvar2,double *r);
void DiamColVect(double *DiamCol);
void fraccion(double *fracc);

void derivada(double t,double *y,double TB, double mp,double k,
              double *colpoints,double **Q,double **R, double *dydt){



   //t: tiempo
   //vector y: y[1]-> radio, y[2]-> velocidad, y[3]-> temperatura , y[4]-> nh2o, y[5]-> TLi
   //nt: cantidad de particulas totales (y[4]+nAr)
   //T: temperatura en el interior de la burbuja
   //TB: temperatura del gas en la superficie de la burbuja
   //mp: es la evaporacion neta en la pared de la burbuja por unidad de area y undad de tiempo (kg/(m^2*seg))
   //k: conductivad termica
   //landa: el camino libre medio
   //dydt: el vector que tiene las derivadas respecto del tiempo del vector y; dydt[1]=velocidad; dydt[2]=aceleracion

   extern double Rg,Na,sigma,mu,rol,n,Pinf,pvap0,cl,sigmap,Cvh2o,CvAr,CvO2,CvN2,Ri,RAr,Dif,f,gama;
   extern double Mh2o,MAr,MO2,MN2,MH2,MH,MO,MOH,MN,MNH,MNH2,MNH3,MN2H,MNO,MNO2,MN2O,MHNO;
   extern double Amp,w,fi,Amp2,w2,fi2,frec2;
   extern double e0,L,kl,cpl,Tinf,T0;
   extern double a_h2o,a_Ar,b_h2o,b_Ar,a_Arh2o,b_Arh2o;
   extern double nt0Ar,nh2o0;
   extern double ngas0,nh2o0,nh20,nh0,no0,no20,noh0,nn0,nnh0,nnh20,nnh30,nn2h0;
   extern double nno0,nno20,nn2o0,nhno0,nn20,nar0,nh2o20,nho20;
   extern double alfaM,ap,alfae,K,Rv;
   extern int N,Npoints;
   extern int rapidez,modelo;
   extern double k1H2,k1H,k1O,k1O2,k1OH,k1h2o,k1N,k1NH,k1NH2,k1NH3,k1N2H,k1NO,k1NO2,k1N2O,k1HNO,k1N2,k1Ar;
   extern double k2H2,k2H,k2O,k2O2,k2OH,k2h2o,k2N,k2NH,k2NH2,k2NH3,k2N2H,k2NO,k2NO2,k2N2O,k2HNO,k2N2,k2Ar;
   extern int NCh,NCh2,NCh3;//numero de polinomios que uso para calcular la T(r,t)
   extern int Nint,NumEsp; //numero de intervalos para integrar las funciones
   extern double nu;
   extern int nn;
   extern int yndis;
   extern double qpmax,tp;
   extern double ngas0;
   extern double Intensidad,duracion;
   extern int yndisoc,yndif;
   extern double Ppert;

   int i,j;
   int cont,icont,iw;
   double sum,aux;

   double pg; //presion en el gas (Vander Waals)
   double ps; //presion acustica
   double pv; //presion de vapor de agua
   double V; //volumen total
   double v; //volumen molar
   double pB; //Presion en liquido en la superficie externa de la bubrbuja
   double rog; //densidad del gas
   double dpB1dt,dpB2dt;//en realidad es el que en el cuaderno llame prima
   double dpgdt,dpg2dt; //derivada de la presion en el gas respecto del tiempo
   double drogdt; //derivada de la densidad del gas respecto del tiempo
   double ngas,dntdt,ntot; //derivada de la cantidad de particulas totales respecto del tiempo
   double dvdt; //derivada del volumen molar respecto del tiempo
   double dTdt; //derivada de la temperatura respecto del tiempo
   double dadt,dbdt; //derivadas de los coeficientes de VW
   double eh2o; //una variable que depende de la temperatura en la superficie de la burbuja
   double M,m; //masa total de gas en la burbuja
   double landa;
   double a,b;
   double lth,Xi;
   double pvap,dpvapdt,dTB1dt,dTB2dt,dT2dt,Gama,Omega,dpv1dt,dpv2dt,dk1dt,dk2dt,dmdt,dnpdt,dlandadt;
   double np;
   double mpp1,mpp2;
   double qpp,dq1dt,dq2dt,dqdTB;
   double dlth1dt,dlth2dt;
   double DPTCh,DSTCh;
   double Pulso;

   int NumEspm1=NumEsp+1;

   double *F=new double[NCh2];
   double *x=new double[NCh2];
   double *Masa=new double[NumEspm1];
   double *k1=new double[NumEspm1];
   double *k2=new double[NumEspm1];

   double *r=new double[102];
   double *y0=new double[NumEspm1];
   double *DiamCol=new double[NumEspm1];
   double *fracc=new double[20];
   long double *InvDi=new long double[20];
   double *ldiffi=new double[20];

   double f1,f2,f3,f4,f5,Pp;

   int Nvar2=(NCh3-1);

   //definicion de las 45 reacciones pertenecientes a los 16 reactivos + gas noble
   //rf reaccion forward, rb reaccion backward
   double rf1,rb1,rf2,rb2,rf3,rb3,rf4,rb4,rf5,rb5,rf6,rb6,rf7,rb7,rf8,rb8,rf9,rb9;
   double rf10,rb10,rf11,rb11,rf12,rb12,rf13,rb13,rf14,rb14,rf15,rb15,rf16,rb16;
   double rf17,rb17,rf18,rb18,rf19,rb19,rf20,rb20,rf21,rb21,rf22,rb22,rf23,rb23;
   double rf24,rb24,rf25,rb25,rf26,rb26,rf27,rb27,rf28,rb28,rf29,rb29,rf30,rb30;
   double rf31,rb31,rf32,rb32,rf33,rb33,rf34,rb34,rf35,rb35,rf36,rb36,rf37,rb37;
   double rf38,rb38,rf39,rb39,rf40,rb40,rf41,rb41,rf42,rb42,rf43,rb43,rf44,rb44;
   double rf45,rb45,rf46,rb46,rf47,rb47,rf48,rb48,rf49,rb49,rf50,rb50,rf51,rb51;

   double VolExcl;
   double MolCov=5.1e-29; //molecular covolumen [m^3]
   double CalorDifusion,CalorDisociacion,CalorTermico;
   double et1,et2,et3,et4,et5,et6,et7,et8,et9,et10,et11,et12,et13,et14,et15,et16,et17,et18,et19;
   double det1,det2,det3,det4,det5,det6,det7,det8,det9,det10,det11,det12,det13,det14,det15,det16,det17,det18,det19;


   double pi=3.14159;


   //***************************************************************************
   //***************************************************************************
   //***************************************************************************

   V=4.0/3.0*pi*pow(y[1],3.0); //Volumen

   MasaMolar(Masa);

   icont=1;

   M=0.0; //Masa total

   while(icont<=NumEsp){
      M=M+(Masa[icont]*y[Nvar2+icont])/Na;
      icont=icont+1;
   }


   ntot=0.0; //Numero total de particulas (condensables + no-condensables)

   icont=1;

   while(icont<=NumEsp){
      ntot=ntot+y[Nvar2+icont];
      icont=icont+1;
   }

   m=M/ntot; //masa media de una particula

   landa=V/(sqrt(2.0)*sigmap*ntot);  //camino libre medio

   v=V/ntot*Na; //volumen molar

   np=ntot/V; //densidad de particulas

   eh2o=Cvh2o/Na*TB; //energia llevada por una molecula que se esta evaporando o codensando


   ngas=ntot-y[Nvar2+6]; //numero de particulas de no-condensables

   a = a_Ar*(ngas/ntot)*(ngas/ntot)+2.0*a_Arh2o*(ngas/ntot)*(y[Nvar2+6]/ntot)+
       a_h2o*(y[Nvar2+6]/ntot)*(y[Nvar2+6]/ntot);
   //const a de VW promediada con la cantidad de particulas

   b = b_Ar*(ngas/ntot)*(ngas/ntot)+2.0*b_Arh2o*(ngas/ntot)*(y[Nvar2+6]/ntot)+
       b_h2o*(y[Nvar2+6]/ntot)*(y[Nvar2+6]/ntot);
   //const b de VW promediada con la cantidad de particulas

   pg=Rg*y[3]/(v-b)-a/(v*v);//presion total

   rog=M/V; //densidad de la mezcla

   pB=pg-2.0*sigma/y[1]-4.0*mu/y[1]*(y[2]-mp/rol)-mp*mp*(1.0/rol-1.0/rog);
   //presion en la interfase del lado del liquido

   pv=y[Nvar2+6]/ntot*pg; //presion parcial del vapor de agua

   //ldiff longitud de difusion

   //***************************************************************************
  //************bloque que calcula el cambio debido a la difusion***************
  //****************************************************************************

   y0[1]=nh20;
   y0[2]=nh0;
   y0[3]=no0;
   y0[4]=no20;
   y0[5]=noh0;
   y0[6]=nh2o0;
   y0[7]=nn0;
   y0[8]=nnh0;
   y0[9]=nnh20;
   y0[10]=nnh30;
   y0[11]=nn2h0;
   y0[12]=nno0;
   y0[13]=nno20;
   y0[14]=nn2o0;
   y0[15]=nhno0;
   y0[16]=nn20;
   y0[17]=nar0;
   y0[18]=nh2o20;
   y0[19]=nho20;

   DiamColVect(DiamCol);


   for(iw=1;iw<20;iw++){
      InvDi[iw]=0.0;

      for(j=1;j<20;j++){

         if((y[Nvar2+iw] != 0) && (y[Nvar2+j] != 0))
            InvDi[iw]+=(y[Nvar2+j]/ntot)/(1.0-y[Nvar2+iw]/ntot)/
            (3.0/8.0*sqrt(fabs(pi*K*Tinf/(2.0*Masa[iw]*Masa[j]/
            (Masa[iw]+Masa[j])))/
            ((y[Nvar2+iw]+y[Nvar2+j])/V*pi*0.5*(DiamCol[iw]+DiamCol[j])*
            0.5*(DiamCol[iw]+DiamCol[j]))));

      }

      if((y[Nvar2+iw] != 0))
         InvDi[iw]-=(y[Nvar2+iw]/ntot)/(1.0-y[Nvar2+iw]/ntot)/
         (3.0/8.0*sqrt(fabs(pi*K*Tinf/(2.0*Masa[iw]*Masa[iw]/
         (Masa[iw]+Masa[iw])))/
         ((y[Nvar2+iw]+y[Nvar2+iw])/V*pi*0.5*(DiamCol[iw]+DiamCol[iw])*
         0.5*(DiamCol[iw]+DiamCol[iw]))));

      if(InvDi[iw]<1.0e-9)
         InvDi[iw]=0.0;

   }


   for(i=1;i<20;i++){

      if(InvDi[i] != 0){
         if(y[2]==0.0)
            ldiffi[i]=y[1]/pi;
         else{
            if(sqrt(y[1]/InvDi[i]/fabs(y[2])) < (y[1]/pi))
               ldiffi[i]=sqrt(y[1]/InvDi[i]/fabs(y[2]));
            else
               ldiffi[i]=y[1]/pi;
         }
      }

      else
         ldiffi[i]=y[1]/pi;


      if(yndif != 0){
         if(InvDi[i] != 0)
            dydt[Nvar2+i]=4.0*pi*y[1]*y[1]/InvDi[i]*(y0[i]/(4.0*pi/3.0*Ri*Ri*Ri)-
            y[Nvar2+i]/(4.0*pi/3.0*y[1]*y[1]*y[1]))/ldiffi[i];
         else
            dydt[Nvar2+i]=0.0;
      }
      else
         dydt[Nvar2+i]=0.0;      

   }

   //***************************************************************************
   //***************************************************************************
   //***************************************************************************

   //Antes de calcular la variacion de particulas por disociacion, necesito solo
   //por difusion

   fraccion(fracc);

   CalorDifusion=0.0;

   for(i=1;i<20;i++)
     CalorDifusion=CalorDifusion+(1.0+fracc[i]/2.0)*K*Tinf*dydt[Nvar2+i];


   //***************************************************************************
   //***********bloque que calcula disociacion de las distintas especies********
   //***************************************************************************

   //y[Nvar2+1] =H2
   //y[Nvar2+2] =H
   //y[Nvar2+3] =O
   //y[Nvar2+4] =O2
   //y[Nvar2+5] =OH
   //y[Nvar2+6] =H2O = y[4] = vapor
   //y[Nvar2+7] =N
   //y[Nvar2+8] =NH
   //y[Nvar2+9] =NH2
   //y[Nvar2+10] =NH3
   //y[Nvar2+11] =N2H
   //y[Nvar2+12] =NO
   //y[Nvar2+13] =NO2
   //y[Nvar2+14] =N2O
   //y[Nvar2+15] =HNO
   //y[Nvar2+16] =N2
   //y[Nvar2+17] =Ar (gas noble)
   //y[Nvar2+18] =H2O2
   //y[Nvar2+19] =HO2


   VolExcl=exp(ntot/V*MolCov/(1.0-ntot/V*MolCov))/(1.0-ntot/V*MolCov);

   //rfi: reacciones forward para la especie Nvar2+i
   //rbi: reacciones backward para la especie Nvar2+i
   reacciones(y,VolExcl,ntot,V,Nvar2,r);

   if(yndisoc !=0){
      rf1=r[0];
      rb1=r[1];
      rf2=r[2];
      rb2=r[3];
      rf3=r[4];
      rb3=r[5];
      rf4=r[6];
      rb4=r[7];
      rf5=r[8];
      rb5=r[9];
      rf6=r[10];
      rb6=r[11];
      rf7=r[12];
      rb7=r[13];
      rf8=r[14];
      rb8=r[15];
      rf9=r[16];
      rb9=r[17];
      rf10=r[18];
      rb10=r[19];
      rf11=r[20];
      rb11=r[21];
      rf12=r[22];
      rb12=r[23];
      rf13=r[24];
      rb13=r[25];
      rf14=r[26];
      rb14=r[27];
      rf15=r[28];
      rb15=r[29];
      rf16=r[30];
      rb16=r[31];
      rf17=r[32];
      rb17=r[33];
      rf18=r[34];
      rb18=r[35];
      rf19=r[36];
      rb19=r[37];
      rf20=r[38];
      rb20=r[39];
      rf21=r[40];
      rb21=r[41];
      rf22=r[42];
      rb22=r[43];
      rf23=r[44];
      rb23=r[45];
      rf24=r[46];
      rb24=r[47];
      rf25=r[48];
      rb25=r[49];
      rf26=r[50];
      rb26=r[51];
      rf27=r[52];
      rb27=r[53];
      rf28=r[54];
      rb28=r[55];
      rf29=r[56];
      rb29=r[57];
      rf30=r[58];
      rb30=r[59];
      rf31=r[60];
      rb31=r[61];
      rf32=r[62];
      rb32=r[63];
      rf33=r[64];
      rb33=r[65];
      rf34=r[66];
      rb34=r[67];
      rf35=r[68];
      rb35=r[69];
      rf36=r[70];
      rb36=r[71];
      rf37=r[72];
      rb37=r[73];
      rf38=r[74];
      rb38=r[75];
      rf39=r[76];
      rb39=r[77];
      rf40=r[78];
      rb40=r[79];
      rf41=r[80];
      rb41=r[81];
      rf42=r[82];
      rb42=r[83];
      rf43=r[84];
      rb43=r[85];
      rf44=r[86];
      rb44=r[87];
      rf45=r[88];
      rb45=r[89];
      rf46=r[90];
      rb46=r[91];
      rf47=r[92];
      rb47=r[93];
      rf48=r[94];
      rb48=r[95];
      rf49=r[96];
      rb49=r[97];
      rf50=r[98];
      rb50=r[99];
      rf51=r[100];
      rb51=r[101];

      //************************************************************************
      //*******************derivadas de las especies****************************
      //************************difusion + disociacion**************************


      dydt[Nvar2+1]=dydt[Nvar2+1]+V*(-1.0*(rf3-rb3)+1.0*(rf5-rb5)-1.0*(rf7-rb7)+1.0*(rf20-rb20)+
      1.0*(rf26-rb26)+1.0*(rf31-rb31)+1.0*(rf37-rb37)+1.0*(rf41-rb41)+1.0*
      (rf43-rb43)+1.0*(rf49-rb49));

      dydt[Nvar2+2]=dydt[Nvar2+2]+V*(-1.0*(rf2-rb2)+1.0*(rf3-rb3)-1.0*(rf4-rb4)-2.0*(rf5-rb5)-
      1.0*(rf6-rb6)+1.0*(rf7-rb7)+1.0*(rf11-rb11)-1.0*(rf14-rb14)-1.0*(rf18-rb18)+
      1.0*(rf19-rb19)-1.0*(rf20-rb20)+1.0*(rf21-rb21)+1.0*(rf25-rb25)+1.0*(rf28-rb28)+
      +1.0*(rf30-rb30)-1.0*(rf31-rb31)+1.0*(rf33-rb33)+1.0*(rf34-rb34)-1.0*(rf37-rb37)-
      1.0*(rf39-rb39)-1.0*(rf41-rb41)-1.0*(rf43-rb43)-1.0*(rf48-rb48)-1.0*(rf49-rb49));

      dydt[Nvar2+3]=dydt[Nvar2+3]+V*(-2.0*(rf1-rb1)-1.0*(rf2-rb2)-1.0*(rf3-rb3)+1.0*(rf4-rb4)
      +1.0*(rf8-rb8)+1.0*(rf9-rb9)+1.0*(rf10-rb10)-1.0*(rf12-rb12)-1.0*(rf13-rb13)+
      1.0*(rf15-rb15)-1.0*(rf16-rb16)-1.0*(rf17-rb17)-1.0*(rf19-rb19)+1.0*
      (rf23-rb23)-1.0*(rf29-rb29)-1.0*(rf30-rb30)-1.0*(rf35-rb35)-1.0*(rf36-rb36)
      -1.0*(rf40-rb40)-1.0*(rf45-rb45)-1.0*(rf50-rb50));

      dydt[Nvar2+4]=dydt[Nvar2+4]+V*(1.0*(rf1-rb1)-1.0*(rf4-rb4)-1.0*(rf10-rb10)+1.0*(rf12-rb12)+
      1.0*(rf17-rb17)-1.0*(rf23-rb23)-1.0*(rf24-rb24)+1.0*(rf46-rb46));

      dydt[Nvar2+5]=dydt[Nvar2+5]+V*(+1.0*(rf2-rb2)+1.0*(rf3-rb3)+1.0*(rf4-rb4)-1.0*(rf6-rb6)
      -1.0*(rf7-rb7)-2.0*(rf8-rb8)-1.0*(rf11-rb11)+1.0*(rf14-rb14)+1.0*(rf18-rb18)
      -1.0*(rf21-rb21)-1.0*(rf22-rb22)+1.0*(rf24-rb24)+1.0*(rf27-rb27)+1.0*(rf29-rb29)
      -1.0*(rf32-rb32)+1.0*(rf35-rb35)-1.0*(rf38-rb38)+1.0*(rf40-rb40)-1.0*(rf42-rb42)
      -1.0*(rf44-rb44)+1.0*(rf45-rb45)+2.0*(rf47-rb47)+1.0*(rf48-rb48)+1.0*(rf50-rb50));

      dydt[Nvar2+6]=dydt[Nvar2+6]+V*(+1.0*(rf6-rb6)+1.0*(rf7-rb7)+1.0*(rf8-rb8)+1.0*(rf22-rb22)
      -1.0*(rf26-rb26)+1.0*(rf32-rb32)+1.0*(rf38-rb38)+1.0*(rf42-rb42)+1.0*(rf44-rb44)+
      1.0*(rf48-rb48)+1.0*(rf51-rb51));

      /*
      dydt[Nvar2+7]=dydt[Nvar2+7]+V*(-1.0*(rf9-rb9)-1.0*(rf10-rb10)-1.0*(rf11-rb11)+1.0*(rf20-rb20)+1.0*(rf22-rb22)
      -1.0*(rf25-rb25));

      dydt[Nvar2+8]=dydt[Nvar2+8]+V*(-1.0*(rf19-rb19)-1.0*(rf20-rb20)-1.0*(rf21-rb21)-1.0*(rf22-rb22)
      -1.0*(rf23-rb23)-1.0*(rf24-rb24)-1.0*(rf25-rb25)-1.0*(rf26-rb26)-1.0*
      (rf27-rb27)-1.0*(rf28-rb28)+1.0*(rf29-rb29)+1.0*(rf31-rb31)+1.0*(rf32-rb32)
      +1.0*(rf36-rb36));

      dydt[Nvar2+9]=dydt[Nvar2+9]+V*(-1.0*(rf29-rb29)-1.0*(rf30-rb30)-1.0*(rf31-rb31)-1.0*(rf32-rb32)
      +1.0*(rf43-rb43)+1.0*(rf44-rb44)+1.0*(rf45-rb45));

      dydt[Nvar2+10]=dydt[Nvar2+10]+V*(-1.0*(rf43-rb43)-1.0*(rf44-rb44)-1.0*(rf45-rb45));

      dydt[Nvar2+11]=dydt[Nvar2+11]+V*(-1.0*(rf33-rb33)-1.0*(rf34-rb34)-1.0*(rf35-rb35)-1.0*(rf36-rb36)
      -1.0*(rf37-rb37)-1.0*(rf38-rb38));

      dydt[Nvar2+12]=dydt[Nvar2+12]+V*(-1.0*(rf9-rb9)+1.0*(rf10-rb10)+1.0*(rf11-rb11)+2.0*(rf13-rb13)
      -1.0*(rf16-rb16)+1.0*(rf17-rb17)+1.0*(rf18-rb18)+1.0*(rf19-rb19)+1.0*(rf24-rb24)
      -1.0*(rf27-rb27)-1.0*(rf28-rb28)+1.0*(rf36-rb36)-1.0*(rf39-rb39)+1.0*(rf40-rb40)
      +1.0*(rf41-rb41)+1.0*(rf42-rb42));

      dydt[Nvar2+13]=dydt[Nvar2+13]+V*(+1.0*(rf16-rb16)-1.0*(rf17-rb17)-1.0*(rf18-rb18));

      dydt[Nvar2+14]=dydt[Nvar2+14]+V*(-1.0*(rf12-rb12)-1.0*(rf13-rb13)-1.0*(rf14-rb14)-1.0*(rf15-rb15)+
      1.0*(rf28-rb28));

      dydt[Nvar2+15]=dydt[Nvar2+15]+V*(+1.0*(rf21-rb21)+1.0*(rf23-rb23)+1.0*(rf26-rb26)+1.0*(rf30-rb30)
      +1.0*(rf39-rb39)-1.0*(rf40-rb40)-1.0*(rf41-rb41)-1.0*(rf42-rb42));

      dydt[Nvar2+16]=dydt[Nvar2+16]+V*(+1.0*(rf9-rb9)+1.0*(rf12-rb12)+1.0*(rf14-rb14)+1.0*(rf15-rb15)+
      1.0*(rf25-rb25)+1.0*(rf27-rb27)+1.0*(rf33-rb33)+1.0*(rf34-rb34)+1.0*
      (rf35-rb35)+1.0*(rf37-rb37)+1.0*(rf38-rb38));
      */


      dydt[Nvar2+17]=dydt[Nvar2+17];

      dydt[Nvar2+18]=dydt[Nvar2+18]+V*(-1.0*(rf47-rb47)-1.0*(rf48-rb48)-1.0*(rf49-rb49)-1.0*(rf50-rb50)-
      1.0*(rf51-rb51)+1.0*(rf46-rb46));

      dydt[Nvar2+19]=dydt[Nvar2+19]+V*(-2.0*(rf46-rb46)+1.0*(rf49-rb49)+1.0*(rf50-rb50)+1.0*(rf51-rb51));

   }

   dydt[Nvar2+6]=dydt[Nvar2+6]+4.0*pi*y[1]*y[1]*(mp*Na/Mh2o); //variacion de la cantidad de particulas
   //de vapor de agua debido a la condensacion, mas por difusion y reacciones quimicas


   dntdt=0.0;

   icont=1;

   while(icont<=NumEsp){
      dntdt=dntdt+dydt[Nvar2+icont];
      icont=icont+1;
   }

   dadt=-2.0/pow(ntot,3.0)*(a_Ar*ngas*ngas+2.0*a_Arh2o*ngas*y[Nvar2+6]-a_Arh2o*ngas*ntot+
   a_h2o*y[Nvar2+6]*y[Nvar2+6]-a_h2o*y[Nvar2+6]*ntot)*dntdt;
   //derivada respecto del tiempo de la const. a de VW

   dbdt=-2.0/pow(ntot,3.0)*(b_Ar*ngas*ngas+2.0*b_Arh2o*ngas*y[Nvar2+6]-b_Arh2o*ngas*ntot+
   b_h2o*y[Nvar2+6]*y[Nvar2+6]-b_h2o*y[Nvar2+6]*ntot)*dntdt;
   //derivada respecto del tiempo de la const. b de VW

   dvdt=-V*Na/(ntot*ntot)*dntdt+4.0*pi*y[1]*y[1]*Na/ntot*y[2];
   //derivada del volumen molar con respecto a t

   Xi=k/(8.0/2.0*nh2o0*K+5.0/2.0*ngas*K);

   if(y[2] == 0)
     lth=y[1]/pi;
   else{
      if((y[1]/pi)<(sqrt(y[1]*Xi/fabs(y[2]))))
	 lth=y[1]/pi;
      else
	 lth=sqrt(y[1]*Xi/fabs(y[2]));
   }


   if(modelo == 1) //Yasui
      qpp=k*(TB-y[3])/n/landa;

   if(modelo == 2) //Toegel
      qpp=k*(Tinf-y[3])/lth;

   if(modelo == 3){ //Hao-Prosperetti
      cont=NCh;
      qpp=0.;

      while(cont >= 0){
         DPTCh=DerPCheb(1.0,2.0*cont);
         qpp=qpp+y[cont+6]*DPTCh;
         cont=cont-1;
      }

      qpp=qpp/y[1];

   }

   if(modelo == 4){
      cont=NCh;
      qpp=0.;

      while(cont >= 2){
         DPTCh=DerPCheb(1.0,2.0*cont);
         qpp=qpp+y[cont+6]*DPTCh;
         cont=cont-1;
      }

      qpp=qpp/y[1]+k*(Tinf-y[3])/lth;;


   }

   if(yndisoc != 0){
      CalorDisociacion=1.0e3*V/Na*((rf1-rb1)*498.0+
      (rf2-rb2)*428.0+(rf3-rb3)*(-8.0)+(rf4-rb4)*(-70.0)+
      (rf5-rb5)*436.0+(rf6-rb6)*498.0+(rf7-rb7)*62.0+(rf8-rb8)*70.0+
      (rf9-rb9)*314.0+(rf10-rb10)*133.0+(rf11-rb11)*203.0+(rf12-rb12)*331.0+
      (rf13-rb13)*150.0+(rf14-rb14)*261.0+(rf15-rb15)*(-167.0)+
      (rf16-rb16)*305.0+(rf17-rb17)*193.0+(rf18-rb18)*123.0+
      (rf19-rb19)*292.0+(rf20-rb20)*97.0+(rf21-rb21)*68.0+
      (rf22-rb22)*159.0+(rf23-rb23)*(-15.0)+(rf24-rb24)*222.0+
      (rf25-rb25)*606.0+(rf26-rb26)*(-70.0)+(rf27-rb27)*403.0+
      (rf28-rb28)*143.0+(rf29-rb29)*43.0+(rf30-rb30)*98.0+
      (rf31-rb31)*51.0+(rf32-rb32)*113.0+(rf33-rb33)*32.0+
      (rf34-rb34)*32.0+(rf35-rb35)*(-485.0)+(rf36-rb36)*56.0+
      (rf37-rb37)*467.0+(rf38-rb38)*531.0+(rf39-rb39)*207.0+
      (rf40-rb40)*221.0+(rf41-rb41)*229.0+(rf42-rb42)*291.0+
      (rf43-rb43)*(-13.0)+(rf44-rb44)*49.0+(rf45-rb45)*(-21.0)+
      (rf46-rb46)*(175.35)+(rf47-rb47)*(-217.89)+(rf48-rb48)*(290.93)+
      (rf49-rb49)*(64.32)+(rf50-rb50)*(56.08)+(rf51-rb51)*(128.67));
   }
   else
      CalorDisociacion=0.0;

   et1=5.0/2.0*K*y[3]+K*6325.0/(exp(6325.0/y[3])-1.0);
   et2=3.0/2.0*K*y[3];
   et3=3.0/2.0*K*y[3];
   et4=5.0/2.0*K*y[3]+K*2273.0/(exp(2273.0/y[3])-1.0);
   et5=5.0/2.0*K*y[3]+K*5370.0/(exp(5370.0/y[3])-1.0);
   et6=6.0/2.0*K*y[3]+K*2295.0/(exp(2295.0/y[3])-1.0)+K*5255.0/(exp(5255.0/y[3])
       -1.0)+K*5400.0/(exp(5400.0/y[3])-1.0);
   et7=3.0/2.0*K*y[3];
   et8=5.0/2.0*K*y[3]+K*3445.0/(exp(3445.0/y[3])-1.0);
   et9=6.0/2.0*K*y[3]+K*4628.0/(exp(4628.0/y[3])-1.0)+K*2152.0/(exp(2152.0/y[3])
       -1.0)+K*4746.0/(exp(4746.0/y[3])-1.0);
   et10=6.0/2.0*K*y[3]+K*4798.0/(exp(4798.0/y[3])-1.0)+K*1366.0/(exp(1366.0/y[3])
        -1.0)+K*4951.0/(exp(4951.0/y[3])-1.0)
	  +K*4951.0/(exp(4951.0/y[3])-1.0)+K*2339.0/(exp(2339.0/y[3])-1.0)+
          K*2339.0/(exp(2339.0/y[3])-1.0);
   et11=5.0/2.0*K*y[3]+K*4658.0/(exp(4658.0/y[3])-1.0)+K*990.0/(exp(990.0/y[3])-1.0)
        +K*3253.0/(exp(3253.0/y[3])-1.0);
   et12=5.0/2.0*K*y[3]+K*2710.0/(exp(2710.0/y[3])-1.0);
   et13=6.0/2.0*K*y[3]+K*1878.0/(exp(1878.0/y[3])-1.0)+K*1064.0/(exp(1064.0/y[3])
        -1.0)+K*2272.0/(exp(2272.0/y[3])-1.0);
   et14=5.0/2.0*K*y[3]+K*3098.0/(exp(3098.0/y[3])-1.0)+K*822.0/(exp(822.0/y[3])-1.0)
        +K*1819.0/(exp(1819.0/y[3])-1.0);
   et15=5.0/2.0*K*y[3]+K*3895.0/(exp(3895.0/y[3])-1.0)+K*2158.0/(exp(2158.0/y[3])
        -1.0)+K*2250.0/(exp(2250.0/y[3])-1.0);
   et16=5.0/2.0*K*y[3]+K*3350.0/(exp(3350.0/y[3])-1.0);
   et17=3.0/2.0*K*y[3];
   et18=6.0/2.0*K*y[3];
   et19=6.0/2.0*K*y[3];

   CalorTermico=-(et1*dydt[Nvar2+1]+et2*dydt[Nvar2+2]+et3*dydt[Nvar2+3]+
	       	et4*dydt[Nvar2+4]+et5*dydt[Nvar2+5]+et6*dydt[Nvar2+6]+et7*dydt[Nvar2+7]+
	       	et8*dydt[Nvar2+8]+et9*dydt[Nvar2+9]+et10*dydt[Nvar2+10]+et11*dydt[Nvar2+11]+
	       	et12*dydt[Nvar2+12]+et13*dydt[Nvar2+13]+et14*dydt[Nvar2+14]+et15*dydt[Nvar2+15]
	       	+et16*dydt[Nvar2+16]+et17*dydt[Nvar2+17]+et18*dydt[Nvar2+18]+et19*dydt[Nvar2+19]);


   det1=5.0/2.0*K+K*pow(6325.0,2.0)/(exp(6325.0/y[3])-1.0)/(exp(6325.0/y[3])-1.0)/y[3]/y[3];
   det2=3.0/2.0*K;
   det3=3.0/2.0*K;
   det4=5.0/2.0*K+K*2273.0*2273.0/(exp(2273.0/y[3])-1.0)/(exp(2273.0/y[3])-1.0)/y[3]/y[3];
   det5=5.0/2.0*K+K*5370.0*5370.0/(exp(5370.0/y[3])-1.0)/(exp(5370.0/y[3])-1.0)/y[3]/y[3];
   det6=6.0/2.0*K+K*2295.0*2295.0/(exp(2295.0/y[3])-1.0)/(exp(2295.0/y[3])-1.0)/y[3]/y[3]
        +K*5255.0*5255.0/(exp(5255.0/y[3])-1.0)/(exp(5255.0/y[3])-1.0)/y[3]/y[3]+
        K*5400.0*5400.0/(exp(5400.0/y[3])-1.0)/(exp(5400.0/y[3])-1.0)/y[3]/y[3];
   det7=3.0/2.0*K;
   det8=5.0/2.0*K+K*3445.0*3445.0/(exp(3445.0/y[3])-1.0)/(exp(3445.0/y[3])-1.0)/y[3]/y[3];
   det9=6.0/2.0*K+K*4628.0*4628.0/(exp(4628.0/y[3])-1.0)/(exp(4628.0/y[3])-1.0)/y[3]/y[3]
        +K*2152.0*2152.0/(exp(2152.0/y[3])-1.0)/(exp(2152.0/y[3])-1.0)/y[3]/y[3]+
        K*4746.0*4746.0/(exp(4746.0/y[3])-1.0)/(exp(4746.0/y[3])-1.0)/y[3]/y[3];
   det10=6.0/2.0*K+K*4798.0*4798.0/(exp(4798.0/y[3])-1.0)/(exp(4798.0/y[3])-1.0)/y[3]/y[3]
         +K*1366.0*1366.0/(exp(1366.0/y[3])-1.0)/(exp(1366.0/y[3])-1.0)/y[3]/y[3]
         +K*4951.0*4951.0/(exp(4951.0/y[3])-1.0)/(exp(4951.0/y[3])-1.0)/y[3]/y[3]
         +K*4951.0*4951.0/(exp(4951.0/y[3])-1.0)/(exp(4951.0/y[3])-1.0)/y[3]/y[3]+
         K*2339.0*2339.0/(exp(2339.0/y[3])-1.0)/(exp(2339.0/y[3])-1.0)/y[3]/y[3]+
         K*2339.0*2339.0/(exp(2339.0/y[3])-1.0)/(exp(2339.0/y[3])-1.0)/y[3]/y[3];
   det11=5.0/2.0*K+K*4658.0*4658.0/(exp(4658.0/y[3])-1.0)/(exp(4658.0/y[3])-1.0)/y[3]/y[3]
         +K*990.0*990.0/(exp(990.0/y[3])-1.0)/(exp(990.0/y[3])-1.0)/y[3]/y[3]+
         K*3253.0*3253.0/(exp(3253.0/y[3])-1.0)/(exp(3253.0/y[3])-1.0)/y[3]/y[3];
   det12=5.0/2.0*K+K*2710.0*2710.0/(exp(2710.0/y[3])-1.0)/(exp(2710.0/y[3])-1.0)/y[3]/y[3];
   det13=6.0/2.0*K+K*1878.0*1878.0/(exp(1878.0/y[3])-1.0)/(exp(1878.0/y[3])-1.0)/y[3]/y[3]
         +K*1064.0*1064.0/(exp(1064.0/y[3])-1.0)/(exp(1064.0/y[3])-1.0)/y[3]/y[3]+
         K*2272.0*2272.0/(exp(2272.0/y[3])-1.0)/(exp(2272.0/y[3])-1.0)/y[3]/y[3];
   det14=5.0/2.0*K+K*3098.0*3098.0/(exp(3098.0/y[3])-1.0)/(exp(3098.0/y[3])-1.0)/y[3]/y[3]+
         K*822.0*822.0/(exp(822.0/y[3])-1.0)/(exp(822.0/y[3])-1.0)/y[3]/y[3]
         +K*1819.0*1819.0/(exp(1819.0/y[3])-1.0)/(exp(1819.0/y[3])-1.0)/y[3]/y[3];
   det15=5.0/2.0*K+K*3895.0*3895.0/(exp(3895.0/y[3])-1.0)/(exp(3895.0/y[3])-1.0)/y[3]/y[3]
         +K*2158.0*2158.0/(exp(2158.0/y[3])-1.0)/(exp(2158.0/y[3])-1.0)/y[3]/y[3]
         +K*2250.0*2250.0/(exp(2250.0/y[3])-1.0)/(exp(2250.0/y[3])-1.0)/y[3]/y[3];
   det16=5.0/2.0*K+K*3350.0*3350.0/(exp(3350.0/y[3])-1.0)/(exp(3350.0/y[3])-1.0)/y[3]/y[3];
   det17=3.0/2.0*K;
   det18=6.0/2.0*K;
   det19=6.0/2.0*K;

   if(t<=(duracion/2.0))
      Pulso=Intensidad*t/(duracion/2.0);

   else{
      if((t>=(duracion/2.0)) && (t<=duracion))
         Pulso=(-2.0*Intensidad/duracion)*t+2.0*Intensidad;
      else
         Pulso=0.0;
   }


   dTdt=(0.*Pulso-4.0*pi*y[1]*y[1]*y[2]*pg+4.0*pi*y[1]*y[1]*mp*eh2o*Na/Mh2o+4.0*pi*y[1]*y[1]*qpp
         +CalorDifusion+CalorDisociacion+CalorTermico)/(det1*y[Nvar2+1]+det2*y[Nvar2+2]+det3*y[Nvar2+3]+
         det4*y[Nvar2+4]+det5*y[Nvar2+5]+det6*y[Nvar2+6]+det7*y[Nvar2+7]+det8*y[Nvar2+8]+
         det9*y[Nvar2+9]+det10*y[Nvar2+10]+det11*y[Nvar2+11]+det12*y[Nvar2+12]+
         det13*y[Nvar2+13]+det14*y[Nvar2+14]+det15*y[Nvar2+15]+det16*y[Nvar2+16]+
         det17*y[Nvar2+17]+det18*y[Nvar2+18]+det19*y[Nvar2+19]);


   dT2dt=0.0;


   dpgdt=Rg/(v-b)*dTdt+(-Rg*y[3]/pow(v-b,2.0)+2.0*a/pow(v,3.0))*dvdt+Rg*y[3]/pow(v-b,2.0)*dbdt-1.0/v/v*dadt;

   dpg2dt=Rg/(v-b)*dT2dt;

   drogdt=0.0;
   icont=1;

   while(icont<=NumEsp){
      drogdt=drogdt+Masa[icont]*dydt[Nvar2+icont];
      icont=icont+1;
   }

   drogdt=drogdt/Na/V-M/V/V*4.0*pi*y[1]*y[1]*y[2];



   dydt[1]=y[2];

   dydt[4]=4.0*pi*y[1]*y[1]*(mp*Na/Mh2o);

   dydt[5]=0.0;

   //***************************************************************************
   //************************calculo de los coeficientes del Chebyshev**********
   //***************************************************************************

   if((modelo == 3) || (modelo == 4)){

      Pp=dpgdt;

      f1=gama/(gama-1.0)*pg/y[3]/k;

      f2=2.0/y[1]/y[1];

      f3=Pp/(3.0*(gama-1.0)*k*y[3]);

      f4=1.0/(y[3]*k*y[1]*y[1]);

      f5=gama*pg*y[2]/(gama-1.0)/y[3]/k/y[1];

      f3=f3+f5;

      for(i=0;i<=NCh;i++)
         F[i]=0.0;

      //definicion del vector independiente

      j=0;
      while(j<=NCh){
         DPTCh=DerPCheb(colpoints[0],2*j);
         DSTCh=DerSCheb(colpoints[0],2*j);

         F[0]=F[0]+1.0/(f1*y[1]*y[1])*y[6+j]*DSTCh;
         j=j+1;
      }

      cont=1;
      aux=4.0;
      sum=0.0;
      while(cont<=NCh){
         sum=sum+aux*y[6+cont];
         aux=pow(-1,cont)*(fabs(aux)+(2*cont+1)*4.0);
         cont=cont+1;
      }

      F[0]=F[0]+Pp/f1+sum*f2/f1;

      for(i=1;i<NCh;i++){
         j=0;
         while(j<=NCh){
            DPTCh=DerPCheb(colpoints[i],2*j);
            DSTCh=DerSCheb(colpoints[i],2*j);

            F[i]=F[i]+1.0/(f1*y[1]*y[1])*y[6+j]*DSTCh+f2/f1/colpoints[i]*y[6+j]*DPTCh+
            f3/f1*colpoints[i]*y[6+j]*DPTCh;

            j=j+1;

         }

         cont=0;
         sum=0.0;
         while(cont<=NCh){
            DPTCh=DerPCheb(colpoints[i],2*cont);
            sum=sum+y[6+cont]*DPTCh;
            cont=cont+1;
         }

         F[i]=F[i]-f4/f1*sum*sum+Pp/f1;


      }

      F[NCh]=0.0;

      for(i=0;i<=NCh;i++)
         x[i]=0.0;

      for(i=0;i<=NCh;i++){
         j=0;
         while(j<=NCh){
            x[i]=x[i]+Q[j][i]*F[j];
            j=j+1;
         }
      }

      for(i=0;i<=NCh;i++)
         dydt[6+i]=0.0;

      dydt[6+NCh]=x[NCh]/R[NCh][NCh];

      for(i=(NCh-1);i>=0;i--){
         sum=0;
         j=i+1;
         while(j<=(NCh)){
            sum=sum+R[i][j]*dydt[j+6];
            j=j+1;
         }
         dydt[i+6]=(x[i]-sum)/R[i][i];
      }
   } //fin del if modelo igual al de Hao y Prosperetti




   //***************************************************************************
   //*******************los mpp*************************************************
   //***************************************************************************

   if(alfaM==0){
      mpp1=0.0;
      mpp2=0.0;
   }
   else{
      pvap=pvap0;
      dpvapdt=0.0;
      Conductividad1(k1);
      Conductividad2(k2);

      dk1dt=0.0;

      icont=1;

      while(icont<=NumEsp){
	 dk1dt=dk1dt+k2[icont]*y[Nvar2+icont]/ntot*dTdt-(k1[icont]+k2[icont]*y[3])*y[Nvar2+icont]/ntot/ntot*dntdt+
	 (k1[icont]+k2[icont]*y[3])*dydt[Nvar2+icont]/ntot;
	 icont=icont+1;
      }

      dk2dt=0.0;

      icont=1;

      while(icont<=NumEsp){
         dk2dt=dk2dt+k2[icont]*y[Nvar2+icont]/ntot*dT2dt;
         icont=icont+1;
      }

      dmdt=0.0;

      icont=1;

      while(icont<=NumEsp){
         dmdt=dmdt+(Masa[icont]*dydt[Nvar2+icont])/ntot/Na;
         icont=icont+1;
      }

      dmdt=dmdt-m/ntot*dntdt;


      dnpdt=dntdt/V-4.0*3.14159*ntot*y[1]*y[1]*y[2]/V/V;

      dlandadt=4.0*3.14159*y[1]*y[1]*y[2]/(sqrt(2.0)*sigmap*ntot)-
			 V/(sqrt(2.0)*sigmap*ntot*ntot)*dntdt;


      if(modelo == 1){
         dqdTB=k/n/landa;
         dq1dt=dk1dt*(TB-y[3])/n/landa-k/n/landa*dTdt-k*(TB-y[3])/n/landa/landa*dlandadt;
         dq2dt=dk2dt*(TB-y[3])/n/landa-k/n/landa*dT2dt;
      }



      if(y[2] == 0){
         dlth1dt=y[2]/pi;
         dlth2dt=0.0;
      }
      else{
         if((y[1]/pi)<(sqrt(y[1]*Xi/fabs(y[2])))){
            dlth1dt=y[2]/pi;
            dlth2dt=0.0;
         }
         else{
            dlth1dt=0.5*sqrt(Xi/y[1]/y[2])*y[2];
            dlth2dt=-0.5*sqrt(y[1]*Xi/y[2]/y[2]/y[2]);
         }

      }

      if(modelo == 2){
         dqdTB=0.0;
         dq1dt=dk1dt*(Tinf-y[3])/lth-k/lth*dTdt-k/lth/lth*(Tinf-y[3])*dlth1dt;
         dq2dt=dk2dt*(Tinf-y[3])/lth-k/lth*dT2dt-k/lth/lth*(Tinf-y[3])*dlth2dt;
      }

      if(modelo == 3){
         dqdTB=0.0;
         dq1dt=0.0;
         dq2dt=0.0;
         j=NCh;

         while(j>=0){
            DPTCh=DerPCheb(1,2*j);
            dq1dt=dq1dt-1.0/y[1]/y[1]*y[2]*y[6+j]*DPTCh+1.0/y[1]*dydt[6+j]*DPTCh;
            j=j-1;
         }

      }

      if(modelo == 4){

         dqdTB=0.0;
         dq1dt=0.0;
         dq2dt=0.0;
         j=NCh;

         while(j>=2){
            DPTCh=DerPCheb(1,2*j);
            dq1dt=dq1dt-1.0/y[1]/y[1]*y[2]*y[6+j]*DPTCh+1.0/y[1]*dydt[6+j]*DPTCh;
            j=j-1;
         }

         dq1dt=dq1dt+dk1dt*(Tinf-y[3])/lth-k/lth*dTdt-k/lth/lth*(Tinf-y[3])*dlth1dt;
         dq2dt=dq2dt+dk2dt*(Tinf-y[3])/lth-k/lth*dT2dt-k/lth/lth*(Tinf-y[3])*dlth2dt;

      }


      dTB1dt=(dydt[5]-(2.0-ap*alfae)/2.0/K/alfae*sqrt(pi/2.0/K)*(-1.0/np/np*dnpdt*sqrt(m/TB)*qpp+
			1.0/2.0/np*sqrt(1.0/m/TB)*dmdt*qpp+1.0/np*sqrt(m/TB)*dq1dt))/
			(1.0+(2.0-ap*alfae)/2.0/K/alfae*sqrt(pi/2.0/K)*
			(-1.0/2.0/np*sqrt(m/TB/TB/TB)*qpp+1.0/np*sqrt(m/TB)*dqdTB));

      dTB2dt=-(2.0-ap*alfae)/2.0/K/alfae*sqrt(pi/2.0/K)/np*sqrt(m/TB)*dq2dt/
		  (1.0+(2.0-ap*alfae)/2.0/K/alfae*sqrt(pi/2.0/K)*
			(-1.0/2.0/np*sqrt(m/TB/TB/TB)*qpp+1.0/np*sqrt(m/TB)*dqdTB));

      dpv1dt=(1.0-y[Nvar2+6]/ntot)*pg/ntot*dntdt+y[Nvar2+6]/ntot*dpgdt;

      dpv2dt=y[Nvar2+6]/ntot*dpg2dt;



      if(rapidez ==1){
         mpp1=alfaM/sqrt(2.0*pi*Rv)*(dpvapdt/sqrt(y[5])-0.5*pvap/pow(y[5],1.5)*dydt[5]-
			1.0/sqrt(TB)*dpv1dt+0.5*pv/pow(TB,1.5)*dTB1dt);

         mpp2=alfaM/sqrt(2.0*pi*Rv)*(-1.0/sqrt(TB)*dpv2dt+0.5*pv/pow(TB,1.5)*dTB2dt);
      }
      else{
         Omega=mp/pv*sqrt(Rv*y[3]/2.0);

         Gama=exp(-Omega*Omega)-Omega*sqrt(3.14159)*(1.0-erf(Omega,N));

         mpp1=alfaM/sqrt(2.0*pi*Rv)*(dpvapdt/sqrt(y[5])-0.5*pvap/pow(y[5],1.5)*dydt[5]-
			mp/pv*(1.0-erf(Omega,N))*dpv1dt*sqrt(pi*Rv*y[3]/2.0/TB)+
			0.5*mp*(1.0-erf(Omega,N))*sqrt(pi*Rv/2.0/TB/y[3])*dTdt-
			Gama/sqrt(TB)*dpv1dt+0.5*Gama*pv/pow(TB,1.5)*dTB1dt)/
			(1.0-alfaM/2.0*(1.0-erf(Omega,N))*sqrt(y[3]/TB));

         mpp2=alfaM/sqrt(2.0*pi*Rv)*(-mp/pv*(1.0-erf(Omega,N))*dpv2dt*sqrt(pi*Rv*y[3]/
			2.0/TB)+0.5*mp*(1.0-erf(Omega,N))*sqrt(pi*Rv/2.0/TB/y[3])*dT2dt-
			Gama/sqrt(TB)*dpv2dt+0.5*Gama*pv/pow(TB,1.5)*dTB2dt)/
			(1.0-alfaM/2.0*(1.0-erf(Omega,N))*sqrt(y[3]/TB));
      }

   }
   //***************************************************************************
   //***************************************************************************
   //***************************************************************************



   dpB1dt=dpgdt+2.0*sigma*y[2]/y[1]/y[1]+4.0*mu*y[2]/y[1]/y[1]*(y[2]-mp/rol)+
		  4.0*mu/y[1]/rol*mpp1-2.0*mp*mpp1*(1.0/rol-1.0/rog)-mp*mp/rog/rog*drogdt;

   dpB2dt=dpg2dt-4.0*mu/y[1]+4.0*mu/y[1]/rol*mpp2-2.0*mp*mpp2*(1.0/rol-1.0/rog);


   ps=Amp*sin(w*(t+y[1]/cl)-fi); //presion acustica

   dydt[4]=dpgdt;

   dydt[2]=(1.0/rol*(1.0+y[2]/cl)*(pB-ps-Pinf)+mp/rol*(y[2]+mp/2.0/rol+mp*y[2]/2.0/cl/rol)-
            3.0/2.0*y[2]*y[2]*(1.0-1.0/3.0*y[2]/cl+2.0*mp/3.0/cl/rol)+y[1]/cl/rol*dpB1dt+
            mpp1*y[1]/rol*(1.0-y[2]/cl+mp/cl/rol))/(y[1]*(1.0-y[2]/cl+mp/cl/rol)-
            mpp2*y[1]/rol*(1.0-y[2]/cl+mp/cl/rol)-y[1]/cl/rol*dpB2dt);


   dydt[3]=dTdt+dT2dt*dydt[2];


   //Calcula las inestabilidades parametricas.
   dydt[NCh3+NumEsp]=y[NCh3+1+NumEsp];

   dydt[NCh3+1+NumEsp]=-y[NCh3+1+NumEsp]*(3.0*y[2]/y[1]+2.0*(nn+2.0)*(2.0*nn+1.0)*nu/y[1]/y[1])-
                       y[NCh3+NumEsp]*(nn-1.0)*((nn+1.0)*(nn+2.0)*sigma/rol/pow(y[1],3.0)+
                       2.0*(nn+2.0)*nu*y[2]/y[1]/y[1]/y[1]-dydt[2]/y[1])+
                       Ppert*sin(w*t)/rol/y[1];


   delete F;
   delete x;
   delete Masa;
   delete k1;
   delete k2;
   delete r;
   delete DiamCol;
   delete y0;
   delete fracc;
   delete InvDi;
   delete ldiffi;


}
