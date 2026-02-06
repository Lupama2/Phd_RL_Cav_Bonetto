#include <math.h>
//Algunas funciones

//************INICIO DE FUNCION QUE ENCUENTRA UN MINIMO*************************

double min(double arg1, double arg2){

   if(arg1 < arg2)
      return arg1;

   else
      return arg2;

}

//************INICIO DE FUNCION PARA CALCULAR EL FACTORIAL**********************

double fact (int n){

   double valor=1;
   int i;

   if (n==0){
      valor=1;
   }
   else{
      for (i=1;i<=n;i++) valor=valor*i;
   }
   return valor;
}

//************INICIO DE FUNCION PARA CALCULAR LA FUNCION ERROR******************
double erf(double omega,int N){

   double sum=0.0;
   int i;
   double valor;

   if (fabs(omega)<4){
      for(i=0;i<=N;i++)
         sum=sum+pow(-1.0,i)*pow(omega,2.0*i+1.0)/((2.0*i+1.0)*fact(i));
      valor=2.0/sqrt(3.14159)*sum;
   }
   else {
      if(omega>=4)
         valor=1.0;
      if(omega<=-4)
         valor=-1.0;
   }

   return valor;

}

//*********INICIO DE FUNCION QUE CALCULA EL mp**********************************

double calculamp(double mpmin,double mpmax,double epsilon, double pvap,
                 double pv,double T,double TLi,double TB,int N){

   //epsilon es la tolerancia

   extern double Rv,alfaM;
   double x,f1,f3;
   double Gama;
   double Omega;


   if(alfaM > 0){
      while ((mpmax-mpmin)>epsilon){

         //mp es la resta de los mp de evaporacion y de condensacion
         Omega=mpmax/pv*sqrt(Rv*T/2.0);
         Gama=exp(-Omega*Omega)-Omega*sqrt(3.14159)*(1.0-erf(Omega,N));

         f1=mpmax-alfaM/sqrt(2.0*3.14159*Rv)*(pvap/sqrt(TLi)-Gama*pv/sqrt(TB));

         x=(mpmin+mpmax)/2.0;

         Omega=x/pv*sqrt(Rv*T/2.0);
         Gama=exp(-Omega*Omega)-Omega*sqrt(3.14159)*(1.0-erf(Omega,N));

         f3=x-alfaM/sqrt(2.0*3.14159*Rv)*(pvap/sqrt(TLi)-Gama*pv/sqrt(TB));

         if(f1*f3>0){
	    mpmax=x;
         }
         else{
            mpmin=x;
         }
      }//fin del while

      return x;
      //devuelve el valor de mp como la raiz de la ecuacion trascendente con tolerancia epsilon


   } //fin del if

   else {
      x=0;
      return x;
   }


}

//*************INICIO DE FUNCION QUE DESPEJA TB*********************************

//metodo de Toegel
double calculaTBT(double Tmin,double Tmax,double epsilon,double TLi,double k,
                  double lth,double np,double T,double m){

 //Para hacer biseccion usa un Tmin y un Tmax, a partir de los cuales hace los calculos
 //epsilon es la tolerancia para encontrar la raiz por el metodo de biseccion
 //TLi es la temperatura en la pared de la burbuja del lado externo
 //k es la conductividad que depende de la temperatura
 //landa es una funcion que depende de la temperatura
 //np es la densidad de numero de moleculas
 //T es la temperatura en el interior de la burbuja

 extern double ap,alfae,K,n,Tinf; //El nodo parametros.cpp tiene que es c/u de estos valores

 double f1,f3,x;


 while((Tmax-Tmin)>epsilon){

  f1=Tmax-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*Tmax))*(2.0-ap*alfae)/alfae*k*(Tinf-T)/(lth);

  x=(Tmax+Tmin)/2.0;

  f3=x-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*x))*(2.0-ap*alfae)/alfae*k*(Tinf-T)/(lth);


  if(f1*f3>0){
	Tmax=x;
  }
  else{
	Tmin=x;
  }


 }

	return x;


}


//metodo de Yasui
double calculaTBY(double Tmin,double Tmax,double epsilon,double TLi,double k,
                  double landa,double np,double T,double m){

 //Para hacer biseccion usa un Tmin y un Tmax, a partir de los cuales hace los calculos
 //epsilon es la tolerancia para encontrar la raiz por el metodo de biseccion
 //TLi es la temperatura en la pared de la burbuja del lado externo
 //k es la conductividad que depende de la temperatura
 //landa es una funcion que depende de la temperatura
 //np es la densidad de numero de moleculas
 //T es la temperatura en el interior de la burbuja

 extern double ap,alfae,K,n,Tinf; //El nodo parametros.cpp tiene que es c/u de estos valores

 double f1,f3,x;


 while((Tmax-Tmin)>epsilon){

  f1=Tmax-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*Tmax))*(2.0-ap*alfae)/alfae*k*(Tmax-T)/(n*landa);

  x=(Tmax+Tmin)/2.0;

  f3=x-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*x))*(2.0-ap*alfae)/alfae*k*(x-T)/(n*landa);


  if(f1*f3>0){
	Tmax=x;
  }
  else{
	Tmin=x;
  }


 }

	return x;


}

//metodo de Ha0-Prosperetti
double calculaTBHP(double Tmin,double Tmax,double epsilon,double TLi,double qpp,
                   double np,double m){

 //Para hacer biseccion usa un Tmin y un Tmax, a partir de los cuales hace los calculos
 //epsilon es la tolerancia para encontrar la raiz por el metodo de biseccion
 //TLi es la temperatura en la pared de la burbuja del lado externo
 //k es la conductividad que depende de la temperatura
 //landa es una funcion que depende de la temperatura
 //np es la densidad de numero de moleculas
 //T es la temperatura en el interior de la burbuja

 extern double ap,alfae,K,n,Tinf; //El nodo parametros.cpp tiene que es c/u de estos valores

 double f1,f3,x;


 while((Tmax-Tmin)>epsilon){

  f1=Tmax-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*Tmax))*(2.0-ap*alfae)/alfae*qpp;

  x=(Tmax+Tmin)/2.0;

  f3=x-TLi+1.0/(2.0*K*np)*sqrt(3.14159*m/(2.0*K*x))*(2.0-ap*alfae)/alfae*qpp;


  if(f1*f3>0){
	Tmax=x;
  }
  else{
	Tmin=x;
  }


 }

	return x;


}

//***********funcion que calcula la concentracion de gases**********************

double calculaCinfCsat2(double Cmin,double Cmax,double epsilon,double f1,
                        double f2,double f3,double f4,double f5,double f6,
                        double f7,double f8, double f9,double nt,double porcgn,
                        double porcn2,double porco2){

   double g1,g3,x;
   double nAr,nO2,nN2;
   double f10,f11,f12,f13,f14,f15;

   while(fabs(Cmax-Cmin)>epsilon){

      nAr=porcgn*Cmax/f7;

      f10=f2-f3;
      f11=(f1-f3)*nAr+f3*nt-f8*f7;
      f12=porco2*f8*Cmax;
      f13=f6-f5;
      f14=(f4-f5)*nAr+f5*nt-f9*f7;
      f15=porcn2*f9*Cmax;

      if(fabs(f8*f7) > 10.0*fabs(f10))
         nO2=-f12/f11;
      else
         nO2=(-f11-sqrt(f11*f11-4.0*f10*f12))/(2.0*f10);


      if(fabs(f9*f7) > 10.0*fabs(f13))
         nN2=-f15/f14;
      else
         nN2=(-f14-sqrt(f14*f14-4.0*f13*f15))/(2.0*f13);

      g1=nt-nAr-nO2-nN2;

      x=(Cmax+Cmin)/2.0;

      nAr=porcgn*x/f7;

      f10=f2-f3;
      f11=(f1-f3)*nAr+f3*nt-f8*f7;
      f12=porco2*f8*x;
      f13=f6-f5;
      f14=(f4-f5)*nAr+f5*nt-f9*f7;
      f15=porcn2*f9*x;

      if(fabs(f8*f7) > 10.0*fabs(f10))
         nO2=-f12/f11;
      else
         nO2=(-f11-sqrt(f11*f11-4.0*f10*f12))/(2.0*f10);

      if(fabs(f9*f7)>10.0*fabs(f13))
         nN2=-f15/f14;
      else
         nN2=(-f14-sqrt(f14*f14-4.0*f13*f15))/(2.0*f13);

      g3=nt-nAr-nO2-nN2;

      if(g1*g3>0){
         Cmax=x;
      }
      else{
         Cmin=x;
      }

   }
   return x;

}

//funcion que retorna el polinomio de Chebyshev, evaluado en y, y de orden NCh
double PolCheb(double y, int NCh){

   double TCh;
   double Tnm2, Tnm1;
   int cont;

   if(NCh == 0)
      TCh=1.0;

   if(NCh == 1)
      TCh=y;
   else{
       cont=2;
       Tnm2=1.0;
       Tnm1=y;
       while(cont<=NCh){

          TCh=2.0*y*Tnm1-Tnm2;
          Tnm2=Tnm1;
          Tnm1=TCh;
          cont=cont+1;

       } //fin while


   }//fin else

   return TCh;
} //fin funcion



//funcion que retorna la derivada primera del polinomio de Chebyshev, evaluado en y, y de orden NCh

double DerPCheb(double y, int NCh){

   double DPTCh;
   double TCh;
   double dTnm1,dTnm2;
   int cont;

    if(NCh == 0)
       DPTCh=0.0;

    if(NCh == 1)
       DPTCh=1.0;
    else{
       cont=2;
       dTnm2=0.0;
       dTnm1=1.0;

       while(cont<=NCh){
          TCh=PolCheb(y,cont-1);
          DPTCh=2.0*TCh+2.0*y*dTnm1-dTnm2;
          dTnm2=dTnm1;
          dTnm1=DPTCh;
          cont=cont+1;
       }

    }//fin else

    return DPTCh;

} //fin funcion derivada


//funcion que retorna la derivada segunda del polinomio de Chebyshev de orden NCh
double DerSCheb(double y, int NCh){

   double DSTCh;
   double DPTCh;
   double dsTnm1, dsTnm2;
   int cont;

   if(NCh == 0)
      DSTCh=0.0;

   if(NCh == 1)
      DSTCh=0.0;
   else{
      cont=2;
      dsTnm1=0.0;
      dsTnm2=0.0;

      while(cont<=NCh){

         DPTCh=DerPCheb(y,cont-1);
         DSTCh=4.0*DPTCh+2.0*y*dsTnm1-dsTnm2;
         dsTnm2=dsTnm1;
         dsTnm1=DSTCh;
         cont=cont+1;

      }


   }

   return DSTCh;

} //fin funcion derivada segunda



//*********funcion que asigna la masa molar a cada especie**********************
void MasaMolar(double *Masa){

   extern double Mh2o,MAr,MO2,MN2,MH2,MH,MO,MOH,MN,MNH,MNH2,MNH3,MN2H,MNO,MNO2,MN2O,MHNO;

   Masa[1]=MH2;
   Masa[2]=MH;
   Masa[3]=MO;
   Masa[4]=MO2;
   Masa[5]=MOH;
   Masa[6]=Mh2o;
   Masa[7]=MN;
   Masa[8]=MNH;
   Masa[9]=MNH2;
   Masa[10]=MNH3;
   Masa[11]=MN2H;
   Masa[12]=MNO;
   Masa[13]=MNO2;
   Masa[14]=MN2O;
   Masa[15]=MHNO;
   Masa[16]=MN2;
   Masa[17]=MAr;
   Masa[18]=0.034;
   Masa[19]=0.033;
}

//funcion que asigna el termino independiente de la conductividad de cada especie
//k=k1+k2*T

void Conductividad1(double *k1){

   extern double k1H2,k1H,k1O,k1O2,k1OH,k1h2o,k1N,k1NH,k1NH2,k1NH3,k1N2H,k1NO,k1NO2,k1N2O,k1HNO,k1N2,k1Ar;

   k1[1]=k1H2;
   k1[2]=k1H;
   k1[3]=k1O;
   k1[4]=k1O2;
   k1[5]=k1OH;
   k1[6]=k1h2o;
   k1[7]=k1N;
   k1[8]=k1NH;
   k1[9]=k1NH2;
   k1[10]=k1NH3;
   k1[11]=k1N2H;
   k1[12]=k1NO;
   k1[13]=k1NO2;
   k1[14]=k1N2H;
   k1[15]=k1HNO;
   k1[16]=k1N2;
   k1[17]=k1Ar;
   k1[18]=k1h2o;
   k1[19]=k1h2o;
}

//funcion que asigna el termino lineal de la conductividad de cada especie
void Conductividad2(double *k2){

   extern double k2H2,k2H,k2O,k2O2,k2OH,k2h2o,k2N,k2NH,k2NH2,k2NH3,k2N2H,k2NO,k2NO2,k2N2O,k2HNO,k2N2,k2Ar;

   k2[1]=k2H2;
   k2[2]=k2H;
   k2[3]=k2O;
   k2[4]=k2O2;
   k2[5]=k2OH;
   k2[6]=k2h2o;
   k2[7]=k2N;
   k2[8]=k2NH;
   k2[9]=k2NH2;
   k2[10]=k2NH3;
   k2[11]=k2N2H;
   k2[12]=k2NO;
   k2[13]=k2NO2;
   k2[14]=k2N2H;
   k2[15]=k2HNO;
   k2[16]=k2N2;
   k2[17]=k2Ar;
   k2[18]=k2h2o;
   k2[19]=k2h2o;

}


void DiamColVect(double *DiamCol){

   DiamCol[1]=2.97e-10;
   DiamCol[2]=3.5e-10;
   DiamCol[3]=3.5e-10;
   DiamCol[4]=3.43e-10;
   DiamCol[5]=3.5e-10;
   DiamCol[6]=2.65e-10;
   DiamCol[7]=3.5e-10;
   DiamCol[8]=3.5e-10;
   DiamCol[9]=3.5e-10;
   DiamCol[10]=2.6e-10;
   DiamCol[11]=3.5e-10;
   DiamCol[12]=3.47e-10;
   DiamCol[13]=3.5e-10;
   DiamCol[14]=3.88e-10;
   DiamCol[15]=3.5e-10;
   DiamCol[16]=3.68e-10;
   DiamCol[17]=3.42e-10;
   DiamCol[18]=3.5e-10;
   DiamCol[19]=3.5e-10;

}

void reacciones(double *y,double VolExcl,double ntot,double V,
					 int Nvar,double *r){

   extern double Na;
   extern double t1,t2,t5,t6,t15,t16,t33,t34,t39;


   r[0]=pow(VolExcl,t1)*1.2e17*(1.0e-12)*(ntot/V)*(y[Nvar+3]/V)*(y[Nvar+3]/V)*
   pow(y[3],-1.0)*exp(-0.0/y[3])/Na/Na;

   r[1]=3.16e19*(1.0e-6)*(ntot/V)*(y[Nvar+4]/V)*pow(y[3],-1.3)*exp(-59893.0/y[3])/Na;

   r[2]=pow(VolExcl,t2)*5.0e17*(1.0e-12)*(ntot/V)*(y[Nvar+3]/V)*(y[Nvar+2]/V)*
   pow(y[3],-1.0)*exp(-0.0/y[3])/Na/Na;

   r[3]=3.54e17*(1.0e-6)*(ntot/V)*(y[Nvar+5]/V)*pow(y[3],-0.9)*exp(-51217.0/y[3])/Na;

   r[4]=3.87e4*(1.0e-6)*(y[Nvar+3]/V)*(y[Nvar+1]/V)*pow(y[3],2.7)*exp(-3150.0/y[3])/Na;

   r[5]=1.79e4*(1.0e-6)*(y[Nvar+2]/V)*(y[Nvar+5]/V)*pow(y[3],2.7)*exp(-2200.0/y[3])/Na;

   r[6]=2.65e16*(1.0e-6)*(y[Nvar+2]/V)*(y[Nvar+4]/V)*pow(y[3],-0.7)*exp(-8576.0/y[3])/Na;

   r[7]=9.0e13*(1.0e-6)*(y[Nvar+3]/V)*(y[Nvar+5]/V)*pow(y[3],-0.3)*exp(83.0/y[3])/Na;

   r[8]=pow(VolExcl,t5)*1.0e18*(1.0e-12)*(ntot/V)*(y[Nvar+2]/V)*(y[Nvar+2]/V)*
   pow(y[3],-1.0)*exp(-0.0/y[3])/Na/Na;

   r[9]=7.46e17*(1.0e-6)*(ntot/V)*(y[Nvar+1]/V)*pow(y[3],-0.8)*exp(-52177.0/y[3])/Na;

   r[10]=pow(VolExcl,t6)*2.2e22*(1.0e-12)*(ntot/V)*(y[Nvar+2]/V)*(y[Nvar+5]/V)*
   pow(y[3],-2.0)*exp(-0.0/y[3])/Na/Na;

   r[11]=3.67e23*(1.0e-6)*(ntot/V)*(y[Nvar+6]/V)*pow(y[3],-2.0)*exp(-59980.0/y[3])/Na;

   r[12]=2.16e8*(1.0e-6)*(y[Nvar+5]/V)*(y[Nvar+1]/V)*pow(y[3],1.5)*exp(-1726.0/y[3])/Na;

   r[13]=5.2e9*(1.0e-6)*(y[Nvar+2]/V)*(y[Nvar+6]/V)*pow(y[3],1.3)*exp(-9529.0/y[3])/Na;

   r[14]=3.57e4*(1.0e-6)*(y[Nvar+5]/V)*(y[Nvar+5]/V)*pow(y[3],2.4)*exp(-1062.0/y[3])/Na;

   r[15]=1.74e6*(1.0e-6)*(y[Nvar+3]/V)*(y[Nvar+6]/V)*pow(y[3],2.2)*exp(-7693.0/y[3])/Na;

   r[16]=2.7e13*(1.0e-6)*(y[Nvar+7]/V)*(y[Nvar+12]/V)*pow(y[3],0.0)*exp(-179.0/y[3])/Na;

   r[17]=5.32e13*(1.0e-6)*(y[Nvar+16]/V)*(y[Nvar+3]/V)*pow(y[3],0.1)*exp(-37989.0/y[3])/Na;

   r[18]=9.0e9*(1.0e-6)*(y[Nvar+7]/V)*(y[Nvar+4]/V)*pow(y[3],1.0)*exp(-3271.0/y[3])/Na;

   r[19]=8.33e8*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+3]/V)*pow(y[3],1.1)*exp(-19115.0/y[3])/Na;

   r[20]=3.36e13*(1.0e-6)*(y[Nvar+7]/V)*(y[Nvar+5]/V)*pow(y[3],0.0)*exp(-194.0/y[3])/Na;

   r[21]=1.16e15*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+2]/V)*pow(y[3],-0.3)*exp(-24715.0/y[3])/Na;

   r[22]=1.4e12*(1.0e-6)*(y[Nvar+14]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(-5440.0/y[3])/Na;

   r[23]=2.05e8*(1.0e-6)*(y[Nvar+16]/V)*(y[Nvar+4]/V)*pow(y[3],1.0)*exp(-44858.0/y[3])/Na;

   r[24]=2.9e13*(1.0e-6)*(y[Nvar+14]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(-11651.0/y[3])/Na;

   r[25]=2.05e8*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+12]/V)*pow(y[3],1.0)*exp(-29105.0/y[3])/Na;

   r[26]=3.87e14*(1.0e-6)*(y[Nvar+14]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-9502.0/y[3])/Na;

   r[27]=1.55e8*(1.0e-6)*(y[Nvar+16]/V)*(y[Nvar+5]/V)*pow(y[3],1.4)*exp(-40259.0/y[3])/Na;

   r[28]=5.59e12*(1.0e-6)*(y[Nvar+14]/V)*(ntot/V)*pow(y[3],-0.7)*exp(-28480.0/y[3])/Na;

   r[29]=pow(VolExcl,t15)*6.48e6*(1.0e-12)*(ntot/V)*(y[Nvar+16]/V)*(y[Nvar+3]/V)*pow(y[3],0.5)*exp(-8002.0/y[3])/Na/Na;

   r[30]=pow(VolExcl,t16)*1.06e20*(1.0e-12)*(ntot/V)*(y[Nvar+12]/V)*(y[Nvar+3]/V)*
   pow(y[3],-1.4)*exp(-0.0/y[3])/Na/Na;

   r[31]=4.98e24*(1.0e-6)*(ntot/V)*(y[Nvar+13]/V)*pow(y[3],-2.2)*exp(-37005.0/y[3])/Na;

   r[32]=3.9e12*(1.0e-6)*(y[Nvar+13]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(121.0/y[3])/Na;

   r[33]=9.73e9*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+4]/V)*pow(y[3],0.6)*exp(-22763.0/y[3])/Na;

   r[34]=1.32e14*(1.0e-6)*(y[Nvar+13]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-181.0/y[3])/Na;

   r[35]=1.85e7*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+5]/V)*pow(y[3],0.9)*exp(-14397.0/y[3])/Na;

   r[36]=4.0e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(-0.0/y[3])/Na;

   r[37]=1.25e15*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+2]/V)*pow(y[3],-0.2)*exp(-35783.0/y[3])/Na;

   r[38]=3.2e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-166.0/y[3])/Na;

   r[39]=6.36e13*(1.0e-6)*(y[Nvar+7]/V)*(y[Nvar+1]/V)*pow(y[3],0.1)*exp(-12388.0/y[3])/Na;

   r[40]=2.0e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+5]/V)*pow(y[3],0.0)*exp(-0.0/y[3])/Na;

   r[41]=6.52e17*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+2]/V)*pow(y[3],-0.9)*exp(-9192.0/y[3])/Na;

   r[42]=2.0e9*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+5]/V)*pow(y[3],1.2)*exp(-0.0/y[3])/Na;

   r[43]=8.93e10*(1.0e-6)*(y[Nvar+7]/V)*(y[Nvar+6]/V)*pow(y[3],1.1)*exp(-20029.0/y[3])/Na;

   r[44]=4.61e5*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+4]/V)*pow(y[3],2.0)*exp(-3271.0/y[3])/Na;

   r[45]=5.41e7*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+3]/V)*pow(y[3],1.1)*exp(-8706.0/y[3])/Na;

   r[46]=1.28e06*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+4]/V)*pow(y[3],1.5)*exp(-50.0/y[3])/Na;

   r[47]=1.08e5*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+5]/V)*pow(y[3],1.7)*exp(-27164.0/y[3])/Na;

   r[48]=1.5e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+7]/V)*pow(y[3],0.0)*exp(-0.0/y[3])/Na;

   r[49]=9.16e14*(1.0e-6)*(y[Nvar+16]/V)*(y[Nvar+2]/V)*pow(y[3],-0.1)*exp(-73602.0/y[3])/Na;

   r[50]=2.0e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+6]/V)*pow(y[3],0.0)*exp(-6970.0/y[3])/Na;

   r[51]=2.92e16*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+1]/V)*pow(y[3],-0.7)*exp(-8361.0/y[3])/Na;

   r[52]=2.16e13*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+12]/V)*pow(y[3],-0.2)*exp(-0.0/y[3])/Na;

   r[53]=3.1e13*(1.0e-6)*(y[Nvar+16]/V)*(y[Nvar+5]/V)*pow(y[3],0.0)*exp(-49091.0/y[3])/Na;

   r[54]=3.65e14*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+12]/V)*pow(y[3],-0.5)*exp(-0.0/y[3])/Na;

   r[55]=1.13e21*(1.0e-6)*(y[Nvar+14]/V)*(y[Nvar+2]/V)*pow(y[3],-1.6)*exp(-18337.0/y[3])/Na;

   r[56]=3.0e12*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(-0.0/y[3])/Na;

   r[57]=2.95e15*(1.0e-6)*(y[Nvar+5]/V)*(y[Nvar+8]/V)*pow(y[3],-0.6)*exp(-14365.0/y[3])/Na;

   r[58]=3.9e13*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(-0.0/y[3])/Na;

   r[59]=3.82e16*(1.0e-6)*(y[Nvar+2]/V)*(y[Nvar+15]/V)*pow(y[3],-0.6)*exp(-14360.0/y[3])/Na;

   r[60]=4.0e13*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-1837.0/y[3])/Na;

   r[61]=1.27e12*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+1]/V)*pow(y[3],0.4)*exp(-7970.0/y[3])/Na;

   r[62]=9.0e7*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+5]/V)*pow(y[3],1.5)*exp(232.0/y[3])/Na;

   r[63]=6.38e7*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+6]/V)*pow(y[3],1.7)*exp(-13698.0/y[3])/Na;

   r[64]=0.0;

   r[65]=0.0;

   r[66]=1.3e14*(1.0e-6)*(y[Nvar+11]/V)*(ntot/V)*pow(y[3],-0.1)*exp(-2506.0/y[3])/Na;

   r[67]=pow(VolExcl,t34)*1.5e13*(1.0e-12)*(ntot/V)*(y[Nvar+16]/V)*(y[Nvar+2]/V)*pow(y[3],0.3)*exp(-6256.0/y[3])/Na/Na;

   r[68]=2.5e13*(1.0e-6)*(y[Nvar+11]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(0.0/y[3])/Na;

   r[69]=1.07e12*(1.0e-6)*(y[Nvar+5]/V)*(y[Nvar+16]/V)*pow(y[3],0.6)*exp(-54976.0/y[3])/Na;

   r[70]=7.0e13*(1.0e-6)*(y[Nvar+11]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(0.0/y[3])/Na;

   r[71]=1.67e12*(1.0e-6)*(y[Nvar+8]/V)*(y[Nvar+12]/V)*pow(y[3],0.4)*exp(-5883.0/y[3])/Na;

   r[72]=5.0e13*(1.0e-6)*(y[Nvar+11]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(0.0/y[3])/Na;

   r[73]=4.61e12*(1.0e-6)*(y[Nvar+1]/V)*(y[Nvar+16]/V)*pow(y[3],0.6)*exp(-55923.0/y[3])/Na;

   r[74]=2.0e13*(1.0e-6)*(y[Nvar+11]/V)*(y[Nvar+5]/V)*pow(y[3],0.0)*exp(0.0/y[3])/Na;

   r[75]=4.14e13*(1.0e-6)*(y[Nvar+6]/V)*(y[Nvar+16]/V)*pow(y[3],0.4)*exp(-63728.0/y[3])/Na;

   r[76]=pow(VolExcl,t39)*4.48e19*(1.0e-12)*(ntot/V)*(y[Nvar+2]/V)*(y[Nvar+12]/V)*
   pow(y[3],-1.3)*exp(-372.0/y[3])/Na/Na;

   r[77]=2.87e22*(1.0e-6)*(ntot/V)*(y[Nvar+15]/V)*pow(y[3],-1.9)*exp(-25004.0/y[3])/Na;

   r[78]=2.15e13*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+3]/V)*pow(y[3],0.0)*exp(0.0/y[3])/Na;

   r[79]=2.4e10*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+5]/V)*pow(y[3],0.7)*exp(-26593.0/y[3])/Na;

   r[80]=9.0e11*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+2]/V)*pow(y[3],0.7)*exp(-332.0/y[3])/Na;

   r[81]=1.05e9*(1.0e-6)*(y[Nvar+1]/V)*(y[Nvar+12]/V)*pow(y[3],1.5)*exp(-27878.0/y[3])/Na;

   r[82]=1.3e7*(1.0e-6)*(y[Nvar+15]/V)*(y[Nvar+5]/V)*pow(y[3],1.9)*exp(478.0/y[3])/Na;

   r[83]=2.91e5*(1.0e-6)*(y[Nvar+12]/V)*(y[Nvar+6]/V)*pow(y[3],2.5)*exp(-24867.0/y[3])/Na;

   r[84]=5.4e5*(1.0e-6)*(y[Nvar+10]/V)*(y[Nvar+2]/V)*pow(y[3],2.4)*exp(-4990.0/y[3])/Na;

   r[85]=1.07e3*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+1]/V)*pow(y[3],2.9)*exp(-2257.0/y[3])/Na;

   r[86]=5.0e7*(1.0e-6)*(y[Nvar+10]/V)*(y[Nvar+5]/V)*pow(y[3],1.6)*exp(-481.0/y[3])/Na;

   r[87]=4.56e6*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+6]/V)*pow(y[3],1.8)*exp(-5543.0/y[3])/Na;

   r[88]=9.4e6*(1.0e-6)*(y[Nvar+10]/V)*(y[Nvar+3]/V)*pow(y[3],1.9)*exp(-3251.0/y[3])/Na;

   r[89]=2.37e4*(1.0e-6)*(y[Nvar+9]/V)*(y[Nvar+5]/V)*pow(y[3],2.3)*exp(435.0/y[3])/Na;

   r[90]=3.0e12*(1.0e-6)*(y[Nvar+19]/V)*(y[Nvar+19]/V)*pow(y[3],0.0)*exp(-700.0/y[3])/Na;

   r[91]=4.53e14*(1.0e-6)*(y[Nvar+18]/V)*(y[Nvar+4]/V)*pow(y[3],-0.39)*exp(-19700.0/y[3])/Na;

   r[92]=1.2e17*(1.0e-6)*(y[Nvar+18]/V)*(ntot/V)*pow(y[3],0.0)*exp(-22900.0/y[3])/Na;

   r[93]=9.0e5*(1.0e-12)*(ntot/V)*(y[Nvar+5]/V)*(y[Nvar+5]/V)*pow(y[3],0.9)*exp(3050.0/y[3])/Na/Na;

   r[94]=3.2e14*(1.0e-6)*(y[Nvar+18]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-4510.0/y[3])/Na;

   r[95]=1.14e9*(1.0e-6)*(y[Nvar+6]/V)*(y[Nvar+5]/V)*pow(y[3],1.36)*exp(-38180.0/y[3])/Na;

   r[96]=4.82e13*(1.0e-6)*(y[Nvar+18]/V)*(y[Nvar+2]/V)*pow(y[3],0.0)*exp(-4000.0/y[3])/Na;

   r[97]=1.41e11*(1.0e-6)*(y[Nvar+1]/V)*(y[Nvar+19]/V)*pow(y[3],0.66)*exp(-12320.0/y[3])/Na;

   r[98]=9.55e6*(1.0e-6)*(y[Nvar+18]/V)*(y[Nvar+3]/V)*pow(y[3],2.0)*exp(-2000.0/y[3])/Na;

   r[99]=4.62e3*(1.0e-6)*(y[Nvar+5]/V)*(y[Nvar+19]/V)*pow(y[3],2.75)*exp(-9277.0/y[3])/Na;

   r[100]=1.00e13*(1.0e-6)*(y[Nvar+18]/V)*(y[Nvar+5]/V)*pow(y[3],0.0)*exp(-900.0/y[3])/Na;

   r[101]=2.8e13*(1.0e-6)*(y[Nvar+6]/V)*(y[Nvar+19]/V)*pow(y[3],0.0)*exp(-16500.0/y[3])/Na;


}
//******************************************************************************

void fraccion(double *fracc){

   fracc[1]=5;
   fracc[2]=3;
   fracc[3]=3;
   fracc[4]=5;
   fracc[5]=5;
   fracc[6]=6;
   fracc[7]=3;
   fracc[8]=5;
   fracc[9]=6;
   fracc[10]=6;
   fracc[11]=5;
   fracc[12]=5;
   fracc[13]=6;
   fracc[14]=5;
   fracc[15]=5;
   fracc[16]=5;
   fracc[17]=3;
   fracc[18]=6;
   fracc[19]=6;
}

//********funcion que calcula el maximo*****************************************
long double Max(long double x1, long double x2){

   if(x1>=x2)
      return x1;
   else
      return x2;

}

//*************funcion que calcula el espectro**********************************
long double espectro(long double RADIUS,long double TEMPERATURE, long double NAR,
                     long double Z,long double LAMBDA,long double EION,
                     long double CTR, long double DTR, long double E2){


   return (4702.048198519318*pow(RADIUS,2.0))/((-1. + pow(2.7182818284590452354,
          1.4387752209501652/(LAMBDA*TEMPERATURE)))*pow(LAMBDA,5.0)) -
          2.351024099259659e15/((-1. + pow(2.7182818284590452354,
          1.4387752209501652/(LAMBDA*TEMPERATURE)))*pow(LAMBDA,5.0)*
          pow((0.0078997375581096*pow(2.7182818284590452354,
          (-7.242963696165495e18*EION)/TEMPERATURE + 1.4387752209501651e-6/
          (TEMPERATURE*Max(1.9864454404374119e-25/E2,10.e-7*LAMBDA)))*
          pow(LAMBDA,3.0)*NAR*TEMPERATURE)/pow(RADIUS,3.0) +
          (5.8839621648085375e-9*pow(LAMBDA,2.0)*NAR*sqrt(NAR/pow(RADIUS,3.0))*
          (CTR + (2.414321232055165e18*DTR)/TEMPERATURE)*
          pow(TEMPERATURE,2.25))/(pow(2.7182818284590452354,
          (3.6214818480827474e18*EION)/TEMPERATURE)*pow(RADIUS,3.0)),2.0)) +
          2.351024099259659e15/(pow(2.7182818284590452354,
          2.e-6*RADIUS*((0.0078997375581096*
          pow(2.7182818284590452354,(-7.242963696165495e18*EION)/TEMPERATURE +
          1.4387752209501651e-6/(TEMPERATURE*Max(1.9864454404374119e-25/
          E2,10.e-7*LAMBDA)))*pow(LAMBDA,3.0)*NAR*TEMPERATURE)/pow(RADIUS,3.0) +
          (5.8839621648085375e-9*pow(LAMBDA,2.0)*pow(NAR/pow(RADIUS,3.0),1.5)*
          (CTR + (2.414321232055165e18*DTR)/TEMPERATURE)*pow(TEMPERATURE,2.25))/
          pow(2.7182818284590452354,(3.6214818480827474e18*EION)/TEMPERATURE)))*
          (-1. + pow(2.7182818284590452354,
          1.4387752209501652/(LAMBDA*TEMPERATURE)))*pow(LAMBDA,5.0)*
          pow((0.0078997375581096*pow(2.7182818284590452354,
          (-7.242963696165495e18*EION)/TEMPERATURE +1.4387752209501651e-6/
          (TEMPERATURE*Max(1.9864454404374119e-25/E2,10.e-7*LAMBDA)))*
          pow(LAMBDA,3.0)*NAR*TEMPERATURE)/pow(RADIUS,3.0) +
          (5.8839621648085375e-9*pow(LAMBDA,2.0)*NAR*sqrt(NAR/pow(RADIUS,3.0))*
          (CTR + (2.414321232055165e18*DTR)/TEMPERATURE)*
          pow(TEMPERATURE,2.25))/(pow(2.7182818284590452354,
          (3.6214818480827474e18*EION)/TEMPERATURE)*pow(RADIUS,3.0)),2.0)) +
          (4.702048198519319e9*RADIUS)/(pow(2.7182818284590452354,
          2.e-6*RADIUS*((0.0078997375581096*pow(2.7182818284590452354,
          (-7.242963696165495e18*EION)/TEMPERATURE +1.4387752209501651e-6/
          (TEMPERATURE*Max(1.9864454404374119e-25/E2,10.e-7*LAMBDA)))*
          pow(LAMBDA,3.0)*NAR*TEMPERATURE)/pow(RADIUS,3.0) +
          (5.8839621648085375e-9*pow(LAMBDA,2.0)*
          pow(NAR/pow(RADIUS,3.0),1.5)*
          (CTR + (2.414321232055165e18*DTR)/TEMPERATURE)*pow(TEMPERATURE,2.25))/
          pow(2.7182818284590452354,(3.6214818480827474e18*EION)/TEMPERATURE)))*
          (-1. + pow(2.7182818284590452354,
          1.4387752209501652/(LAMBDA*TEMPERATURE)))*pow(LAMBDA,5.0)*
          ((0.0078997375581096*pow(2.7182818284590452354,
          (-7.242963696165495e18*EION)/TEMPERATURE +1.4387752209501651e-6/
          (TEMPERATURE*Max(1.9864454404374119e-25/E2,10.e-7*LAMBDA)))*
          pow(LAMBDA,3.0)*NAR*TEMPERATURE)/pow(RADIUS,3.0) +
          (5.8839621648085375e-9*pow(LAMBDA,2.0)*NAR*sqrt(NAR/pow(RADIUS,3.0))*
          (CTR + (2.414321232055165e18*DTR)/TEMPERATURE)*pow(TEMPERATURE,2.25))/
          (pow(2.7182818284590452354,
          (3.6214818480827474e18*EION)/TEMPERATURE)*pow(RADIUS,3.0))));

}   
