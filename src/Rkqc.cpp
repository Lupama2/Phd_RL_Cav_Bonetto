#include <conio.h>
#include <stdio.h>
#include <math.h>
#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


void derivada(double t,double *y,double TB, double mp,double k,
				  double *colpoints,double **Q,double **R, double *dydt);


void rk4(double *y,double *dydx,double x,double h, double TB, double mp,double k,
			double *colpoints, double **Q, double **R,double *yout);



void rkqc(double *y,double *dydx,double *x,double htry,double eps,double *yscal,
			 double *hdid,double *hnext,double TB, double mp, double k,
			 double *colpoints,double **Q,double **R, double *yout){

   int i;

   extern int Nvar;
   extern double hmin;

   double xsav,hh,h,temp,errmax;


   double *dysav=new double[Nvar];
   double *ysav=new double[Nvar];
   double *ytemp=new double[Nvar];
   double *dydxp=new double[Nvar];


   xsav=(*x);

   for(i=1;i<Nvar;i++){
      ysav[i]=y[i];
      dysav[i]=dydx[i];
   }



   h=htry;

   while (h>hmin){

      hh=0.5*h;

      rk4(ysav,dysav,xsav,hh,TB,mp,k,colpoints,Q,R,ytemp);

      *x=xsav+hh;

      derivada(*x,ytemp,TB,mp,k,colpoints,Q,R,dydxp);

      rk4(ytemp,dydxp,*x,hh,TB,mp,k,colpoints,Q,R,yout);

      *x=xsav+h;

      if (*x == xsav) printf("Step size too small in routine RKQC \n");

      rk4(ysav,dysav,xsav,h,TB,mp,k,colpoints,Q,R,ytemp);

      errmax=0.0;

      for (i=1;i<Nvar;i++){
         ytemp[i]=yout[i]-ytemp[i];
         temp=fabs(ytemp[i]/yscal[i]);
         if (errmax < temp)
            errmax=temp;
      }	/*for*/

      errmax= errmax/eps;

      if (errmax <= 1.0){
         *hdid=h;
         *hnext=(errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);

         break;

      }	/*if*/

      h=SAFETY*h*exp(PSHRNK*log(errmax));

   }

   for(i=1;i<Nvar;i++)
      yout[i] = yout[i]+ytemp[i]*FCOR;


   delete dysav;
   delete ysav;
   delete ytemp;
   delete dydxp;

}

#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON

