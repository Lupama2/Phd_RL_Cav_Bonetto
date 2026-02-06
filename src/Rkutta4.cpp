void derivada(double t,double *y,double TB,
              double mp,double k,double *colpoints,double **Q,
              double **R, double *dydt);


void rk4(double *y,double *dydx,double x,
         double h, double TB, double mp,double k,double *colpoints,
         double **Q, double **R,double *yout){


   extern int Nvar;

   int i;
   double xh,hh,h6;
   double *dym=new double[Nvar];
   double *dyt=new double[Nvar];
   double *yt=new double[Nvar];


   hh=h*0.5;

   h6=h/6.0;

   xh=x+hh;


   for(i=1;i<Nvar;i++)
      yt[i]=y[i]+hh*dydx[i];

   derivada(xh,yt,TB,mp, k,colpoints,Q,R,dyt);

   for(i=1;i<Nvar;i++)
      yt[i]=y[i]+hh*dyt[i];

   derivada(xh,yt,TB,mp,k,colpoints,Q,R,dym);

   for(i=1;i<Nvar;i++){
      yt[i]=y[i]+h*dym[i];
      dym[i] += dyt[i];
   }


   derivada(x+h,yt,TB,mp,k,colpoints,Q,R,dyt);


   for(i=1;i<Nvar;i++)
      yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);


   delete dym;
   delete dyt;
   delete yt;

}  