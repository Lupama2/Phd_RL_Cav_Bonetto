%RT
%FJB
%02nov2020

clear all
clear global
clc

%1_t 2_radio 3_v 4_a 5_T 6_pg nH2 nH nO nO2 nOH nH2O nN nNH nNH2 nNH3 nN2H nNO nNO2 nN2O nHNO nN2 nAr nH2O2 nHO2 
r6file = 'bfa03.dat';
r6 = load(r6file);
rc = r6(:,2);
[rmin,imin]=min(rc);
tc = r6(:,1);
tc(imin);
rc = r6(:,2);
rc(imin);
rp = r6(:,3);
r2 = r6(:,4);
te = r6(:,5);
te(imin);
n=2;
j=0;
last=length(tc);
tspan = [tc(1) tc(imin)];
eps = 0.0025;
ri = r6(1,2);
ic = [eps*ri 0.];
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
%opts = odeset('Refine',10)
j=0;
ti=0.;
ratio = eps;
temp=r6(1,5);
save('savefile','ti','ratio','temp','j')
[t,y] = ode45(@(t,y) rtode(t,y,tc,rc,rp,r2,te,j), tspan, ic,opts);

loglog((tc(imin)*ones(size(t))-t)*1.e6,abs(y(:,1)./interp1(tc,rc,t)),'o')