function dydt = rtode(t,y,tc,rc,rp,r2,te,j)
r0p = interp1(tc,rc,t); % Interpolate the data set (r) at time t
r1p = interp1(tc,rp,t); % Interpolate the data set (rp) at time t
r2p = interp1(tc,r2,t); % Interpolate the data set (r2p) at time t
tem = interp1(tc,te,t); % Interpolate the data set (tem) at time t
load ('savefile')
j=j+1;
ti(j)=t;
ratio(j)=y(1)/r0p;
temp(j)=tem;
if abs(y(1))>r0p
     timebreak=interp1(ratio,ti,1.)
     tempbreak=interp1(ratio,temp,1.)
     stop
end
save('savefile','ti','ratio','temp','j')
n=2;
sigma=0.08;
rhol=1.863e3;
mul=0.14;
nul=mul/rhol;
an = (n-1)*r2p/r0p-(n-1)*(n+1)*(n+2)*sigma/(rhol*r0p^3)-2*nul*r1p/r0p^3*(...
    -(n-1)*(n+1)*(n+2)+n*(n-1)*(n+2)/(1+1/n) );
    
bn = 3*r1p/r0p+2*nul/r0p^2*(-(n-1)*(n+1)*(n+2)+n*(n+2)^2/(1+1/n));
dydt = [y(2); -bn*y(2) + an*y(1)]; % Evaluate ODE at time t
end
