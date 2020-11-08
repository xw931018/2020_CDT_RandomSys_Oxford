%
% Test weak convergence of Euler method
% for European call option
%
% Test problem:   dS   = r*S dt + sig*S dW
%

clear all; close all
randn('state',0)

addpath ('..')

%
% problem parameters and exact solution
%

r   = 0.05; 
sig = 0.5; 
T   = 1; 
S0  = 100; 
K   = 110;

Ve  = european_call(r,sig,T,S0,K,'value');

%
% Monte Carlo simulation comparing to exact solution
%

M  = 1e+7;   % total number of Monte Carlo paths
M2 = 1e+4;   % number of paths at a time 

for p = 1:6
  N = 2^p;
  h = T/N;

  sum1 = 0;
  sum2 = 0;

  for m = 1:M2:M
    m2 = min(M2,M-m+1); 

    S = S0*ones(1,m2);

    for n = 1:N
      dW = sqrt(h)*randn(1,m2);
      S  = S.*(1+r*h+sig*dW);    
    end

    P = exp(-r*T)*max(S-K,0);

    sum1 = sum1 + sum(P);
    sum2 = sum2 + sum(P.^2);
  end

  V  = sum1/M;
  sd = sqrt((sum2/M - (sum1/M)^2)/(M-1));

  hh(p)   = h;
  err1(p) = V-Ve;
  err2(p) = 3*sd;
end

figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.8]; set(gcf,'pos',pos);
loglog(hh,abs(err1),'b-*',hh,err2,'r-*')
title('Weak convergence -- comparison to exact solution')
xlabel('h'); ylabel('Error');
legend(' Weak error',' MC error','Location','NorthWest')
axis([.01 0.5 0.01 0.2])
print('-deps2c','lec7a.eps')


%
% Monte Carlo simulation comparing to 2h simulation
%

M  = 1e+7;   % total number of Monte Carlo paths
M2 = 1e+4;   % number of paths at a time 

for p = 1:7
  N = 2^p;
  h = T/N;

  sum1 = 0;
  sum2 = 0;

  for m = 1:M2:M
    m2 = min(M2,M-m+1); 

    S  = S0*ones(1,m2);
    S2 = S0*ones(1,m2);

    for n = 1:N/2
      dW1 = sqrt(h)*randn(1,m2);
      S  = S.*(1+r*h+sig*dW1);    
      dW2 = sqrt(h)*randn(1,m2);
      S  = S.*(1+r*h+sig*dW2);

      S2 = S2.*(1+r*2*h+sig*(dW1+dW2));
    end

    P  = exp(-r*T)*max(S-K,0);
    P2 = exp(-r*T)*max(S2-K,0);

    sum1 = sum1 + sum(P-P2);
    sum2 = sum2 + sum((P-P2).^2);
  end

  Vd = sum1/M;
  sd = sqrt((sum2/M - (sum1/M)^2)/(M-1));

  hh(p)   = h;
  err1(p) = Vd;
  err2(p) = 3*sd;
end

figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.8]; set(gcf,'pos',pos);
loglog(hh,abs(err1),'b-*',hh,err2,'r-*')
title('Weak convergence -- difference from 2h approximation')
xlabel('h'); ylabel('Error');
legend(' Weak error',' MC error','Location','NorthWest')
axis([.005 0.5 0.001 0.2])
print('-deps2c','lec7b.eps')
