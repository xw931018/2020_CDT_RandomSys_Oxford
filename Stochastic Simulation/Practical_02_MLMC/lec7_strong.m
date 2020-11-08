%
% Test strong convergence of Euler method
% for solving stochastic o.d.e.'s
%
% Test problem:   dS   = r*S dt + sig*S dW
% Exact solution: S(1) = S(0)*exp((lambda-0.5*sig^2)+sig*W(1))
%

clear all; close all
randn('state',0)

addpath('..')

%
% problem parameters
%

r   = 0.05; 
sig = 0.5; 
T   = 1; 
S0  = 100;

%
% Monte Carlo simulation comparing to 2h simulation
%

M  = 1e+5;   % total number of Monte Carlo paths
M2 = 1e+4;   % number of paths at a time 

for p = 1:7
  N = 2^p;
  h = T/N;

  sum1 = 0;
  sum2 = 0;
  sum3 = 0;
  sum4 = 0;

  for m = 1:M2:M
    m2 = min(M2,M-m+1); 

    S  = S0*ones(1,m2);
    S2 = S0*ones(1,m2);
    W  = 0;

    for n = 1:N/2
      dW1 = sqrt(h)*randn(1,m2);
      S  = S.*(1+r*h+sig*dW1);    
      dW2 = sqrt(h)*randn(1,m2);
      S  = S.*(1+r*h+sig*dW2);

      S2 = S2.*(1+r*2*h+sig*(dW1+dW2));

      W = W + dW1 + dW2;
    end

    Se = S0.*exp((r-0.5*sig^2)*T+sig*W);

    del  = (Se-S).^2;
    sum1 = sum1 + sum(del);
    sum2 = sum2 + sum(del.^2);

    del  = (S2-S).^2;
    sum3 = sum3 + sum(del);
    sum4 = sum4 + sum(del.^2);
  end

  hh(p)   = h;

  Vd = sum1/M;
  sd = sqrt((sum2/M - (sum1/M)^2)/(M-1));
  err1(p) = sqrt(Vd);
  err2(p) = (0.5/sqrt(Vd)) * 3*sd;

  Vd = sum3/M;
  sd = sqrt((sum4/M - (sum3/M)^2)/(M-1));
  err3(p) = sqrt(Vd);
  err4(p) = (0.5/sqrt(Vd)) * 3*sd;
end

figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.8]; set(gcf,'pos',pos);
loglog(hh,abs(err1),'b-*',hh,err2,'r-*',hh,abs(err3),'b--o',hh,err4,'r--o')
title('Strong convergence -- difference from exact and 2h approximation')
xlabel('h'); ylabel('Error');
legend(' exact error',' MC error',' relative error',' MC error','Location','NorthWest')
axis([.005 0.5 0.01 50])
print('-deps2c','lec7c.eps')
