%
% function V = european_call(r,sigma,T,S,K,opt)
%
% Black-Scholes European call option solution
% as defined in equation (3.17) on page 48 of 
% The Mathematics of Financial Derivatives
% by Wilmott, Howison and Dewynne
%
% r     - interest rate
% sigma - volatility
% T     - time interval
% S     - asset value(s)
% K     - strike price(s)
% opt   - 'value', 'delta', 'gamma' or 'vega'
% V     - option value(s)
%

function V = european_call(r,sigma,T,S,K,opt)

if nargin ~= 6
  error('wrong number of arguments');
end

S  = max(1e-100,S);     % avoids problems with S=0
K  = max(1e-100,K);     % avoids problems with K=0

d1 = ( log(S) - log(K) + (r+0.5*sigma^2)*T ) / (sigma*sqrt(T));
d2 = ( log(S) - log(K) + (r-0.5*sigma^2)*T ) / (sigma*sqrt(T));

switch opt
  case 'value'
    V = S.*N(d1) - exp(-r*T)*K.*N(d2);

  case 'delta'
    V = N(d1);

  case 'gamma'
    V = exp(-0.5*d1.^2) ./ (sigma*sqrt(2*pi*T)*S);

  case 'vega'
    V =           S.*(exp(-0.5*d1.^2)./sqrt(2*pi)).*( sqrt(T)-d1/sigma) ...
      - exp(-r*T)*K.*(exp(-0.5*d2.^2)./sqrt(2*pi)).*(-sqrt(T)-d2/sigma);

  otherwise
    error('invalid value for opt -- must be ''value'', ''delta'', ''gamma'', ''vega''')
end


%
% Normal cumulative distribution function
%

function ncf = N(x)

%ncf = 0.5*(1+erf(x/sqrt(2)));

xr = real(x);
xi = imag(x);

if abs(xi)>1e-10
  error 'imag(x) too large in N(x)'
end

ncf = 0.5*(1+erf(xr/sqrt(2))) ...
    + i*xi.*exp(-0.5*xr.^2)/sqrt(2*pi);
