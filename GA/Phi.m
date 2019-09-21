function [y]=Phi(x)

alpha=-0.4527;beta=0.0218;gama=0.86;
  for kk=1:length(x)
if double( x(kk))<=10
y(kk)=exp(alpha*(x(kk)^gama)+beta);
  else
y(kk)=sqrt(pi/x(kk))*exp(-x(kk)/4)*(1-10/7/x(kk));
  end
  end