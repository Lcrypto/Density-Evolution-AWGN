function [y] = phi(x);
if x<10;
    y = exp(-.4527*x^.86+.0218);
else;
    y = (1-10/7/x)*sqrt(pi/x)*exp(-x/4);
end;