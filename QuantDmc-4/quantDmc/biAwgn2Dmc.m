function [P,Quantizer,Boundary] = BiAwgn2Dmc(var,quantStep)
%BIAWGN2DMC  Creates a two-input DMC from a quantized AWGN channel
%
% P = BIAWGN2DMC(VAR,M)
%
%   Uses AWGN noise with variance VAR (defaults to 0.4)
%   quantized to M levels between +2 and -2 (defaults to 20).
%
%   P is a 2-by-M matrix, where:
%      P(j,m) = Pr( Y=m | X=j )
%   for a DMC with inputs X and outputs Y.  
%
%
% P = BIAWGN2DMC(VAR,QUANT)
%
%   Instead uses QUANT bin centers, where QUANT is a vector
%
% If VAR is a 1x2 vector then VAR(1) is used for -1 and VAR(2)
% is used for +1.
%
% [P,QUANTIZER,BOUNDARY] = BIAWGN2DMC
%  QUANTIZER is the 1-by-M vector used for the quantizer bin centers
%  BOUNDARY is the 1-by-M+1 vector used for the quantizer bin edges
%
%
% QuantDMC is (c) 2010-2012 Brian Kurkoski
% Distributed under an MIT-like license; see the file LICENSE
%

input = [-1 1];

if nargin < 1
    var = 0.4;
end

if nargin < 2
    quantStep = 20;
end

if length(quantStep) == 1
    Quantizer = linspace(-2,2,quantStep);
else
    Quantizer = sort(quantStep(:)');
end

if length(var) == 1
    var = [var var];
end

Boundary = (Quantizer(1:end-1) + Quantizer(2:end))/2;
Boundary = [-Inf Boundary Inf];

for ii = 1:2
    for jj = 1:length(Boundary)-1;
        t = Qf( (Boundary(jj) - input(ii) ) / sqrt(var(ii)) ) ...
            - Qf( (Boundary(jj+1) - input(ii) ) / sqrt(var(ii)) );
        P(ii,jj) = t;
    end
end


%Gaussian Q function
function out = Qf(x)

out = 0.5 .* erfc(x ./ sqrt(2) );
