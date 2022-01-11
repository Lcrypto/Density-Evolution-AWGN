% [Q,MI] = QuantBiDmc(P,K)
%
% Finds the optimal quantizer of DMC given by P, quantized to K values,
% in the sense of maximizing mutual information.
% P is a 2-by-M matrix, where:
%   P(j,m) = Pr( Y=m | X=j ),
% for a DMC with inputs X and outputs Y.  K is an integer, generally less
% than M. 
%
% Q is one of the optimal quantizers. Q(m,k) is a 1 if DMC output m is quantized
% to k, and otherwise is a 0.
%
% MI is the mutual information between the channel input and the quantizer
% output.
%
% QuantDMC is (c) 2010-2012 Brian Kurkoski
% Distributed under an MIT-like license; see the file LICENSE
%

error('Run the command MEXIFY to compile QUANTMIDMC before using it');