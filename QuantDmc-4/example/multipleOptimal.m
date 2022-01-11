% multipleOptimal - Simple example DMC with two optimal quantizers 
%
% Example using QUANTMIDMCMULTI, which produces multpile optimal quantizers, 
% if they exist, for the binary-input DMC.
%
% The DMC below has two optimal but distinct quantizers, it is Example 1 from:
%
%   B. Kurkoski and H. Yagi, "Concatenation of a discrete memoryless channel 
%   and a quantizer," in Proceedings of the IEEE Information Theory Workshop, 
%  (Cairo, Egypt), pp. 160-164, January 2010.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

P(1,:) = [0.001  0.01  0.02  0.04  0.2   0.729];
P(2,:) = [0.729  0.2   0.04  0.02  0.01  0.001];

[Q,mi] = quantBiDmcMulti(P,4);

fprintf('Channel:\n');
P

fprintf('This channel has %d optimal quantizers:\n',length(Q));
for ii = 1:length(Q)
    fprintf('Q{%d}=\n',ii);
    disp(Q{ii});
    fprintf('with mutual information %g.\n\n',mi(ii) );
end
