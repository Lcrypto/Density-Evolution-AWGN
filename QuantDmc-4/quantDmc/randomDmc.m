function P = randomDMC(J,M)
% randomDMC - genreates a random DMC
%
% P = RandomDMC(J,M) generates a randomly generated discrete memoryless 
% channel (DMC) with J inputs and M outputs. P is a J-by-M matrix, where:
%   P(j,m) = Pr(Y=m | X=j ),
% for a DMC with inputs X and outputs Y.
%
% If not specified, J defaults to 2, and M defaults to 3
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

if nargin < 2
    M = 3;
end

if nargin < 1
    J = 2;
end

P = rand(J,M);
for jj = 1:J
    P(jj,:) = P(jj,:) ./ sum( P(jj,:) );
end

%sort when J=2
if J == 2
    [~,t] = sort(P(1,:) ./ P(2,:));
    P=P(:,t);
end
    