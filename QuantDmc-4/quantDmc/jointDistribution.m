function pxy = jointDistribution(pygx,px)
%jointDistribution - compute joint distribution from P(Y|X) and P(X)
%
%   If PYGX is a conditional probability distribution with rows summing to
%   1, and PX is a vector, then JOINTDISTRIBUTION(PYGX,PX) is the joint
%   distribution P(X,Y).
%
%   QuantDMC (c) Brian Kurkoski and contributors
%   Distributed under an MIT-like license; see the file LICENSE

pxy = repmat(px(:),1,size(pygx,2)) .* pygx;