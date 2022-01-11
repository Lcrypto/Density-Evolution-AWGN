function [I] = joint2MI(P)
%joint2MI - Compute mutual information 
%
%   P is a joint distribution between X and Y, then
%   JOINT2MI(P) gives the mutual information I(X;Y).
%
%   QuantDMC (c) Brian Kurkoski and contributors
%   Distributed under an MIT-like license; see the file LICENSE

P(find(P==0)) = 1000*eps;

py = sum(P,1);
py = py / sum(py);
pya = repmat(py,size(P,1),1);

px = sum(P,2);
px = px / sum(px);
pxa = repmat(px,1,size(P,2));

pxy = P/sum(sum(P));

I = sum(sum(   pxy .* log2( pxy ./ (pxa .* pya) ) ) );

return


%testing

%p. 189-190 of Cover and Thomas
m = [0.3 0.2 0.5; 0.5 0.3 0.2; 0.2 0.5 0.3];
%following should be equal to 0:
joint2MI(m) - ( log2(3) - H([0.5 0.3 0.2])) 
