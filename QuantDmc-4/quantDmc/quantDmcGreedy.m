function [Q,IXZ] = QuantMiDmcGreedy(pygx,K,px)
%greedyCombining Greedy combining quantization of a discrete memoryless channel.
%   [Q,IXZ] = GREEDYCOMBINING(PYGX,PX,K) Quantizes a discrete memoryless channel
%   (DMC) PYGX with input distribution PX to K levels.  PYGX has J rows and M
%   columns where the columns sum to one.  If K is greater than the number
%   of columns M, then no quantization is performed.
%
%   Q is an M-by-K quantization matrix such that the conditional distribution 
%   on the quantizer output Z Pr(Z=z|X=x) = Pr(Y=x|X=x) * Q(y,z) is given by
%      PZGX = PYGX * Q
%
%   IXZ is the mutual information I(X;Z).
%
% Brian Kurkoski
% Distributed under an MIT-like license; see the file LICENSE
%

if nargin < 2
    error('two input arguments are required');
end

[J M] = size(pygx);

if nargin < 3
    px = ones(1,J)/J;
end

pxy = repmat(px(:),1,M) .* pygx;

py = sum(pxy,1);
px = sum(pxy,2);
 
iself = zeros(1,M);
for yy = 1:M
    iself(yy) = sum( pxy(:,yy) .* log2( pxy(:,yy) ./ (px(:)*py(yy)) ) );
end



%Initial distance computation
dist = ones(M,M)*Inf;
for ii = 1:M-1
    for jj = ii+1:M
        t = sum((pxy(:,ii) + pxy(:,jj)) .* log2( (pxy(:,ii)+pxy(:,jj)) ./ (px(:)*(py(ii) + py(jj))) ) );
        dist(ii,jj) = iself(ii) + iself(jj) - t;
        %dist(ii,jj) =  t;
    end
end

Q = eye(M);

%BEGIN LOOP
while(M > K)
    [imin,jmin]  = find(dist == min(min(dist)),1);
    
    %combining (jmin will be deleted)
    pxy(:,imin) = pxy(:,imin) + pxy(:,jmin);
    
    %update distances
    py(imin) = sum(pxy(:,imin));
    iself(imin) = sum( pxy(:,imin) .* log2( pxy(:,imin) ./ (px(:)*py(imin)) ) );
    
    %update distances in row imin
    for j = imin+1:length(dist)
        t            = sum( (pxy(:,imin) + pxy(:,j)) .* log2(  (pxy(:,imin)+pxy(:,j)) ./ (px(:)*(py(imin) + py(j))) ) );
        dist(imin,j) = iself(imin) + iself(j) - t;
    end
    
    %update distance in column imin
    for i = 1:imin-1
        t            = sum( (pxy(:,i) + pxy(:,imin)) .* log2(  (pxy(:,i)+pxy(:,imin)) ./ (px(:)*(py(i) + py(imin))) ) );
        dist(i,imin) = iself(i) + iself(imin) - t;
    end
    
    %build the quantizer
    Qthis = eye(M-1);
    c     = Qthis(:,imin);
    Qthis = [ Qthis(:,1:jmin-1) , c , Qthis(:,jmin:end) ];
    Q     = Qthis * Q;

    %delete the merged term
    r = [1:jmin-1 jmin+1:length(dist)];
    pxy = pxy(:,r);
    py = py(r);
    iself = iself(r);
    dist = dist(r,:);
    dist = dist(:,r);
    
    [J,M] = size(pxy);
end

Q = Q';

pzgx = pygx*Q;
pxz  = repmat(px(:),1,M) .* pzgx;
pz   = sum(pxz,1);
IXZ  = sum(sum( pxz .* log2(pxz ./ (px * pz) )));

return
