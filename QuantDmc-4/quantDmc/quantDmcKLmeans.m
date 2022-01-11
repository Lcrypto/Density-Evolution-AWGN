function [Q,IXZ,means,clusters,clusterHistory] = QuantMiDmcKLmeans(pygx,K,px)
%klmeans  KL means clustering quantization of a discrete memoryless channel.
%
%   Q is an M-by-K quantization matrix such that the conditional distribution 
%   on the quantizer output Z Pr(Z=z|X=x) = Pr(Y=x|X=x) * Q(y,z) is given by
%      PZGX = PYGX * Q
%
%   IXZ is the mutual information I(X;Z).
%
%   If empty clusters are found, then columns are removed from Q, and 
%   Q may have fewer than K columns.
%
% Alan Zhang and Brian Kurkoski
% Distributed under an MIT-like license; see the file LICENSE
%

if nargin < 2
    error('two input arguments are required');
end

maxItr = 100; % maximum number of iteration 

if nargout > 3
    trackHistory = true;
else
    trackHistory = false;
end

% Calculate Pr(X|Y)
J    = size(pygx,1); % number of inputs
M    = size(pygx,2); % number of outputs
if nargin < 3
    px   = ones(1,J) / J;    %default is uniform input distribution
end
pxy  = repmat(px(:),1,M) .* pygx;
py   = sum(pxy,1);
pxgy = pxy' ./ repmat(py(:)',J,1)';

% Random mean initialization
% sample = randperm(M);
% means  = zeros(K,J);
% for i = 1:K
%     means(i,:) = pxgy(sample(i),:);
% end

% Seeding mean initialization
means = seed(pxgy,K);

thisMeans = zeros(size(means));

if trackHistory
    clusterHistory = cell(1,maxItr);
end

% Assign-update iteration until no change occurs
for i = 1:maxItr
    
    % assign %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z = means;
    y = pxgy;
    
    clusters = zeros(1,M);
    for m = 1:M
        dist = zeros(K,1); % distance between y(i) and the K means
        p = y(m,:);
        for k = 1:K
            q = z(k,:);
            dist(k) = sum(p .* log2(p./q) );
        end
        [~, cIdx] = min(dist);
        clusters(m) = cIdx;
    end
    
    if trackHistory
        clusterHistory{i} = clusters;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:K
        idx = find(clusters == jj);
        thisMeans(jj,:) = mean(pxgy(idx,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    newMeans = thisMeans;
    
    stop = isequal(newMeans,means); % stop if no change occurs
    if stop
        fprintf('stable after %d iterations\n',i)
        break;
    end
    means = newMeans; % if means change then repeat
end

%Output quantizer Pr[Z|Y]
Q = zeros(M,K);
for ii = 1:K
    idx = find(clusters == ii);
    Q(idx,ii) = 1;
end

%remove empty columns from Q
emptyClusters = find(sum(Q) == 0);
if ~isempty(emptyClusters)
    nonEmpty = setdiff(1:K,emptyClusters);
    Q = Q(:,nonEmpty);
    fprintf('found %d empty clusters\n',length(emptyClusters));
    K = size(Q,2);
end

%compute mutual information    
pzgx = pygx*Q;
pxz  = repmat(px(:),1,K) .* pzgx;
pz   = sum(pxz,1);
IXZ  = sum(sum( pxz .* log2(pxz ./ (px(:) * pz) )));

if trackHistory
    clusterHistory = clusterHistory(1:i);
end
