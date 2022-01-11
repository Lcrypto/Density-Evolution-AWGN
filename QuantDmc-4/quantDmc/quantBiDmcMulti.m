function [Qout,mi] = quantBiDmcMulti(Pin,K)
% [Q,MI] = quantBiDmcMulti(P,K)
%
% Finds the optimal quantizer of DMC given by P, quantized to K values,
% in the sense of maximizing mutual information.
% P is a 2-by-M matrix, where:
%   P(j,m) = Pr( Y=m | X=j ),
% for a DMC with inputs X and outputs Y.  K is an integer, generally less
% than M. 
%
% Q is a cell array containing M-by-K matrices, each one is an optimal quantizer.
% For each matrix Q{i}, Q{i}(m,k) is a 1 if DMC output m is quantized
% to k, and otherwise is a 0.
%
% MI is the mutual information between the channel input and the quantizer
% output, which is the same for all optimal quantizers.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE
%

J = size(Pin,1);
I = size(Pin,2);
if J ~= 2
    error('This only works with binary matrices');
end

if K >= I
    T = Pin;
    mi = NaN;
    Qout = {eye(size(Pin,2))};
    return
end

%Input Pin is conditional, but this function assumes Pin is joint.
Pin = Pin / 2;

%Pin is joint, construct Pcond to sort.
Pcond = Pin ./ repmat( sum(Pin,2),1,I);
LLR = log( Pcond(1,:) ./ Pcond(2,:) );
[t,sortorder] = sort(LLR);

%Sorted P is used from now on
P = Pin(:,sortorder);



%initial distance computation
dist = zeros(I,I);
pj = sum(P,2);
for ii = 1:I
    for kk = ii:I
        t = sum(P(:,ii:kk),2);
        s = sum(t);
        %if s > 1; s =1; end
        dist(ii,kk) = sum( (t .* log2(t ./ (s * pj) ) ) +eps );
    end
end

SM = zeros(I,K);
ps = cell(I,K);
SM(:,1) = dist(1,:);
fl = 0;
for kk = 2:K
    for ii = kk:I
        t = zeros(size([kk-1:ii-1]) );
        for ell = kk-1:ii-1
            t(ell - (kk-2) ) =  SM(ell,kk-1) + dist(ell+1,ii) ;
        end
        [SM(ii,kk),ps{ii,kk}] =  max(t);
        %ps(ii,kk) = ps(ii,kk) + kk - 2;
        ps{ii,kk} =  find(t == max(t));
        ps{ii,kk} = ps{ii,kk} + kk - 2;
        %if length(ps{ii,kk}) > 1
        %    fl = 1;
        %end
    end
end

%build quantizer list
Q = [I];
for kk = K:-1:1
    Qnew = [];
    for ii = 1:size(Q,1)
        s = Q(ii,K-kk+1);
        t = ps{s,kk};
        if length(t) == 0;
            t = 0;
        end
        if length(t) == 1
            Q(ii,K-kk+2) = t;
        else
            Q(ii,K-kk+2) = t(1);
            Qt = repmat(Q(ii,1:K-kk+1),length(t)-1,1);
            tp = t(:);
            Qt = [Qt tp(2:end)];
            Qnew = [Qnew; Qt];
        end
    end
    Q = [Q; Qnew];
end

%build the quantizer from the quantizer list
Qlist = Q;
Q={};
for ii = 1:size(Qlist,1)
    Q{ii} = zeros(K,I);
    for kk = 1:K
        Q{ii}(K-kk+1, Qlist(ii,kk+1)+1:Qlist(ii,kk)) = 1;
    end
end

%compute mutual information associated with each quantier.
for ii = 1:length(Q)
    T = P * Q{ii}' ;
    p=sum(T,2);
    q=sum(T,1);
    mi(ii)=0;
    for kk = 1:K
        for jj = 1:J
            mi(ii) = mi(ii) + T(jj,kk) * log2( T(jj,kk) / (q(kk) * p(jj))) ;
        end;
    end
end
%disp([ SM(I,K) mi])

%reverse the sort order to agree with Pin input
for ii = 1:length(Q);
    Q{ii}(:,sortorder) = Q{ii};
end    

Qout = Q;

return 
