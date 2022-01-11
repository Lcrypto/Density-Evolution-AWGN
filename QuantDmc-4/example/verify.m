% VERIFY - Compares the output of QUANTBIDMC and QUANTBIDMCMULTI
%
% QuantDMC contains two implementations of the quantizer for binary-input
% DMCs: QUANTBIDMC is MEX-based and faster, but outputs only one optimal 
% quantizer.  QUANTBIDMCMULTI is Matlab-based and slower, but outputs all
% optimal quantizers.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

ntest = 100;

Mlist = [3 4 5 6 7];
Mo = min(Mlist)-1;
Mo = 0;

%identicalQuantizer = zeros(max(Mlist));
maxMiError = -Inf * ones(max(Mlist)-Mo,max(Mlist)-1);
clf

for M = Mlist
    lg = {'none'};
    for K = 2:M-1
        lg{end+1} = ['K=', num2str(K)];
        for ii = 1:ntest
            
            P = randomDmc(2,M);
            
            [Q1,mi1] = quantBiDmcMulti(P,K);
            Q1 = Q1{1};  %just select the first quantizer
            mi1 = mi1(1); 
            [Q2,mi2] = quantBiDmc(P,K);
            
            if any(Q1(:) ~= Q2(:))
                error('Quantizers do not match');
            end
            
            maxMiError(M-Mo,K) = max( abs(mi1-mi2) , maxMiError(M-Mo,K) );
        end
        
        plot(maxMiError,'o-')
        xlabel('DMC outputs M');
        ylabel('maximum error in mutual information')
        legend(lg,2)
        xlim([min(Mlist) , max(Mlist)])
        set(gca,'xtick',Mlist);
        drawnow;
        
    end
end

fprintf('All quantizers are identical (%d tests performed).\n',ntest)

maxMiError
