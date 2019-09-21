function [dc,rho,dv,lambda,newRate,thres,lambda2star,EbN0mindB] = profgen(targetRate,max_deg_dv,dc_init);
% PROFGEN   Generates optimal LDPC code profile 
%   PROFGEN(RATE,MAX_DEG_DV,DC_INIT) generates an optimal block code profile (dv,dc) of rate defined by RATE.
%   DV, LAMBDA and RHO are optimized for given DC_INIT (concentration theorem).
%   A constraint on the maximal degree for the variable nodes can be added
%   which guarantees that MAX_DEG_DV = max(dv). The default value is 20.
%   Optimal THRESHOLD, LAMBDA_STAR(2) and EBN0_MIN (dB) are given as well.
%   MATLAB COMMAND EXAMPLE
%   >> targetRate = 1/2;max_deg_dv = 10;dc_init = [7];[dc,rho,dv,lambda,newRate,thres,lambda2star,EbN0mindB] = profgen(targetRate,max_deg_dv,dc_init);

% ---------- Initialization -----------------------

warning off MATLAB:fzero:UndeterminedSyntax; %to suppress this warning
TRUE = 1;
FALSE = 0;
f1 = inline('Phi-(1-3/x)*sqrt(pi/x)*exp(-x/4)','x','Phi');
f2 = inline('(1+1/(7*x))*sqrt(pi/x)*exp(-x/4)-Phi','x','Phi');
end_while = FALSE;
if isempty(dc_init);
    dc_init = floor(3/(1-targetRate));
else;
    dc_init = dc_init(1);
end;
dc = [dc_init,dc_init+1];       % Concentration theorem
step = 1e-1;
dv = [2:max_deg_dv].';
lambda = zeros(size(dv));
C = (1./dv).';
d = .5;
sigma_opt_dc = [0,Inf];
step = 0.1;
rho1_min = 0;
rho1_max = 1;

% ---------- Optimization---------------------

while (abs(sigma_opt_dc(1)-sigma_opt_dc(2)) > 1e-5);
    rho1 = [(rho1_min+rho1_max)/2,(rho1_min+rho1_max)/2+step];
    for lrho1 = 1:length(rho1);
        rho(1) = rho1(lrho1);rho(2) = 1-rho(1);
        sigma_min = .5;
        sigma_max = 1.17/targetRate;
        sigma = sigma_min;
        newRate = Inf;
        while ( abs(targetRate-newRate) > 1e-5 );
            s = 2/sigma^2;
            %lambda2max = max([0,(sum(rho./dc)/(1-targetRate)-1/3)/(1/2-1/3),min([(sum(rho./dc)/(1-targetRate)-1/dv(end))/(1/2-1/dv(end)),exp(1/2/sigma^2)/sum(rho.*(dc-1)),1])]);
	    lambda2max = max([0,min([exp(1/2/sigma^2)/sum(rho.*(dc-1)),1])])
            r = linspace(eps,phi(s),100);r = reshape(r,length(r),1);fr_tmp = zeros(size(r));
            for n = 1:length(r);rn = r(n);tmp = 0;for jc = 1:length(dc);tmp = tmp+rho(jc)*phi_1(1-(1-rn)^(dc(jc)-1),f1,f2);end;fr_tmp(n) = tmp;end;
            A = zeros(length(r),length(dv));for lr = 1:length(r);for iv = 1:length(dv);A(lr,iv) = phi(s+(dv(iv)-1)*fr_tmp(lr));end;end;
            A = [A;-1./(dv.')];
            b = [r;-sum(rho./dc)];
            Aeq = ones(1,length(dv));
            beq = 1;
            LB = zeros(length(dv),1);
            UB = [lambda2max;ones(length(dv)-1,1)-eps];
            lambda=lsqlin(C,d,A,b,Aeq,beq,LB,UB);
            newRate = 1-sum(rho./dc)/sum(lambda./dv);
            if newRate < targetRate;
                sigma_max = sigma;
                sigma = (sigma_min+sigma)/2;
            else;
                sigma_min = sigma;
                sigma_opt_dc(lrho1) = sigma;
                sigma = (sigma_max+sigma)/2;
            end;
        end
    end;
    if sigma_opt_dc(1)>sigma_opt_dc(2);
        rho1_max = rho1(2);
    else;
        rho1_min = rho1(1);
    end;
    step = step/2;
end;
thres = mean(sigma_opt_dc);
EbN0mindB = 10*log10(1/(2*newRate*thres^2));
lambda2star = exp(1/2/sigma^2)/sum(rho.*(dc-1));

