%binaryDmcPerformanceData - Generates data comparing three algorithms quantizing a binary-input DMC
%
%  Generates the data file BINARYDMCPERFORMANCE.MAT.
% 
%  See BINARYDMCPERFORMANCE for details.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

M = 128;
px = [0.5 0.5];

Klist = 5:10;
iter = 100;

Klist = 2:20;
iter = 1000;

for kk = 1:length(Klist)
    K        = Klist(kk);
    
    pygx     = biAwgn2Dmc(1,M);
    pxy      = jointDistribution(pygx,px);
    IXY      = joint2MI(pxy);
    [Q,IXZ1] = quantBiDmc(pygx,K);
    [Q,IXZ2] = quantDmcGreedy(pygx,K);
    
    IXZ3 = zeros(1,iter);
    for ii = 1:iter
        [Q,IXZ3(ii)] = quantDmcKLmeans(pygx,K);
    end
    
    IXZopt(kk)     = IXY - IXZ1;
    IXZgreedy(kk)  = IXY - IXZ2;
    IXZklmeans(kk) = IXY - mean(IXZ3);
end

save binaryDmcPerformance