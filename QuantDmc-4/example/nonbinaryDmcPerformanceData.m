%nonbinaryDmcPerformanceData - Generates data comparing three algorithms quantizing a binary-input DMC
%
% See NONBINARYDMCPERFORMANCEPLOT for more details
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

Jlist = [3 6];
M = 512;

Klist = 6:3:54; %suitable for M=128
Klist = 6:4:70; %suitable for M=512
iter  = 50;

IXY  = cell( length(Jlist), length(Klist) );
IXZ2 = cell( length(Jlist), length(Klist) );
IXZ3 = cell( length(Jlist), length(Klist) );

for ii = 1:iter
    tic
    for jj = 1:length(Jlist)
        J = Jlist(jj)
        px = ones(1,J)/J;
        for kk = 1:length(Klist)
            
            K = Klist(kk);
            pygx = rand(J,M);
            for j = 1:J
                pygx(j,:) = pygx(j,:) / sum(pygx(j,:));
            end
            
            [Q,IXZ2{jj,kk}(end+1)] = quantDmcGreedy(pygx,K);
            [Q,IXZ3{jj,kk}(end+1)] = quantDmcKLmeans(pygx,K);
            pxy                    = jointDistribution(pygx,px);
            IXY{jj,kk}(end+1)      = joint2MI(pxy);
        end
        
    end
    toc
end
save nonbinaryDmcPerformance


