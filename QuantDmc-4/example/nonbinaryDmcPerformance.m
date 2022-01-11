%nonbinaryDmcPerformanceData - Generates data comparing two algorithms quantizing DMCs
%
%   Plots performance comparisons of QUANTDMCKLMEANS and QUANTDMCGREEDY,
%   using randomly generated DMCs with J=3 and J=6 outputs.
%
%   QuantDMC (c) Brian Kurkoski and contributors
%   Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

load nonbinaryDmcPerformance

for jj = 1:length(Jlist)
    J = Jlist(jj);
    px = ones(1,J)/J;
    for kk = 1:length(Klist)
        limit             = mean(IXY{jj,kk});
        IXZgreedy(jj,kk)  = limit - mean(IXZ2{jj,kk});
        IXZklmeans(jj,kk) = limit - mean(IXZ3{jj,kk});
    end
end

figure(3);
clf;
hold on
for jj = 1:length(Jlist)
    plot(Klist,IXZgreedy(jj,:),'rx-');
    if jj == 1
        ha = addlabel(sprintf('Greedy Quantizer, J=%d',Jlist(jj)),1,38);
    else
        ha = addlabel(sprintf('J=%d',Jlist(jj)),1,33);
    end
    ha.FontName = 'Arial';
    
    plot(Klist,IXZklmeans(jj,:),'bs-');
    if jj == 1;
        ha = addlabel(sprintf('KL Means Quantization, J=%d',Jlist(jj)),3,30);
    else
        ha = addlabel(sprintf('J=%d',Jlist(jj)),3,33);
    end
    ha.FontName = 'Arial';
end
nicecolor

xlabel('Number of Quantizer Outputs, K')
ylabel('Mutual Information Gap')
ylabel('Mutual Information Gap \Delta, log scale')

title(sprintf('quantization of a random DMC with %d outputs',M));

yl = ylim;
ylim([1E-3 3E-1])
set(gca,'xtick',Klist);
xlim([Klist(1)-2 Klist(end)+2]);

set(gca,'Color','none')
set(gca,'FontName','Arial')
set(gca,'YScale','log')

