%binaryDmcPerformance - Comapres three algorithms quantizing a binary-input DMC
%
%  This function plots the mutual information gap for three quantization 
%  algorithms, versus number of outputs K. Data is loaded 
%  form the file BINARYDMCPERFORMANCE.MAT.
%
%  A DMC is constructed by fine quantization of binary-input AWGN channel.
%  Quantization this DMC using the following three algorithms are compared:
%     * quantBiDmc, Dynamic programming, optimal for binary-input channels
%     * quantDmcGreedy, the greedy combining algorithm 
%     * quantDmcKLmeans, KL-means based quantization.
%
%  Refer to the function BINARYDMCPERFORMANCEDATA contains the  will geneate 
%  BINARYDMCPERFORMANCEPLOT.MAT, but this is time consuming.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

clear all
addpath('../quantDmc');

load binaryDmcPerformance

figure(1); 
clf

plot(Klist,IXZopt,'go-');
p = 6;
ha = text(Klist(p),IXZopt(p),'Optimal Quantizer','Color',[0 1 0],'FontName','Arial','HorizontalAlignment','right','VerticalAlignment','top')
hold on

plot(Klist,IXZgreedy,'rx-');
p = 2;
ha = text(Klist(p),IXZgreedy(p),'Greedy Quantizer','Color',[1 0 0],'FontName','Arial','HorizontalAlignment','left','VerticalAlignment','bottom')

plot(Klist,IXZklmeans,'bs-');
p = 3;
ha = text(Klist(p),IXZklmeans(p),'  KL Means Quantization','Color',[0 0 1],'FontName','Arial','HorizontalAlignment','left','VerticalAlignment','bottom')

xlabel('Number of Quantizer Outputs, K')
ylabel('Mutual Information Gap \Delta, log scale')
title(sprintf('Quantization of a Random Binary-Input DMC with %d Outputs',M));

yl = ylim;
ylim([1E-3 yl(2)]);
xl = xlim;
set(gca,'xtick',xl(1):1:xl(2));
set(gca,'Color','none')
set(gca,'FontName','Arial')
set(gca,'YScale','log')
grid on
