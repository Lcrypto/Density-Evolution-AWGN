% biAwgnQuantization - Quantizes a DMC created from a binary-input AWGN channel
%
% Makes a plot of how a binary-input AWGN channel is quantized to K=8 levels,
% when the AWGN channel has first been finely quantized to M=30 levels.
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

function biAwgnQuantization

addpath('../quantDmc');

M = 30;
K = 8;
var = [0.1 0.35 0.6];

clf

for ii = 1:length(var)
    %Create a fine DMC with M outputs from an AWGN channel
    [P,FineQuantizer,Boundary] = biAwgn2Dmc(var(ii),M);
    
    %perform quantization
    Q = quantBiDmcMulti(P,K);
    Q = Q{1};
    
    
    %Make the plot
    %
    %The narrow vertical lines show the fine quantizer boundaries
    %
    %The colors indicate fine channel outputs that quantize to the same
    %quantizer channel output
    
    subplot(3,1,ii);
    hold on
    
    t = linspace(-3,3,1001);
    Delta = t(2) - t(1);
    
    Boundary(1)   = t(1);   %replace +/- Inf with real numbers
    Boundary(end) = t(end);
    for jj = 1:length(Q)
        x0 = Boundary(jj);
        x1 = Boundary(jj+1);
        
        s = linspace(x0,x1,round((x1-x0)/Delta));
        patchY = [gaussmax(s,var(ii)) 0 0];
        patchX = [s x1 x0];
        
        colors = {'r' , 'g', 'm', 'b'};
        quantizerNumber = find(Q(:,jj) == 1);
        useColor = mod(quantizerNumber,length(colors)) + 1;
        patch(patchX,patchY,colors{useColor})
    end
    
    plot(t,gauss(t,-1,var(ii)),'k-','linewidth',2);
    hold on
    plot(t,gauss(t, 1,var(ii)),'k-','linewidth',2);
    ylim([0 1.4]);
    grid on
    title(sprintf('AWGN var = %g, M=%d DMC outputs, K=%d quantizer outputs',var(ii),M,K))
    ylabel('Prob. distribution on AWGN output')
end
xlabel('AWGN channel output')


function g=gauss(t,mean,var)
k = 1 / (sqrt(2*pi*var) );
g = k * exp( - (t - mean) .^ 2 / (2 * var) ) ;


function g=gaussmax(t,var)

k = 1 / (sqrt(2*pi*var) );
mean = 1;
g1 = k * exp( - (t - mean) .^ 2 / (2 * var) ) ;
mean = -1;
g0 = k * exp( - (t - mean) .^ 2 / (2 * var) ) ;
g = max([g0 ; g1]);



