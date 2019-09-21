function   [sigma_n]= price_evalution(var,iter_GA) 
global var_degree 
global var_index 
global fixed_Rate 
global ch_index 
global ch_degree 
global sigma_n 
global count 
 
tt=1-sum(var);var=[tt,var]; 
var_degree(var_index)=var; 
[ch_degree,ch_index]=rate_redress2(var_degree(var_index),fixed_Rate,var_index); 
sigma_PRECISION=0.1; 
[sigma_n]=ir_GA(sigma_n,sigma_PRECISION,iter_GA); 
sigma_n 
 
 
 
 
 
 
 
 
count=count+1; 
% Eb_N0_dB=0.4822;     %SNR per bit, Eb/N0=1/R/N0=1/R/(2var)=1/(2R*var)-->1/var=2R*(Eb/N0) 
% var_ch=8*R*10^(Eb_N0_dB/10);    %LLR��� var_ch=4/var=4*2R*(Eb/N0)=8R*Eb/N0 
% sigma_n=sqrt(4/var_ch);         %sigma_n=sqrt(var) 
% sigma_n=0.9460;      %Eb_N0_dB=10*log10(1/(sigma_n*sigma_n*2*R)) 
% 1-sum(ch_degree(ch_index)./ch_index)/sum(var_degree(var_index)./var_index);