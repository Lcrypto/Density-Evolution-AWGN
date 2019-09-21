function [var_people]=initial(R,IT_MAX,pe_limit)
global fixed_Rate
global var_index
global ch_index
global sigma_n
global struct

digits(100);
x=[10 200];
y(1)=log10(sqrt(pi/x(1))*exp(-x(1)/4)*(1-10.0/7.0/x(1)));
y(2)=log10(sqrt(pi/x(2))*exp(-x(2)/4)*(1-10.0/7.0/x(2)));
struct.x=x;
struct.y=y;
%IT_MAX=40; %Number of iteration in decoder
%pe_limit= 1e-10; % BER level of error
fixed_Rate=R; % Rate

first=2;last=20;number=6;% Maximal, Minimal Column Weight, Number of weight 
[var_index0]=var_nodes(first,last,number);
[var_index]=search_index(var_index0,R);
save var_index04 var_index R
load var_index04
number=length(var_index);
[var_people]=var(number);sigma_n=0.4237;






% var_no=[0.2339  0.2123  0.1468  0.1027  0.3043];
%var_index=[2 3 4 6 7 14 15 16 39 40];
% var_index=[2 3 6 7 19 20];

