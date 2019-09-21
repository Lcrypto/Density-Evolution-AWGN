clear all
clc
global var_degree
global ch_degree
global count
global fixed_Rate
global var_index
global ch_index
global IT_MAX
global pe_limit
R=0.5; % rate of code
IT_MAX=40; %Number of iteration in BP decoder
pe_limit= 1e-10; % BER level of error
count=0;
[var_people]=initial(R,IT_MAX,pe_limit);

[I_NP,I_D]=size(var_people);
I_itermax = 50;iter_GA=4; % Setting for differential evolution iteration of DE
% F_weight    DE-stepsize F_weight ex [0, 2]
F_weight = 0.80;
% F_CR      crossover probabililty constant ex [0, 1]
F_CR = 0.8;

% Type of  differential evolution (DE) optimization strategy 
% https://en.wikipedia.org/wiki/Differential_evolution
%1 DE/rand/1 ; 
%2 DE/local-to-best/1
%3 DE/best/1 with jitter
%4 DE/rand/1 with per-vector-dither
%for detail 94 line in deopt.m     
I_strategy = 3; 
S_struct.I_NP    = I_NP;
S_struct.F_weight  = F_weight;
S_struct.F_CR    = F_CR;
S_struct.I_D     = I_D;
S_struct.I_itermax  = I_itermax;
S_struct.I_strategy = I_strategy;
S_struct.iter_GA = iter_GA;

[FVr_x,sigma,I_nf,FM_pop] = deopt(S_struct,var_people)



