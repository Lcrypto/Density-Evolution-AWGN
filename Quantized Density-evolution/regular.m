function regular 
chan = zeros(1,6001);
ext = [-30 0.01 6001]; 
mapping = [-10 0.0002 50000];
dv = 3;
dc = 6;
iter = 20;
stop_pe = 1e-5;
% squre root of variance for Guass
noise=1.0; 
chan=chan_mess(ext,noise);
result_pe = de_regular(chan,iter,ext,mapping,stop_pe,dv,dc)

   