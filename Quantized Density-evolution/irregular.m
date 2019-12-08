function irregular 
chan = zeros(1,6001);
ext = [-30 0.01 6001]; 
mapping = [-10 0.0002 50000];
dv = 6;
dc = 6; 
iter = 100;
stop_pe = 1e-5; 
% vard(1,:)=[0.000000 0.201350 0.102124 0.024210 0.432443 0.139250 0.015776 0.084848];
% chkd(1,:)=[0.000000 0.000000 0.000000 0.000000 0.000000 0.244334 0.061851 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.693815];
% vard(1,:)=[0 0.062500 0.495536 0.017857 0.424107 ];
% chkd(1,:)=[0 0 0 0 0 0 1.000000 0.000000 0.000000 ];
% vard(1,:)=[0 0.261905 0.000000 0.380952 0.000000 0.357143];
% chkd(1,:)=[0 0 0 0 0 1.000000];
vard(1,:)=[0 0.2895    0.3158 0 0  0.3947];
chkd(1,:)=[ 0 0 0 0 0 0.9032  0.0968];
% vard=zeros(1,dv);
% chkd=zeros(1,dc);
% vard(2)=0.38354;
% vard(3)=0.04237;
% vard(4)=0.57409; 
% chkd(5)=0.24123; 
% chkd(6)=0.75877;
% vard(2)=0.2524580000; 
% vard(3)=0.0748275000;
% vard(5)=0.6727150000; 
% chkd(14)=1;
% squre root of variance for Guass  
% vard(3)=1;
% chkd(6)=1;
noise=0.872345;%from GA
ns_str_low=0.9;   
ns_str_high=0.94;     
ns_stp=0.01; 
% chan=chan_mess(ext,noise); 
% result_pe = de_irregular(chan,iter,ext,mapping,stop_pe,vard,chkd,dv,dc);
threshold=threshold_tst(vard,chkd,ns_str_low,ns_str_high,ns_stp,ext,mapping,iter,stop_pe);
fprintf('%8f', threshold);

   