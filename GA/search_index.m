function [var_index_out]=search_index(var_index0,fixed_Rate)
global var_degree
global var_index
global ch_degree
global ch_index
[M,N]=size(var_index0);sigma_n=0.2;
degree=rand(1,N);var_index=var_index0(1,:);var_degree=[];var_degree(var_index)=degree/sum(degree);
[ch_degree,ch_index,dc_av]=rate_redress2(var_degree(var_index),fixed_Rate,var_index);
sigma_PRECISION=0.1;
[sigma_n]=ir_GA(sigma_n,sigma_PRECISION,4);
for loop=1:10
for i=2:M
  degree=rand(1,N);
  var_degree=[];var_degree(var_index)=degree/sum(degree);var_index=var_index0(i,:);
  [ch_degree,ch_index,dc_av]=rate_redress2(var_degree(var_index),fixed_Rate,var_index);
  sigma_n2=sigma_n+0.0001;
  iter_mu0=ir_GA2(sigma_n2);
  if iter_mu0==1
  sigma_PRECISION=0.1;
  [sigma_n]=ir_GA(sigma_n,sigma_PRECISION,5)
  var_index_out=var_index;
  end
end
for i=2:M
  degree=abs(randn(1,N));
  var_degree=[];var_degree(var_index)=degree/sum(degree);var_index=var_index0(i,:);
  [ch_degree,ch_index,dc_av]=rate_redress2(var_degree(var_index),fixed_Rate,var_index);
  sigma_n2=sigma_n+0.0001;
  iter_mu0=ir_GA2(sigma_n2);
  if iter_mu0==1
  sigma_PRECISION=0.1;
  [sigma_n]=ir_GA(sigma_n,sigma_PRECISION,5)
  var_index_out=var_index;
  end
end
end

