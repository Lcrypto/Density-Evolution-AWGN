function [H,final_column_weights,final_row_weights]=rand_proto(n,m,lambda,rho)
%Function construct protograph or subgraph(subpart) of protograph with degree
%distribution defined by column (lambda, closes as it posible) and row
%(rho) distributions
%without eliminating balance cycles(TS) and grown of code distance pattern
%base on Implementing the Belief Propagation Algorithm in MATLAB
%R?ffer, B. S. and Kellett, C. M.
%Technical report. Department of Electrical Engineering and Computer Science, University of Newcastle, Australia, November 2008.

%Example  of parameters 
% n=32;%columns in protograph
% m=22;%row in protograph
% lambda=[0,0.1,0.2,0,0,0.7]; %Column distributions polynomial
% rho=[[0,0,0,0.1,0.05,0.15,0.2,0.1,0.1,0.2,0.1];%Row distribution polynomial
v=lambda ;
h=rho;
%Initialisation
%m=floor(n*(1-r));
H=zeros([m,n]);
alpha=[];%alphawillcontainthecolumnweightfor
%eachcolumn
for i=1:length(v)
for j=1:(floor(v(i)*n))%alwaysunderfillandthenaddextras
%later
alpha=[alpha i];
end
end
while(length(alpha)~=n)
alpha=[alpha i];
end

beta=[];%betawillcontaintherowweightforeachrow
for i=1:length(h)
for j=1:(floor(h(i)*m))%alwaysunderfillandthenaddextras
beta=[beta,i];
end
end
while(length(beta)~=m)
 beta=[beta i];
end
%Construction
for i=1:n
%construct column i
c=[];
beta_temp=beta;
for j=1:alpha(i)
%temp_row=randint(1,1,[1,m]);
% 
%  X=randint(1,Nused*Nframe,M)
%  X = randi(M, 1, Nused*Nframe) - 1;
 tic;
 temp_row=randi(m,1);
while(((beta_temp(temp_row)==0)&&...
    (max(beta_temp)>0))||...
    ((beta_temp(temp_row)<=-10)))
rng('shuffle');
time=toc;
if time>2
    error('start again several time or use different size of protograph and degree distibution');
   
end


temp_row=mod(temp_row+1,m)+1
end
c=[c temp_row];
beta_temp(temp_row)=-10;
end
%decremententriesinbeta
for k=1:length(c)
 beta(c(k))=beta(c(k))-1;
end
%populate H
for j=1:alpha(i)
H(c(j),i)=1;
end

end

%Calculateactualcolumndistribution
column_weights=H'*ones(m,1);
for i=1:max(column_weights)
count=0;
for j=1:length(column_weights)
if(column_weights(j)==i)
count=count+1;
end
end
  final_column_weights(i)=count/length(column_weights);
end
%Calculateactual rowweights
row_weights=H*ones(n,1);
for i=1:max(row_weights')
     count=0;
     for j=1:length(row_weights)
     if(row_weights(j)==i)
     count=count+1;
     end
     end
final_row_weights(i)=count/length(row_weights);
%spy(H)
end
