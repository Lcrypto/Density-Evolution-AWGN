function [var_profile]=var(number) 
L=40*(number+1); 
var_profile=zeros(L,number-1); 
for i=1:L 
    var=rand(1,number);var=var/sum(var); 
    var_profile(i,:)=var(2:end); 
    if mod(i,2) 
    var=randn(1,number);var=abs(var);var=var/sum(var); 
    var_profile(i,:)=var(2:end); 
    end 
end