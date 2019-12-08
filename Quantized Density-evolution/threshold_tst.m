% given degree distibution find the threshold
function  ns_tmp= threshold_tst(vard,chkd,ns_str_low,ns_str_high,ns_stp,ext,mapping,iter,stop_pe)
% vard variable node degree
% chkd check nodedegree
% ns_str start noise
% ns_stp decrease noise step 
% chan channle initial message
% ext variable and channel form of the message 
% mapping check node messange domain
% iter iterative number
% stop_pe if error probability smaller than, stop without doing anything
flag=1;
% ns_tmp=ns_str;
[row,dv]=size(vard);
[row,dc]=size(chkd);  

while (ns_str_high-ns_str_low>=0.01)
  ns_tmp=(ns_str_high+ns_str_low)/2;
   chan=chan_mess(ext,ns_tmp);
z = chan;
c = 0;
pe = 0.5;
pe_former=0;
while (c < iter) 
c = c + 1;
y_ave=0;   
for i=2:dc
    if chkd(i)~=0
        y = new_xchk(ext, z, i-1, mapping);
        y = y / (sum(y)*ext(2));
        y_ave=y_ave+chkd(i)*y;
    end
end
xvar_ave=0;
for j=2:dv
    if vard(j)~=0
        z = new_xvar(chan, y_ave, j-1, ext, ext);
        xvar_ave=xvar_ave+vard(j)*z;
    end
end 
z=xvar_ave;
pe = sum(z(1:round((ext(3)-1)/2 + 1)))*ext(2)

if(abs(pe_former-pe)<1e-5)
    break;
end

pe_former=pe;
end
if pe>stop_pe
    ns_str_high=ns_tmp;
else
    ns_str_low=ns_tmp;
end
end


