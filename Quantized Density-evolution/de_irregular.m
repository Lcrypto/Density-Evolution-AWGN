function result_pe = de_irregular(chan,iter,ext,mapping,stop_pe,vard,chkd,dv,dc)
z = chan;
c = 0;
pe = 0.5;
result_pe = zeros(1,iter); 
round(0.3);
while ((c < iter) & (pe > stop_pe))
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
pe = sum(z(1:round((ext(3)-1)/2 + 1)))*ext(2);
result_pe(c) = pe;
fprintf('%f\n',pe); 
end 
