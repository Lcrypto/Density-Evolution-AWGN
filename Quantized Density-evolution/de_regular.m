function result_pe = de_regular(chan,iter,ext,mapping,stop_pe,dv,dc)
z = chan;
c = 0;
pe = 0.5;
result_pe = zeros(1,iter); 

while ((c < iter) & (pe > stop_pe))
c = c + 1;
y = new_xchk(ext, z, dc-1, mapping);
y = y / (sum(y)*ext(2));
z = new_xvar(chan, y, dv-1, ext, ext);
pe = sum(z(1:round((ext(3)-1)/2 + 1)))*ext(2);
result_pe(c) = pe;
fprintf('%f\n',pe);
end 
