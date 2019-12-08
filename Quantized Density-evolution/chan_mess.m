%channel output signal 'distribution is N(2/deta^2,4/deta^2)
% apply error function to evalute the value of distribution funtion 
% then the probability in a unit area can be found 
function chan_ini=chan_mess(ext,deta)
mean=2/deta^2;
var_sqrt=2/deta;
chan_ini=zeros(1,ext(3));
step=[ext(2) ext(2) (ext(3)-1)/2];
for i=1:step(3)-1
    chan_ini(i+step(3)+1)=(erf(((i+0.5)*step(2)-mean)/(sqrt(2)*var_sqrt))-erf(((i-0.5)*step(2)-mean)/(sqrt(2)*var_sqrt)))/2;
    chan_ini(i+1)=(erf(((i-step(3)+0.5)*step(2)-mean)/(sqrt(2)*var_sqrt))-erf(((i-step(3)-0.5)*step(2)-mean)/(sqrt(2)*var_sqrt)))/2;
end
    chan_ini(1)=(erf(((-step(3)+0.5)*step(2)-mean)/(sqrt(2)*var_sqrt))+1)/2;
    chan_ini(step(3)+1)=(erf(((+0.5)*step(2)-mean)/(sqrt(2)*var_sqrt))-erf(((-0.5)*step(2)-mean)/(sqrt(2)*var_sqrt)))/2;
    chan_ini(ext(3))=1-(erf(((step(3)-0.5)*step(2)-mean)/(sqrt(2)*var_sqrt))+1)/2;
    chan_ini=chan_ini/ext(2); 