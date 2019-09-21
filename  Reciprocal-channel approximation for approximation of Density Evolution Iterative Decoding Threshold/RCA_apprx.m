function [flag] = RCA_apprx(HB,SNR,p,R,snr_R,infor,iter)   
  
s = 2*10^(SNR/10);

iter_num = iter;
mvout_threshold = 30;

[m,n] = size(HB);
HBTmp = HB;

mu0 = ones(1,n) * s;
mu0 (1:p) = zeros(1,p);
mu = zeros(m,n);
mv = zeros(m,n);

for j = 1:n
    mv(:,j) = mu0(j) * HBTmp(:,j);
end
mvout = mu0;
muout = zeros(m,1);
flag = 0;
for l = 1:iter_num
    
    for i = 1:m
        idx = find(HB(i,:) ~= 0);
        mvidx = mv(i,idx);
        mvtmp = interp1(snr_R,R,mvidx);
        muout(i) = sum(HB(i,idx).*mvtmp);
        for j = 1:length(idx)
            mu(i,idx(j)) = muout(i) -  mvtmp(j);
        end
    end
    mu(mu>50) = 50;
    
    for j = 1:n
        idx = find(HB(:,j) ~= 0);
        muidx = mu(idx,j);
        mutmp = interp1(snr_R,R,muidx);
        mvout(j) = mu0(j) + sum(HB(idx,j).*mutmp);
        for i = 1:length(idx)
            mv(idx(i),j) = mvout(j) - mutmp(i);
        end
    end
    mv(mv>50) = 50;
    
    if(min(mvout(1:infor))>mvout_threshold)
        flag = 1;
        break;
    end
end
end
