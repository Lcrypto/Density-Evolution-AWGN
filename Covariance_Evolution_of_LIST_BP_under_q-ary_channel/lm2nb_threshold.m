function thr=lm2nb_threshold(lambda,rho,N,flgg,delta_p)
% function thr = lm2nb_threshold(lambda,rho,N,flgg,delta_p)
%
%   This program computes the error threshold for
%  node-based LM2 decoding of irreguar LDPC codes
%  by solving the differential equation associated
%  the equivalent peeling decoder.
%
%  "lambda" is the symbol degree distribution from
%    the edge perspective starting from degree 1.
%
%  "rho" is the check degree distribution from
%    the edge perspective starting from degree 1.
%
%  "delta_p" is threshold tolerance for the bisection
%    threshold search.
%
%  Written by Fan Zhang 2/23/10
%
%(3,6) N=12000 0.2593, N=8000, 0.2593
%(4,8) N=12000 0.2398,
%(5,10) N=12000 0.2183 
%rate 0.5 [0 0.2 0.8] [0 0 0 0 0.5 0.5] 0.2623

if nargin<3
    N=12000;
end
%show figures
if nargin<4
    flgg=0;
end
if nargin<1
    lambda = [0 0.2 0.8];
    rho = [0 0 0 0 0.5 0.5];
end
    pos_n=1;
    neg_n=-1;
if nargin<5    
    delta_p=0.001;
end

%time units that the differential equation runs
total_time=N+1;
dv=length(lambda);
dc=length(rho);

%calculate the rate
sum1=0;
for i=2:dc
	sum1=sum1+rho(i)/i;
end
sum2=0;
for i=2:dv
	sum2=sum2+lambda(i)/i;
end 
rate=1-sum1/sum2
%bisection parameters
p1=0.0001;
p2=0.9999;
%channel error prob
p=0.22;

%l(k,t): fraction of correct edges connected to 
%degree k CVN at time t i 
%n(i,j,k): fraction of edges connected to CN's of type n(i,j) 
% at time i
%i,j,k start from 0
%er(t) is the fraction of remaining edges connected to IVN's
%el(t) is the fraction of remaining edges connected to CVN's
% a(t) average degree of CVN hit by CER edges at time t
% 
dec_succ=0;
while(p2-p1>delta_p)
%total system parameters    
l=zeros(total_time,dv+1);
r=zeros(dv+1,dv+1,dv+1,total_time);
n=zeros(dc+1,dc+1,total_time);
er=zeros(1,total_time);
el=zeros(1,total_time);
a=zeros(1,total_time); 
%CER system parameters
l1=zeros(total_time,dv+1);
r1=zeros(dv+1,dv+1,dv+1,total_time);
n1=zeros(dc+1,dc+1,total_time);
er1=zeros(1,total_time);
el1=zeros(1,total_time);
a1=zeros(1,total_time); 
%IER2 system parameters
l2=zeros(total_time,dv+1);
r2=zeros(dv+1,dv+1,dv+1,total_time);
n2=zeros(dc+1,dc+1,total_time);
er2=zeros(1,total_time);
el2=zeros(1,total_time);
a2=zeros(1,total_time); 
%IER1 system parameters
l3=zeros(total_time,dv+1);
r3=zeros(dv+1,dv+1,dv+1,total_time);
n3=zeros(dc+1,dc+1,total_time);
er3=zeros(1,total_time);
el3=zeros(1,total_time);
a3=zeros(1,total_time); 
%total number of edges
E=0;
ave_v_deg=0;
for i=1:dv
    ave_v_deg=ave_v_deg+lambda(i)/i;
end
E=N/ave_v_deg;
l(1,:)=E*(1-p)*[0,lambda]; 
el(1)=E*(1-p);
er(1)=E*p;
temp=0;
for i=2:dv+1
    temp=temp+l(1,i)*(i-1)/el(1);
end
a(1)=temp; 
for k=1:dc
    for i=0:k
        n(i+1,k-i+1,1)=n(i+1,k-i+1,1)+E*rho(k)*nchoosek(k,i)*(1-p)^i*p^(k-i);
    end
end
%IER1 edges
eta1=n(1,2,1);
%IER2 edges
eta2=0;
for i=2:dc
    eta2=eta2+n(i,2,1)/i;
end 
%NIE edges
eta0=er(1)-eta2-eta1;
%initial r_ijk, i NIE j IER2 k IER1
r=zeros(1+dv,1+dv,1+dv,total_time);
for i=1:dv+1
    for j=1:dv+1-i+1
        for k=1:dv+1-i-j+2
            if i==1 && j==1 && k==1
                r(i,j,k,1)=0;
            else
                r(i,j,k,1)=E*p*lambda(i+j+k-3)*nchoosek(i+j+k-3,i+j-2)*...
                nchoosek(i+j-2,i-1)*(eta0/er(1))^(i-1)*(eta2/er(1))^(j-1)...
                *(eta1/er(1))^(k-1);
            end
        end
    end
end

temp=0;
   for i=1:dv+1
       for j=1:dv+1-i+1
           for k=1:dv+1-i-j+2
               if i==1
                   continue;
               end
               temp=temp+r(i,j,k,1)/(i+j+k-3)*(i-1);
               
           end
       end
   end
%self-consistency test   
%    if abs(temp-eta0)>0.01
%        i=i;
%    end   
for t=2:total_time   
    temp=0;
    for i=1:dc+1
        for j=1:dc+1-i+1
            if i==1 && j==1 
                continue;
            end
            temp=temp+n(i,j,t-1)/(i+j-2)*(i-1);
        end
    end
%self-consistency test       
%     if isnan(el(t-1)-temp)
%         i=i;
%     end
%self-consistency test   
%     if abs(el(t-1)-temp)>0.01
%         i=i;
%     end
%self-consistency test   
%     temp=0;
%     for i=1:dc+1
%         for j=1:dc+1-i+1
%             if i==1 && j==1 
%                 continue;
%             end
%             temp=temp+n(i,j,t-1)/(i+j-2)*(j-1);
%         end
%     end 
%     abs(er(t-1)-temp);
%     if abs(er(t-1)-temp)>0.00001
%         i=i;
%     end
    %number of IER1 edges
    eta1=n(1,2,t-1);
    %number of IER2 edges
    eta2=0;
    for i=2:dc
        eta2=eta2+n(i,2,t-1)/i;
    end
%self-consistency test   
%     temp=0;
%     for j=2:dv+1
%         for i=1:dv+1-j+1
%             for k=1:dv+1-i-j+2
%                 temp=temp+r(i,j,k,t-1)/(i+j+k-3)*(j-1);
%             end
%         end
%     end 
%  
%     abs(temp-eta2);
%     if abs(temp-eta2)>0.001
%         i=i;
%     end
%     eta2=temp;
    %number of NIE edges
    eta0=er(t-1)-eta1-eta2;
%self-consistency test       
%     temp=0;
%    for j=3:dc+1
%        for i=1:dc+1-j+1
%            if i==1 && j==1
%                continue;
%            end
%            temp=temp+n(i,j,t-1)/(i+j-2)*(j-1);
%            
%        end
%    end
% abs(temp-eta0);
%    if abs(temp-eta0)>=0.00001
%        i=i;
%    end

    %number of CER nodes
    s4(t)=0;
    s4(t)=sum(n(:,1,t-1));
    
    %number of IER1 nodes
    s1(t)=0;
    for k=2:dv+1
        for i=1:dv+1-k+1
            for j=1:dv+1-i-k+2
                s1(t)=s1(t)+r(i,j,k,t-1)/(i+j+k-3);
            end
        end
    end
    %number of IER2 nodes
    s2(t)=0;
    for j=3:dv+1
        for i=1:dv+1-j+1
            for k=1:dv+1-i-j+2
                s2(t)=s2(t)+r(i,j,k,t-1)/(i+j+k-3);
            end
        end
    end
    %number of NIE nodes
    s0(t)=0;
    k=1;
    for j=1:2
        for i=1:dv+1-j+1
            if i==1 && j==1 
                continue;
            end
            s0(t)=s0(t)+r(i,j,k,t-1)/(i+j+k-3);
        end
    end
    %probability that a randomly chosen NIE node is type r_ijk
    pr_nie=zeros(dv+1,dv+1,dv+1);
    k=1;
    for j=1:2
        for i=1:dv+1-j+1
            if i==1 && j==1 
                continue;
            end            
            pr_nie(i,j,k)=r(i,j,k,t-1)/(i+j+k-3)/s0(t);
        end
    end
    %probability that a randomly chosen IER2 node is type r_ijk
    pr_ier2=zeros(dv+1,dv+1,dv+1);
    for j=3:dv+1
        for i=1:dv+1-j+1
            for k=1:dv+1-i-j+2
            if i==1 && j==1 && k==1
                continue;
            end                
                pr_ier2(i,j,k)= r(i,j,k,t-1)/(i+j+k-3)/s2(t);
            end
        end
    end
    %probability that a randomly chosen IER1 node is type r_ijk
    pr_ier1=zeros(dv+1,dv+1,dv+1);
    for k=2:dv+1
        for i=1:dv+1-k+1
            for j=1:dv+1-i-k+2
                pr_ier1(i,j,k)= r(i,j,k,t-1)/(i+j+k-3)/s1(t);
            end
        end
    end
    if s1(t)==0
        pr_ier1=zeros(dv+1,dv+1,dv+1);
    end
    
    %s1 = IER1 =c3
    %s2 = IER2 =c2
    %s4 = CER = c1
        c1=s4(t)/(s4(t)+s1(t)+s2(t));
        c2=s2(t)/(s4(t)+s1(t)+s2(t));
        c3=s1(t)/(s4(t)+s1(t)+s2(t));

%  if s4(t)==max([s1(t),s2(t),s4(t)])
%     c1=1;
%     c2=0;
%     c3=0;
% else
%     if s2(t)==max([s1(t),s2(t),s4(t)])
%         c1=0;
%         c2=1;
%         c3=0;
%     else 
%             c1=0;
%             c2=0;
%             c3=1;
%     end
% end

   %u_ipjp is the d(n_ipjp)/d(t) caused by removing one NIE edge
   u_ipjp=zeros(dc+1,dc+1);
   for ip=1:dc+1
       for jp=1:dc+1-ip+1
           if jp>=3 && jp<=dc
               u_ipjp(ip,jp)=-(jp-1)*n(ip,jp,t-1)/eta0+...
                   jp/(ip+jp-1)*(ip+jp-2)*n(ip,jp+1,t-1)/eta0;
           else
               if jp==2
                   u_ipjp(ip,jp)=jp/(ip+jp-1)*(ip+jp-2)*n(ip,jp+1,t-1)/eta0;
               else
                   if jp==dc+1
                        u_ipjp(ip,jp)=-(jp-1)*n(ip,jp,t-1)/eta0;
                   else
                        u_ipjp(ip,jp)=0;
                   end
               end
           end
       end
   end
   
                   
   %v_ipjp is the d(n_ipjp)/d(t) caused by removing one IER2 edge  
   v_ipjp=zeros(dc+1,dc+1);
   %if ip=0, v_ipjp=0
   for ip=2:dc+1
       for jp=1:dc+1-ip+1
           if jp>=3
               v_ipjp(ip,jp)=0;
           else
               if jp==2
                   v_ipjp(ip,jp)=-n(ip,2,t-1)/eta2;
               else
                   v_ipjp(ip,jp)=n(ip,2,t-1)/(ip)*(ip-1)/eta2;
               end
           end
               
       end
   end
   sum(sum(v_ipjp));
   if abs(sum(sum(v_ipjp))+1)>0.001
       i=i;
   end
       
   %w_ipjp is the d(n_ipjp)/d(t) caused by removing one IER1 edge
   w_ipjp=zeros(dc+1,dc+1);
   w_ipjp(1,2)=-1;
   
   
   %up_ijk is the d(r_ijk)/d(t) caused by the reflecting edges of one
   %removed IER2 edge
    up_ijk=zeros(dv+1,dv+1,dv+1);
    %prob that an NIE edge hit n(i,2) CN's
    pr=0;
    for i=2:dc+1-2
        pr=pr+2*n(i,3,t-1)/(i+1)/eta0;
    end

for i=1:dv+1
       for j=1:dv+1-i+1
           for k=1:dv+1-i-j+2
                     if i==1 && j==1 && k==1
                         continue;
                     end
            if i==dv+1  || j==1
                   up_ijk(i,j,k)=pr*(   (-(i-1))*r(i,j,k,t-1)/eta0+...
                     0 );
                        
            else
                   up_ijk(i,j,k)=pr*(   (-(i-1))*r(i,j,k,t-1)/eta0+...
                       i*r(i+1,j-1,k,t-1)/eta0    );
            end
           end 
       end
end

for i=1:dv+1
       for j=1:dv+1-i+1
           for k=1:dv+1-i-j+2
                     if i==1 && j==1 && k==1
                         continue;
                     end
            if i==dv+1  || k==1
                   up_ijk(i,j,k)=up_ijk(i,j,k)+ n(1,3,t-1)/eta0*...
                       (  -(i-1)*r(i,j,k,t-1)/eta0...
                      +0);
                        
            else
                  up_ijk(i,j,k)=up_ijk(i,j,k)+ n(1,3,t-1)/eta0*...
                       (  -(i-1)*r(i,j,k,t-1)/eta0...
                      + i*r(i+1,j,k-1,t-1)/eta0   ); 
            end
           end 
       end
end

%self-consistency test   
%    abs(sum(sum(sum(up_ijk)))-0);
% if abs(sum(sum(sum(up_ijk)))-0)>0.001
%     i=i;
% end

%vp_ijk is the d(r_ijk)/d(t) caused by the reflecting edges of one 
   %removed IER2 edge   
   vp_ijk=zeros(dv+1,dv+1,dv+1);
   %wp_ijk is the d(r_ijk)/d(t) caused by the reflecting edges of one
   %removed IER2 edge
   wp_ijk=zeros(dv+1,dv+1,dv+1);
   
% differential equations for CER  
    for i=1:dv+1
       for j=1:dv+1-i+1
           for k=1:dv+1-i-j+2
               if i==1 && j==1 && k==1
                   r1(i,j,k,t)=r(i,j,k,t-1);
                   continue;
               end                   
               if j==dv+1 || k==1
                   r1(i,j,k,t)=r(i,j,k,t-1)+ (a(t-1)-1)*n(2,2,t-1)/2/(el(t-1))*...
                       (-r(i,j,k,t-1)*(j-1)/eta2);
               else
                   r1(i,j,k,t)= r(i,j,k,t-1)+ (a(t-1)-1)*n(2,2,t-1)/2/(el(t-1))*...
                       (-r(i,j,k,t-1)*(j-1)/eta2+r(i,j+1,k-1,t-1)*j/eta2);
               end
           end
       end
    end
 %avoid numerical issue   
      r1(:,:,:,t)=r1(:,:,:,t)./abs(sum(sum(sum(r1(:,:,:,t))))).*abs(sum(sum(sum(r(:,:,:,t-1)))));
%self-consistency test   
%       if abs(sum(sum(sum(r1(:,:,:,t)-r(:,:,:,t-1)))))>=0.00001
%        i=i;
%    end

%diff eqn for by CER
        for i=2:dv+1
            l1(t,i)=l(t-1,i)/el(t-1)*(-i+1)+l(t-1,i);
        end
        pij1=zeros(dc+2,dc+1);
        temp=0;
        for i=1:dc+1
            for j=1:dc+1
                if(i==1 && j==1)
                    continue;
                end
                pij1(i,j)=n(i,j,t-1)*(i-1)*(a(t-1)-1)/(i+j-2)/el(t-1);
            end
        end
        pij1(dc+2,:)=zeros(1,dc+1);
        q=zeros(1,dc+2);
        temp=0;
        for i=1:dc+1
            temp=temp+n(i,1,t-1);
        end
        for i=1:dc+1
            q(i)=n(i,1,t-1)/temp;
        end
        q(dc+2)=0;
        for i=1:dc+1
            for j=1:dc+1
                if(j>1)
                    n1(i,j,t)=n(i,j,t-1)+(pij1(i+1,j)-pij1(i,j))*(i+j-2);
                end
                if(j==1)
                    n1(i,j,t)=n(i,j,t-1)+(pij1(i+1,1)-pij1(i,1))*(i-1)+(q(i+1)-q(i))*(i-1);
                end
            end
        end

        el1(t)=sum(l1(t,2:dv+1));
        er1(t)=sum(sum(sum(r1(:,:,:,t)))); 
        a1(t)=0;
        
        temp=0;
        for i=2:dv+1
            temp=temp+l1(t,i)*(i-1)/el1(t);
        end
        a1(t)=temp;
         
%% differential equations for IER2

        l2(t,:)=l(t-1,:);
        r2(:,:,:,t)=zeros(dv+1,dv+1,dv+1);
        for i=1:dv+1
           for j=1:dv+1-i+1
               for k=1:dv+1-i-j+2
                   temp=0;
                       for ip=1:dv+1
                           for jp=1:dv+1-ip+1
                               for kp=1:dv+1-ip-jp+2
                                   temp=temp+pr_ier2(ip,jp,kp)*...
                                       ((ip-1)*up_ijk(i,j,k)+(jp-1)*vp_ijk(i,j,k)...
                                       +(kp-1)*wp_ijk(i,j,k));
                               end
                           end
                       end
                   if j>=3 
                       r2(i,j,k,t)=r(i,j,k,t-1)+pr_ier2(i,j,k)*(-(i+j+k-3))...
                           +temp;
                   else
                       r2(i,j,k,t)=r(i,j,k,t-1)+ temp;
                   end
               end
           end
        end
        
        n2(:,:,t)=zeros(dc+1,dc+1);
        for ip=1:dc+1
           for jp=1:dc+1-ip+1 
               temp=0;
               for i=1:dv+1
                   for j=1:dv+1-i+1
                       for k=1:dv+1-i-j+2
                           temp=temp+pr_ier2(i,j,k)*((i-1)*u_ipjp(ip,jp)+...
                               (j-1)*v_ipjp(ip,jp)+(k-1)*w_ipjp(ip,jp));
                       end
                   end
               end
               n2(ip,jp,t)=n(ip,jp,t-1)+temp;
           end
        end
        
        el2(t)=sum(l2(t,2:dv+1));
        er2(t)=sum(sum(sum(r2(:,:,:,t))));
        a2(t)=0;
        
        temp=0;
        for i=2:dv+1
            temp=temp+l2(t,i)*(i-1)/el2(t);
        end
        a2(t)=temp;
          

% differential equations for IER2       
        l3(t,:)=l(t-1,:);
        r3(:,:,:,t)=zeros(dv+1,dv+1,dv+1);
        for i=1:dv+1
           for j=1:dv+1-i+1
               for k=1:dv+1-i-j+2
                   temp=0;
                       for ip=1:dv+1
                           for jp=1:dv+1-ip+1
                               for kp=1:dv+1-ip-jp+2
                                   temp=temp+pr_ier1(ip,jp,kp)*...
                                       ((ip-1)*up_ijk(i,j,k)+(jp-1)*vp_ijk(i,j,k)...
                                       +(kp-1)*wp_ijk(i,j,k));
                               end
                           end
                       end
                   if k>=2 
                       r3(i,j,k,t)=r(i,j,k,t-1)+pr_ier1(i,j,k)*(-(i+j+k-3))...
                           +temp;
                   else
                       r3(i,j,k,t)=r(i,j,k,t-1)+ temp;
                   end
               end
           end
        end
        
        n3(:,:,t)=zeros(dc+1,dc+1);
        for ip=1:dc+1
           for jp=1:dc+1-ip+1 
               temp=0;
               for i=1:dv+1
                   for j=1:dv+1-i+1
                       for k=1:dv+1-i-j+2
                           temp=temp+pr_ier1(i,j,k)*((i-1)*u_ipjp(ip,jp)+...
                               (j-1)*v_ipjp(ip,jp)+(k-1)*w_ipjp(ip,jp));
                       end
                   end
               end
               n3(ip,jp,t)=n(ip,jp,t-1)+temp;
           end
        end
         
        el3(t)=sum(l3(t,2:dv+1));
        er3(t)=sum(sum(sum(r3(:,:,:,t))));
        
        a3(t)=0;
        
        temp=0;
        for i=2:dv+1
            temp=temp+l3(t,i)*(i-1)/el3(t);
        end
        a3(t)=temp;
%avoid numerical issue          
    if c1==0
        l1(t,:)= l(t-1,:);
        r1(:,:,:,t)=r(:,:,:,t-1);
        el1(t)=el(t-1);
        er1(t)= er(t-1);
        a1(t)= a(t-1) ; 
        n1(:,:,t)= n(:,:,t-1);
    end
    
    if c2==0
        l2(t,:)= l(t-1,:);
        r2(:,:,:,t)=r(:,:,:,t-1);
        el2(t)=el(t-1);
        er2(t)= er(t-1);
        a2(t)= a(t-1) ; 
        n2(:,:,t)= n(:,:,t-1);
        
    end
    
    if c3==0
        l3(t,:)= l(t-1,:);
        r3(:,:,:,t)=r(:,:,:,t-1);
        el3(t)=el(t-1);
        er3(t)= er(t-1);
        a3(t)= a(t-1) ; 
        n3(:,:,t)= n(:,:,t-1);
    end
    if c1==0 && c2==0 && c3==0
        l(t,:)= l(t-1,:);
    r(:,:,:,t)=r(:,:,:,t-1);
    el(t)=el(t-1);
    er(t)= er(t-1);
    a(t)=a(t-1); 
    n(:,:,t)=n(:,:,t-1);
    else
%weighted sum        
    l(t,:)=c1*l1(t,:)+c2*l2(t,:)+c3*l3(t,:);
    r(:,:,:,t)=c1*r1(:,:,:,t)+c2*r2(:,:,:,t)+c3*r3(:,:,:,t);
    el(t)=c1*el1(t)+c2*el2(t)+c3*el3(t);
    er(t)=c1*er1(t)+c2*er2(t)+c3*er3(t);
    a(t)=c1*a1(t)+c2*a2(t)+c3*a3(t); 
    n(:,:,t)=c1*n1(:,:,t)+c2*n2(:,:,t)+c3*n3(:,:,t);
    end
%self-consistency test   
%     temp=0;
% for i=1:dc+1
%     for j=1:dc+1-i+1
%         if i==1 && j==1 
%             continue;
%         end
%         n(i,j,t);
%        temp=temp+n(i,j,t)/(i+j-2)*(i-1);
%     end
% end
% 
% 
% temp1=sum(l(t,:));
% if abs(temp-temp1)>0.01
%     
% i=i;
% end
end
    
 

if all(all(all(n(:,:,1:(total_time-1))>=neg_n))) && sum(sum(n(:,:,total_time)))<=pos_n
    dec_succ=1;
    
else
    dec_succ=0;
end 

'p'
p
if dec_succ==1
    p1=p;
    'good'
    p=(p1+p2)/2;
    n11=n;
else 
    'bad'
    p2=(p1+p2)/2;
    p=(p1+p2)/2;
end
end


for t=1:total_time
    x12(t)=n11(1,2,t);
    x13(t)=n11(1,3,t);
    x14(t)=n11(1,4,t);
    x15(t)=n11(1,5,t);
    x16(t)=n11(1,6,t);
    x17(t)=n11(1,7,t);
    x21(t)=n11(2,1,t);
    x22(t)=n11(2,2,t);
    x23(t)=n11(2,3,t);
    x24(t)=n11(2,4,t);
    x25(t)=n11(2,5,t);
    x26(t)=n11(2,6,t);
    x31(t)=n11(3,1,t);
    x32(t)=n11(3,2,t);
    x33(t)=n11(3,3,t);
    x34(t)=n11(3,4,t);
    x35(t)=n11(3,5,t);
    x41(t)=n11(4,1,t);
    x42(t)=n11(4,2,t);
    x43(t)=n11(4,3,t);
    x44(t)=n11(4,4,t);
    x51(t)=n11(5,1,t);
    x52(t)=n11(5,2,t);
    x53(t)=n11(5,3,t);
    x61(t)=n11(6,1,t);
    x62(t)=n11(6,2,t);
    x71(t)=n11(7,1,t);
    x11(t)=E-sum(sum(n11(:,:,t)));
end

if flgg==1
figure(1);
clf;
title(['IER1 curve p=',num2str(p)]);
hold
plot(x12./E,'-r');
legend('n_{0,1}');
xlabel('time');
ylabel('fractions of edges with different types');
grid


figure(2);
clf;
title(['CER curves p=',num2str(p)]);
hold
plot(x21./E,'-b');
grid
plot(x31./E,'-b');
plot(x41./E,'-b');
plot(x51./E,'-b');
plot(x61./E,'-b');
plot(x71./E,'-b');
legend('n_{1,0}','n_{2,0}','n_{3,0}','n_{4,0}','n_{5,0}','n_{6,0}');
xlabel('time');
ylabel('fractions of edges with different types');
figure(3);
clf;
title(['IER2 curves p=',num2str(p)]);
hold
plot(x22./E,'-k');
xlabel('time');
ylabel('fraction');
grid

plot(x32./E,'-k');
plot(x42./E,'-k');
plot(x52./E,'-k');
plot(x62./E,'-k');

legend('n_{1,1}','n_{2,1}','n_{3,1}','n_{4,1}','n_{5,1}');
figure(4);
clf;
title(['other curves p=',num2str(p)]);
hold
xlabel('time');
ylabel('fractions of edges with different types');
plot(x13./E,'-r');
grid
plot(x14./E,'-r');
plot(x15./E,'-r');
plot(x16./E,'-r');
plot(x17./E,'-r');

plot(x23./E,'-r');
plot(x24./E,'-r');
plot(x25./E,'-r');
plot(x26./E,'-r'); 

plot(x33./E,'-r');
 plot(x34./E,'-r');
 plot(x35./E,'-r');
 
 plot(x43./E,'-r');
 plot(x44./E,'-r');
 
 plot(x53./E,'-r');

end
i=i;
thr=p;
