function [x]=inv_Phi_2(y,start)
start=10;
  x=start;inv_delta= 10;inv_limit= 10000;
  while x<inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref<y) break; end 
x=x+inv_delta;
    end
    
   while x>inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref>y) break;  end
x=x+inv_delta/10;
    end
    
    while x<inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref<y) break;  end 
x=x+inv_delta/100;
    end
    
    while x>inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref>y) break; end 
x=x+inv_delta/1000;
    end
    
    while x<inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref<y) break;  end 
x=x+inv_delta/10000; 
    end
    
    while x>inv_limit
y_ref=sqrt(pi/x)*exp(-x/4)*(1-10.0/7.0/x);
if(y_ref>y) break; end 
x=x+inv_delta/100000;
    end