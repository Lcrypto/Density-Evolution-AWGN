function [m_last]=decompose(vec,deta)

LL=length(vec);
x=vec(1);y=vec(2);factor=1.1;
j=1;
matrix=[];tmp=vec;
while x>factor*deta
x=x-deta;
y=y+deta;
tmp(j:j+1)=[x y];
matrix=[matrix;tmp];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:LL-2
[m,n]=size(matrix);
for ii=1:m
  x=matrix(ii,j);
  y=matrix(ii,j+1);tmp=matrix(ii,:);
  while x>factor*deta
   x=x-deta;
   y=y+deta;
   tmp(j:j+1)=[x y];
   matrix=[matrix;tmp];
  end
end
end


j=LL-1;m_last=[];
[m,n]=size(matrix);
for ii=1:m
  x=matrix(ii,j);
  y=matrix(ii,j+1);tmp=matrix(ii,:);
  while x>factor*deta
   x=x-deta;
   y=y+deta;
   tmp(j:j+1)=[x y];
   matrix=[matrix;tmp];m_last=[m_last;tmp];
  end
end
