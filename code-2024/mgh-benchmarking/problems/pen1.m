% Penalty I  function
% ------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [fvec,J]=pen1(n,m,x,option)
% Dimensions -> n=variable, m=n+1
% Problem no. 23             
% Standard starting point -> x=(s(j)) where 
%                            s(j)=j 
% Minima -> f=2.24997...10^(-5)  if n=4
%           f=7.08765...10^(-5)  if n=10            
%                                     
% 11/21/94 by Madhu Lamba  
% Modified by Daniel McKenzie 14/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fvec,J] = pen1(n,m,x,option)  

J=zeros(m,n);
for i=1:n

   if (option==1 | option==3)
        fvec(i)=sqrt(1.e-5)*(x(i)-1);
    else fvec='?';
   end;        

   if (option==2 | option==3)
	J(i,i) = sqrt(1.e-5);
else J='?';
   end;
end; 
   
if (option==1 | option==3)
   sum=0;
   for j=1:n-1
       sum=sum+x(j)'*x(j);
   end;
   fvec(n)=sum-(1/4);
else fvec='?';
end;

if (option==2 | option==3)
   for j=1:n-1
        J(n+1,j) = 2*x(j);
   end;
else J='?';
end; 

if (option==1 | option==3)
   fvec=fvec';
end;

%

