function [s,p]=sparsity(g,theta,epsilon)

% This code estimates the numerical (effective) sparsity of vector g

[n m]=size(g);
s1=0;
for i=1:n
    if abs(g(i))>theta*epsilon
        s1=s1+1;
    end
end
g1=abs(g);
v=sort(g1, 'descend');
b=zeros(n,1);
for i=1:n
    b(i)=log(v(i));
end
A=ones(n,2);
for i=1:n
    A(i,1)=log(i);
end
z=A\b
p=1;
if z(1)~=0
    p=-1/z(1);
end
s2=n;
if and(p>0,p<1)
    a=sqrt((2/p)-1);
    b=((1/p)-1);
    c=(2*p)/(2-p);
    s2=((4/theta)*((13.2/a)+(11/b)))^(c);
end
s=min([s1 s2]);

