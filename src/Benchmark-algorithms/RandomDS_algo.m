%
%  RDS (STP) Algorithm
function [                                                              ...
    x_opt  ,                                                            ...
    f_opt  ,                                                            ...
    fin    ,                                                            ...
    k,                                                                   ...
    vecfun,                                                               ...
    fcount                                                               ...
    ]                                                                   ...
    = RandomDS_algo(une_f,                                                ...
    f0,                                                                  ...
    un_x0          ,                                                    ...
    fm,                                                                 ...
    un_nit_max,max_feval,epsilon,version,muP,L,lamda)
%**********************
%    CODE             *
%**********************
t0 = cputime;
fcount    = 1                                                            ;                                                                ;
% time will be a colone vector    %
x          = un_x0                                                        ;
n          = length(x)                                                    ;
II  =eye(n);
k          =   1                                                          ;
fin        =   0                                                          ;
f =  f0; % matters when noisy case
% feval(une_f,x); % fcount must be 1 beacuse of this
% f0 = f;
vecfun = [f0];
%vecfun = [vecfun,f];
%f_count    = f_count    +1;
alpha =1;
while(fin==0)
    if (strcmp(version,'our-vs0')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e0/sqrt(k);
        %ind = unidrnd(n,1,1);
        %s = II(:,ind);
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'our-vs1')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e1/sqrt(k);
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'our-vs2')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e-1/sqrt(k);
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'our-fs0')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = epsilon;
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'our-fs1')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e-1*epsilon;
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'our-fs2')) % our method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e1*epsilon;
        ss = randn(n,1);
        %s= ss/norm(ss);
        s = 1/sqrt(n)*ss;
        %         ind = unidrnd(n,1,1);
        %         s = II(:,ind);
        s_c = alpha*s;
        x1 = x+s_c;
        x2 = x-s_c;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f);
        if (ff(1) == 1)
            x = x1;
            f = f1;
        end
        if (ff(2) == 1)
            x = x2;
            f = f2;
        end
    end
    if (strcmp(version,'roy')) % Royer method
        ss1 = randn(n,1);
        s1= ss1/norm(ss1);
        s_c1 = alpha*s1;
        ss2 = randn(n,1);
        s2= ss2/norm(ss2);
        s_c2 = alpha*s2;
        x1 = x+s_c1;
        x2 = x+s_c2;
        f1 =  feval(une_f,x1);
        f2 =  feval(une_f,x2);
        fcount = fcount +2;
        ff = ([f1 f2]<= f - 1e0*alpha^2);
        if (ff(1) == 1)
            alpha =2*alpha;
            x = x1;
            f = f1;
        elseif (ff(2) == 1)
            alpha =2*alpha;
            x = x2;
            f = f2;
        else
            alpha =alpha/2;
        end
    end
    if (strcmp(version,'nes0')) % Nesterov method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e0/(4*(n+4));
        ss = randn(n,1);
        s= ss/norm(ss);
        ff = feval(une_f,x+1e-4*s);
        gs = (ff - f)/1e-4;
        s_c = alpha*gs*s;
        x = x-s_c;
        f =  feval(une_f,x);
        fcount = fcount +2;
    end
    if (strcmp(version,'nes1')) % Nesterov method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e-1/(4*(n+4));
        ss = randn(n,1);
        s= ss/norm(ss);
        ff = feval(une_f,x+1e-4*s);
        gs = (ff - f)/1e-4;
        s_c = alpha*gs*s;
        x = x-s_c;
        f =  feval(une_f,x);
        fcount = fcount +2;
    end
    if (strcmp(version,'nes2')) % Nesterov method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e-2/(4*(n+4));
        ss = randn(n,1);
        s= ss/norm(ss);
        ff = feval(une_f,x+1e-6*s);
        gs = (ff - f)/1e-6;
        s_c = alpha*gs*s;
        x = x-s_c;
        f =  feval(une_f,x);
        fcount = fcount +2;
    end
    if (strcmp(version,'nes3')) % Nesterov method
        %alpha = muP/L*sqrt(2*lamda*(f-fm));
        alpha = 1e1/(4*(n+4));
        ss = randn(n,1);
        s= ss/norm(ss);
        ff = feval(une_f,x+1e-6*s);
        gs = (ff - f)/1e-6;
        s_c = alpha*gs*s;
        x = x-s_c;
        f =  feval(une_f,x);
        fcount = fcount +2;
    end
    if (strcmp(version,'dds')) % DDS method
        j =1;
        while( j < 2*n)
            ss = zeros(n,1);
            if (j<=n)
                ss(j,1) = 1;
            else
                ss(j-n,1) = -1;
            end
            s_c = alpha*ss;
            x1 = x-s_c;
            f1 =  feval(une_f,x1);
            fcount = fcount +1;
            if (f1<= f - 1e0*alpha^2)
                alpha =2*alpha;
                x = x1;
                f = f1;
                j=0;
            end
            j = j+1;
        end
        if (j== 2*n)
            alpha =alpha/2;
        end
    end
    vecfun = [vecfun,f];
    if(((f-fm)/(f0-fm)) < epsilon)
        fin=1;
    end
    if(k> un_nit_max)
        fin=2;
    end
    if(fcount > max_feval)
      fin=3;
    end
    k = k +1;
end
f_opt=f;
x_opt=x;
end
