if strcmp(fname, 'band')
    m = n; x0 = -(10^s)*ones(n,1); fmin = 0; %true_sol = [1.098e-5, 9.106]';
    f_M = @(x) band(n,m,x,1);
elseif strcmp(fname, 'bv')
     m = n;
    tmp = (1:n)/(n+1); x0 = zeros(n,1);
    for test_idx=1:n
        x0(test_idx) = tmp(test_idx)*(tmp(test_idx)-1);
    end
    x0 = (10^s)*x0;
    fmin = 0;
    f_M = @(x) bv(n,m,x,1);
elseif strcmp(fname, 'ie')
    m = n;
    tmp = (1:n)/(n+1); x0 = zeros(n,1);
    for test_idx=1:n
        x0(test_idx) = tmp(test_idx)*(tmp(test_idx)-1);
    end
    x0 = (10^s)*x0;
    fmin = 0;
    f_M = @(x) ie(n,m,x,1);
elseif strcmp(fname, 'lin')
    m = n; x0 = ones(n,1); fmin = m-n;
    x0 = (10^s)*x0;
    f_M = @(x) lin(n,m,x,1);
elseif strcmp(fname, 'lin0')
    m = n; x0 = ones(n,1); fmin = (m*m+3*m-6)/2/(2*m-3);
    x0 = (10^s)*x0;
    f_M = @(x) lin0(n,m,x,1);
elseif strcmp(fname, 'lin1')
    m = n; x0 = ones(n,1); fmin = (m*(m-1))/2/(2*m+1);
    x0 = (10^s)*x0;
    f_M = @(x) lin1(n,m,x,1);
elseif strcmp(fname,'pen1')
    m = n+1; x0 = (1:n)'; fmin = 7.08765e-5;
    x0 = (10^s)*x0;
    f_M = @(x) pen1(n,m,x,1);
elseif strcmp(fname, 'pen2')
    m = 2*n; x0 = 0.5*ones(n,1); fmin = 2.93660e-4;
    x0 = (10^s)*x0;
    f_M = @(x) pen2(n,m,x,1);
elseif strcmp(fname, 'rosex')
    m = n; 
    x0 = zeros(n,1);
    for test_idx=1:n
        if rem(test_idx,2)==0
            x0(test_idx) = 1;
        else
            x0(test_idx) = -1.2;
        end
    end
    x0 = (10^s)*x0;
    fmin = 0;
    f_M = @(x) rosex(n,m,x,1);
elseif strcmp(fname, 'singx')
    m = n; 
    x0 = zeros(n,1);
    for test_idx=1:n
        if rem(test_idx,4)==1
            x0(test_idx) = 3;
        elseif rem(test_idx,4)==2
            x0(test_idx) = -1;
        elseif rem(test_idx,4)==3
            x0(test_idx) = 0;
        elseif rem(test_idx,4)==4
            x0(test_idx) = 1;
        end
    end
    x0 = (10^s)*x0;
    fmin = 0;
    f_M = @(x) singx(n,m,x,1);
elseif strcmp(fname, 'trid')
    m = n; x0 = -ones(n,1); fmin = 0;
    x0 = (10^s)*x0;
    f_M = @(x) trid(n,m,x,1);
elseif strcmp(fname, 'trig')
    m = n; x0 = 1/n*ones(n,1); fmin = 0;
    x0 = (10^s)*x0;
    f_M = @(x) trig(n,m,x,1);
elseif strcmp(fname, 'vardim')
    m = n+2; x0 = 1 - (1:n)'/n; fmin = 0;
    x0 = (10^s)*x0;
    f_M = @(x) vardim(n,m,x,1);
elseif strcmp(fname, 'bl')
    m = n; x0 = 0.5*ones(n,1); fmin = 0;
    x0 = (10^s)*x0;
    f_M = @(x) bl(n,m,x,1);
end


fx0 = sum(f_M(x0).^2);
f = @(x) sum(f_M(x).^2);
fparam = struct;
fparam.f = f;
fparam.fmin = fmin;
fparam.name = fname;

param = struct;
param.n = n;
param.x0 = x0;
param.maxit = maxit;
param.budget = (n+1)*budget;