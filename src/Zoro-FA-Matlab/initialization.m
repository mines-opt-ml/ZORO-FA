% Script that initializes all parameters for the algorithms in this
% benchmarking suite. Constructing a separate script ensures consistency

Result = struct;
Result.algname = algname;
n = param.n;

if isfield(param, 'target_reduction_frac')
    target_reduction_frac = param.target_reduction_frac;
end
if isfield(param, 'x0')
    x = param.x0;
else 
    x = zeros(n,1); % initial guess.
end
if isfield(param, 'early_stopping')
    early_stopping = param.early_stopping;
else
    early_stopping = true;
end
if isfield(param, 'maxit')
    maxit = param.maxit;
else
    maxit = 100;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
else
    verbose = false;
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
if isfield(param, 'b')
    b = param.b;
else
    b = 1.0; % default value
end

% Next if statement allows for functions which require further 
% parameters at evaluation.
if isfield(fparam, 'requires_params')
    requires_params = fparam.requires_params;
else
    requires_params = false;
end
f = fparam.f;

% Some arrays to store results
objval_seq = zeros(maxit+1,1);
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;  % evaluated f(x) at initialization
if requires_params
    objval_seq(1) = f(x, fparam); %initialization
else
    objval_seq(1) = f(x); %initialization
end