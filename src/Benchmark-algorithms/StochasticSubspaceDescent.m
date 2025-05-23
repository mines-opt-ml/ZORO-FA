function Result = StochasticSubspaceDescent(fparam, param)
%        Result = StochasticSubspaceDescent(fparam, param)
% An implementation of the stochastic subspace descent algorithm analyzed
% in "A stochastic subspace approach to gradient-free optimization in high
% dimensions" by Kozak et al (2021).

algname = 'StochasticSubspaceDescent';
%% INITIALIZATION
Result = struct;
Result.algname = algname;
n = param.n;
eps = 1e-10; 

if isfield(param, 'eps')
    eps = param.eps;
end
if isfield(param, 'x0')
    x = param.x0;
else
    x = zeros(n,1); % default initial solution.
end
if isfield(param, 'maxit')
    maxit = param.maxit;
else
    maxit = 100;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
else
    verbose = false; %by default turn verbose mode off.
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
if isfield(param, 'early_stopping')
    early_stopping = param.early_stopping;
else
    early_stopping = true;
end
if isfield(param, 'num_samples')
    num_samples = param.num_samples;
else
    num_samples = n;
end
if isfield(param, 'step_size_SSD')
    alpha = param.step_size_SSD;
    alpha = (num_samples/n)*alpha; % adjust for subsampling.
else
    alpha = 0.01;
end
if isfield(param, 'delta') %sampling radius
    mu = param.delta;
else
    mu = 10^-5;
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
gamma_seq = zeros(maxit+1,1);
if requires_params
    objval_seq(1) = f(x, fparam); %initialization
else
    objval_seq(1) = f(x); %initialization
end
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;

for k=1:maxit
    num_queries(k+1) = num_queries(k);
    fx = objval_seq(k); % use saved value
    delta = zeros(length(x),1); % the gradient approximation
    % Sample num_samples orthogonal directions. We use the approach
    % suggested by Kozak et al on pg. 348.
    Z = randn(n, num_samples);
    [Q,~] = qr(Z);
    for i=1:num_samples
        u = sqrt(n/num_samples)*Q(:, i);
        if requires_params
            fp = f(x+mu*u, fparam);
        else
            fp = f(x+mu*u); 
        end
        num_queries(k+1) = num_queries(k+1) + 1;
        d = (fp - fx)/mu; % directional derivative
        delta = delta + d*u;
    end
    
    delta = -alpha*delta; % move to the next iterate
    
    % Update
    if requires_params
        fxnew = f(x+delta, fparam);
    else
        fxnew = f(x+delta);
    end
    num_queries(k+1) = num_queries(k+1) + 1; %single query here
    
    x = x + delta;

    objval_seq(k+1) = fxnew;
    disp(['SSD. Obj val is ', num2str(fxnew), ' and num queries is ', num2str(num_queries(k+1))])
    
    if isfield(fparam, 'fmin') && early_stopping == true
        if (fxnew < fmin + eps) 
            if verbose
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
                num_iter = k;
            end
            converged = true;
            break;
        end
    end
    if (num_queries(k+1)>param.budget)
        if verbose
            disp('Max queries hit!')
        end
        converged = false;
        num_iter = k;
        break;
    end
end

if (k>=maxit)
    if verbose
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

objval_seq = objval_seq(1:num_iter+1);
gamma_seq = gamma_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);

% put into a struct for output
Result.objval_seq = objval_seq;
Result.gamma_seq = gamma_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.converged = converged;
Result.sol = x;
end
