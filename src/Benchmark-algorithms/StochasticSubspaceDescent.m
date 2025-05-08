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
eps = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
randAlg = 'U'; % uniform
verbose = false;


if isfield(param, 'eps')
    eps = param.eps;
end
if isfield(param, 'x0')
    x = param.x0;
end
if isfield(param, 'maxit')
    maxit = param.maxit;
end
if isfield(param, 'randAlg')
    randAlg = param.randAlg;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
if isfield(param, 'num_samples')
    num_samples = param.num_samples;
else
    num_samples = 1;
end
if isfield(param, 'step_size')
    alpha = param.step_size;
else
    alpha = 0.01;
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
damped_cnt = 0;

%% ITERATION
mu_init = 1e-4;
mu = mu_init;
mu_cnt = 0;
gd_cnt = 0;

mu_seq = zeros(maxit+1,1);
mu_seq(1) = mu;

for k=1:maxit
    num_queries(k+1) = num_queries(k);
    mu =  1/sqrt(k+1);
    fx = objval_seq(k); % use saved value
    delta = zeros(length(x),1); % the gradient approximation
    % Sample num_samples orthogonal directions. We use the approach
    % suggested by Kozak et al on pg. 348.
    Z = randn(n, num_samples);
    [Q,R] = qr(Z);
    for i=1:num_samples
        u = sqrt(n/num_samples)*Q(:, i); %check scaling factor is right
        if requires_params
            fp = f(x+mu*u, fparam);
        else
            fp = f(x+mu*u); %initialization
        end
        num_queries(k+1) = num_queries(k+1) + num_samples;
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
    disp(['Current objective function value is ', num2str(fxnew)])
    
    if isfield(fparam, 'fmin')
        if (fxnew < fmin + eps) 
            if verbose
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
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
        break;
    end
    mu_seq(k+1) = mu;
end

if (k>=maxit)
    if verbose
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
gamma_seq = gamma_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);
mu_seq = mu_seq(1:num_iter+1);

% put into a struct for output
if isfield(param,'save_x')
    if param.save_x
%         Result.sol = sol_seq;
    end
end
Result.objval_seq = objval_seq;
Result.gamma_seq = gamma_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.converged = converged;
Result.damped_cnt = damped_cnt;
Result.mu_cnt = mu_cnt;
Result.gd_cnt = gd_cnt;
Result.mu_seq = mu_seq;
Result.sol = x;
end
