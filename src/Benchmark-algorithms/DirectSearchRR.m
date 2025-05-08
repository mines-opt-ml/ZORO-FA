function Result = DirectSearchRR(fparam, param)
%        Result = DirectSearchRR(fparam, param)
% An implementation of the probabilistic direct search algorithm of 
% "Direct Search Based on Probabilistic Descent in Reduced Subspaces" by 
% Roberts and Royer (2023).
% We implement the most performant variant: r=1 and Gaussian P_k

algname = 'DirectSearchRR';
%% INITIALIZATION
Result = struct;
Result.algname = algname;
n = param.n;
eps = 1e-6; maxit = 100;

if isfield(param, 'eps')
    eps = param.eps;
end
if isfield(param, 'x0')
    x = param.x0;
else
    x = zeros(n,1);
end
if isfield(param, 'maxit')
    maxit = param.maxit;
end
if isfield(param, 'randAlg')
    randAlg = param.randAlg;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
else
    verbose = false;
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
if isfield(param, 'num_samples')
    num_samples = param.num_samples;
else
    num_samples = 1;
end
if num_samples > 1
    disp('Sorry, we have not implemented that functionality')
    return
end
if isfield(param, 'alpha0')
    alpha0 = param.alpha0;
else
    alpha = 1.0; %value used in experiments of Roberts and Royer
end
if isfield(param, 'gamma_inc')
    gamma_inc = param.gamma_inc;
else
    gamma_inc = 2; %value used in experiments of Royer and Roberts
end
if isfield(param, 'gamma_dec')
    gamma_dec = param.gamma_dec;
else
    gamma_dec = 0.5; % value used in experiments of Roberts and Royer.
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
if requires_params
    objval_seq(1) = f(x, fparam); %initialization
else
    objval_seq(1) = f(x); %initialization
end
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;

for k=1:maxit
    num_queries(k+1) = num_queries(k);
    mu =  1/sqrt(k+1);
    fx = objval_seq(k); % use saved value
    delta = zeros(length(x),1); % the gradient approximation
    u = randn(n,1);
    x_plus = x + alpha_k*u;
    x_minus = x - alpha_k*u;
    target = c*alpha_k^2*norm(u,2)^2/2;
    if requires_params
        f_plus = f(x_plus, fparam);
        r_minus = f(x_minus, fparam);
    else
        f_plus = f(x_plus);
        f_minus = f(x_minus, fparam);
    end
    % DM stopped here.
    
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
