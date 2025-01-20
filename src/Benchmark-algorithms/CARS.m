function Result = CARS(fparam, param)
%        Result = CARS(fparam, param)
% The CARS algorithm, as proposed in "Curvature-aware Derivative-free
% optimization" by Kim et al (2023).
% Bumsu Kim 2021
% Edited and simplified by Daniel McKenzie 2023

algname = 'CARS';
%% INITIALIZATION
Result = struct;
Result.algname = algname;
n = param.n;
target_reduction_frac = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
randAlg = 'U'; % uniform
verbose = false;


if isfield(param, 'target_reduction_frac')
    target_reduction_frac = param.target_reduction_frac;
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
f0 = objval_seq(1);

num_queries = zeros(maxit+1,1);
num_queries(1) = 1;
damped_cnt = 0;

%% ITERATION

mu_init = 1e-2; % More
mu = mu_init;
mu_cnt = 0;
gd_cnt = 0;

mu_seq = zeros(maxit+1,1);
mu_seq(1) = mu;
alpha = 0.5;
CARScounter = zeros(4,1);

for k=1:maxit
    num_queries(k+1) = num_queries(k);
    mu = mu* sqrt(k)/sqrt(k+1);
    fx = objval_seq(k); % use saved value
    u = PickRandDir(1, n, randAlg)'; % first choose std normal
    u = u/norm(u); % normalize u

    if requires_params
        fp = f(x+mu*u, fparam);
        fm = f(x-mu*u, fparam);
    else
        fp = f(x+mu*u);
        fm = f(x-mu*u);
    end
    num_queries(k+1) = num_queries(k+1) + 2;
    d = (fp - fm)/(2*mu); % directional derivative
    
    h = (fp + fm - 2*fx)/mu^2; % 2nd order dir deriv
    delta = -alpha*d/h*u; % move to the next iterate
    if requires_params
        fxcars = f(x+delta, fparam);
    else
        fxcars = f(x+delta);
    end
    num_queries(k+1) = num_queries(k+1) + 1; %single query here
    
    fs = [fx, fp, fm, fxcars];
    [fxnew, midx] = min(fs);
    CARScounter(midx) = CARScounter(midx)+1;
    if midx == 1
        % no update
        delta = 0;
    else
        if midx == 2
            delta = mu*u;
        elseif midx == 3
            delta = -mu*u;
        elseif midx == 4
            % delta not changed
        end 
    end
    
    x = x + delta;
    objval_seq(k+1) = fxnew;
     
    if isfield(fparam, 'fmin')
        tolerance = target_reduction_frac * (f0 - fmin);
        if (fxnew < fmin + tolerance) 
            if verbose
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
            end
            converged = true;
            break;
        end
    end
    if (num_queries(k+1)>param.budget)
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
