function Result = CARS_adaptive_mu(fparam, param, NQ)
algname = 'CARS';
%% INITIALIZATION
Result = struct;
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
f = fparam.f;
objval_seq = zeros(maxit+1,1);
gamma_seq = zeros(maxit+1,1);
objval_seq(1) = f(x); %initialization
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;
damped_cnt = 0;

%% ITERATION
mu_init = 1;
mu_max = 1;
mu_min = 1e-8;
mu = mu_init;
mu_cnt = 0;
gd_cnt = 0;

mu_seq = zeros(maxit+1,1);
mu_seq(1) = mu;
alpha = 0.5;

n_dir_grad = 1; ceil(n/2);
been_there = false;
stuck_count = 0;
CARScounter = zeros(4,1);
for k=1:maxit
    num_queries(k+1) = num_queries(k);
    if mu<mu_min
        mu = mu_min;
    else
        mu = mu* sqrt(k/(k+1));
    end
    fx = objval_seq(k); % use saved value
    %     if stuck_count > n*2
    %         been_there = false;
    %         if mu<mu_max
    %             mu = mu*2;
    %         end
    %         stuck_count = 0;
    %     end
    %     if n_dir_grad == 1 && ~been_there
    u = PickRandDir(1, n, randAlg)'; % first choose std normal
    %         us = zeros(0,n);
    %         fvals = zeros(0,2);
    %     else
    %         if ~been_there
    %             us = PickRandDir(n_dir_grad, n, randAlg);
    %             fvals = Queries(f, x, us, mu);
    %         else
    %             us_1 = PickRandDir(n_dir_grad, n, randAlg);
    %             fvals_1 = Queries(f, x, us_1, mu);
    %             us = [us; us_1];
    %             fvals = [fvals; fvals_1];
    %         end
    %         u = EstGrad( fvals, us, mu, 'P');
    %         num_queries(k+1) = num_queries(k+1) + 2*n_dir_grad;
    %     end
    u = u/norm(u); % normalize u
    if NQ==0 % Regular CARS
        fp = f(x+mu*u);
        fm = f(x-mu*u);
        num_queries(k+1) = num_queries(k+1) + 2;
        d = (fp - fm)/(2*mu); % directional derivative
        
        h = (fp + fm - 2*fx)/mu^2; % 2nd order dir deriv
        
        delta = -alpha*d/h*u; % move to the next iterate
        fxnewton = f(x+delta);
        num_queries(k+1) = num_queries(k+1) + 1; %single query here
        
        fs = [fx, fp, fm, fxnewton];
        [fxnew, midx] = min(fs);
        CARScounter(midx) = CARScounter(midx)+1;
        if midx == 1
            % no update
            delta = 0;
            been_there = true;
            mu = mu/2;
        else
            been_there = false;
            %         mu = mu_init;
            if midx == 2
                delta = mu*u;
                mu = mu/2;
            elseif midx == 3
                delta = -mu*u;
                mu = mu/2;
            elseif midx == 4
                % delta not changed
            end
        end
    else % CARS-NQ
        [delta, fxnew] = Directional_Newton(f, x, u, mu, NQ, true, fx, false);
        num_queries(k+1) = num_queries(k+1) + NQ; % (NQ-1) queries + 1 query
        if any(isnan(delta)) || isnan(fxnew)
            mu = mu/2;
            delta = 0;
            fxnew = fx;
        end
    end
    
    
    
    x = x + delta;
    %     if norm(delta) < 0.99*mu
    %         been_there = true;
    %     end
    %     if been_there || (fxnew-objval_seq(k) < 1e-6)
    %         stuck_count = stuck_count + 1;
    %     else
    %         mu = mu_init;
    %         stuck_count = 0;
    %     end
    %     if mu>1e-8
    %         mu = 0.5* mu + norm(delta) * mu;
    %     end
    %     mu= 1/sqrt(k+2);
    
    if isnan(fxnew) % should never happen..
        disp('errrrrrr');
    end
    %     x = x+delta; % update the solution
    objval_seq(k+1) = fxnew;
    
    if isfield(fparam, 'fmin')
        if (fxnew < fmin + eps)
            if verbose>1
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
            end
            converged = true;
            break;
        end
    end
    if (num_queries(k+1)>param.MAX_QUERIES)
        break;
    end
    mu_seq(k+1) = mu;
end

if (k>=maxit) || (num_queries(k+1)>param.MAX_QUERIES)
    if verbose>1
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
% sol_seq = sol_seq(1:num_iter+1);
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
% Result.k_lambda = k_lambda_condition;
if verbose>1 && NQ==0
    disp(CARScounter');
end
end
