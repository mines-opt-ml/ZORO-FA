function Result = VRNS(fparam, param, NQ)
algname = 'VRNS';
%% INITIALIZATION
Result = struct;
n = param.n;
eps = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
randAlg = 'U'; % uniform distribution
verbose = 0; % default setting

%% Variance Reduction (VR) options
N_vr = n;
vr_threshold = 10*n;
%%
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
if isfield(fparam, 'N_vr')
    fmin = fparam.N_vr;
end
f = fparam.f;
objval_seq = zeros(maxit+1,1);
gamma_seq = zeros(maxit+1,1);
objval_seq(1) = f(x); %initialization
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;
damped_cnt = 0;

%% ITERATION

% "mu" is denoted "r" in the paper

% mu_init = 1e-2; % More
mu_init = 5e-1; % osc
mu = mu_init;

mu_seq = zeros(maxit+1,1);
mu_seq(1) = mu;
alpha = 0.5; % step size param, 1/Lhat in the paper

CARScounter = zeros(4,1);

VR_fd_method = 2; % 1 for forward, 2 for central difference
g_vr = zeros(n,1);
for k=1:maxit
    num_queries(k+1) = num_queries(k);
    mu =  1/sqrt(k+1);
    fx = objval_seq(k); % use saved value
    
    if mod(k,N_vr) == 1 && k > vr_threshold
        [g_vr, f_queries, x_queries] = gradest(f, x, mu, fx, VR_fd_method);
        num_queries(k+1) = num_queries(k+1) + VR_fd_method*n;
        delta = -alpha*g_vr;
    else
        u = PickRandDir(1, n, randAlg)';
        u = u/norm(u); % normalize u
        fp = f(x+mu*u);
        num_queries(k+1) = num_queries(k+1) + 1;
        d = (fp - fx)/mu; % directional derivative
        f_queries = fp;
        x_queries = mu*u;
        delta = -alpha*( (d-dot(g_vr,u))*u + g_vr); % move to the next iterate
    end
    fxnew = f(x+delta);
    num_queries(k+1) = num_queries(k+1) + 1; %single query here
    
    fs = [fx, f_queries, fxnew];
    [fxnew, midx] = min(fs);
    if midx == 1
        % no update
        delta = 0;
    else
        if midx == 2
            delta = x_queries;
        elseif midx == 3
            % delta not changed
        end
    end
    x = x + delta;
    
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

if (k>=maxit) || (num_queries(k+1)>param.MAX_QUERIES)
    if verbose>1
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

% polish the output
num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
gamma_seq = gamma_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);
mu_seq = mu_seq(1:num_iter+1);

% put into a struct for output
Result.objval_seq = objval_seq;
Result.gamma_seq = gamma_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.converged = converged;
Result.damped_cnt = damped_cnt;
Result.mu_seq = mu_seq;
Result.sol = x;
if verbose>1 && NQ==0
    disp(CARScounter');
end
end

function [g, fbest, xbest] = gradest(f, x, r, fx, npts)
d = length(x); % 1-d vector is allowed for x
E = eye(d);
g = zeros(d,1);
fp = zeros(d,1);
if npts>1
    fm = zeros(d,1);
end
fbest = Inf;
for i=1:d
    fp(i) = f(x + r*E(:,i));
    if fp(i) < fbest
        xbest = r*E(:,i);
        fbest = fp(i);
    end
    if npts > 1
        fm(i) = f(x-r*E(:,i));
        g(i) = (fp(i)-fm(i))/2/r;
        if fm(i) < fbest
            xbest = -r*E(:,i);
            fbest = fm(i);
        end        
    else
        g(i) = (fp(i)-fx)/r;
    end
end
end
