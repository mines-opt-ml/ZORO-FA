function Result = SMTP(fparam, param, algname, momentum)
% Search direction found in the orthogonal complement of the previous
% direction(s)

%% INITIALIZATION
Result = struct;
n = param.n;
% eps = 1e-6; maxit = 100;

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
objval_seq(1) = f(x); %initialization
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;


%% ITERATION


mu_seq = zeros(maxit+1,1);

beta = 0.5;
v = zeros(n,1);
EPS = param.eps_bergou;
if (strcmp(algname, 'STP-vs'))
    version = 'our-vs0';
elseif (strcmp(algname, 'STP-fs'))
    version = 'our-fs1';
end
if ~strcmp(algname,'SMTP')
    [x_opt, f_opt, fin, k, vecfun, fcount] ...
        = RandomDS_algo(f,objval_seq(1), x, fmin, maxit, param.MAX_QUERIES, EPS, version, 0, 0, 0);
    x = x_opt; % sol
    objval_seq = vecfun';
    num_iter = k;
    num_queries = fcount;
    converged = (fin==1);
    
else % SMTP
    for k=1:maxit
        
        num_queries(k+1) = num_queries(k);
        
        fx = objval_seq(k); % use saved value
        u = PickRandDir(1, n, randAlg)'; % first choose std normal
        u = u/norm(u);
        if param.fixed_step
            gamma = 0.1*param.eps_bergou;
        else
            gamma = 1/sqrt(k+1);
        end
        
        %     gamma = 1e-2/sqrt(k+1);
        
        if momentum
            vp = beta*v + u;
            vm = beta*v - u;
            xp = x - gamma*vp;
            xm = x - gamma*vm;
            zp = xp - gamma*beta/(1-beta)*vp;
            zm = xm - gamma*beta/(1-beta)*vm;
            
            fp = f(zp);
            fm = f(zm);
            
            [fx, idx] = min([fx, fp, fm]);
            if idx==1
                % no updates
                %             sol_seq{k+1} = sol_seq{k}; % same as previous
            elseif idx==2
                x = xp;
                v = vp;
                %             sol_seq{k+1} = zp; % record new iterate
            elseif idx==3
                x = xm;
                v = vm;
                %             sol_seq{k+1} = zm; % record new iterate
            end
            
        else
            fp = f(x+gamma*u);
            fm = f(x-gamma*u);
            [fx, idx] = min([fx, fp, fm]);
            if idx==1
                delta = 0;
            elseif idx==2
                delta = gamma*u;
            elseif idx==3
                delta = -gamma*u;
            end
            x = x+delta;
        end
        mu_seq(k+1) = gamma;
        fxnew = fx;
        objval_seq(k+1) = fx;
        num_queries(k+1) = num_queries(k+1) + 2;
        
        if isfield(fparam, 'fmin')
            if (fx < fmin + eps)
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
    end
    
    if (k==maxit) || (max(num_queries)>param.MAX_QUERIES)
        if verbose>1
            disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
        end
        converged = false;
    end
    num_iter = k;
    objval_seq = objval_seq(1:num_iter+1);
    % sol_seq = sol_seq(1:num_iter+1);
    % gamma_seq = gamma_seq(1:num_iter+1);
    num_queries = num_queries(1:num_iter+1);
    mu_seq = mu_seq(1:num_iter+1);
    % ndelta = ndelta(1:num_iter);
    
    % put into a struct for output
    if isfield(param,'save_x')
        if param.save_x
            %         Result.sol = sol_seq;
        end
    end
end

Result.objval_seq = objval_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.converged = converged;
Result.mu_seq = mu_seq;
Result.sol = x;
% Result.ndelta = ndelta;
end





