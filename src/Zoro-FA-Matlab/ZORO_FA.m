function Result = ZORO_FA(fparam, param)
% Implementation of ZORO with Line Search designed to work with
% MGH benchmarking code.
% Daniel McKenzie and Geovani Nunes Grapiglia
% July 2023

algname = 'ZORO-FA';

%% INITIALIZATION
Result = struct;
Result.algname = algname;
n = param.n;
target_reduction_frac = 1e-100; maxit = 100;

x = zeros(n,1); % initial sol (x0)
verbose = false;

if isfield(param, 'target_reduction_frac')
    target_reduction_frac = param.target_reduction_frac;
end
if isfield(param, 'x0')
    x = param.x0;
end
if isfield(param, 'early_stopping')
    early_stopping = param.early_stopping;
else
    early_stopping = true;
end
if isfield(param, 'maxit')
    maxit = param.maxit;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
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
if requires_params
    objval_seq(1) = f(x, fparam); %initialization
else
    objval_seq(1) = f(x); %initialization
end
f0 = objval_seq(1);

sparse_seq = zeros(maxit+1,1); % track effective sparsity
sparse_seq(1) = 1; % just declare = 1 at beginning
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;  % evaluated f(x) at initialization

% Algorithm specific parameters.
epsilon = param.epsilon;
theta = param.theta;
sigma0 = param.sigma0;
s0 = param.sparsity;

tau = 11; % Constant coming from the compressed sensing theory.
rho = 0.5; % From compressed sensing theory.
% Step 0
Z = 2*((rand(n,n) > 0.5) - 1); % generate Rademacher vectors

for k=1:maxit
    % Step 1
    j = 0;
    flag = false;
    num_queries(k+1) = num_queries(k);
    while flag == false
        % == Step 2.1
        s_j = (2^j)*s0;
        sigma_j = (2^j)*sigma0;
        m_j = ceil(b*s_j*log2(n)); 
        disp(['mj is ', num2str(m_j)])

        % == Step 2.2 and 2.3
        if  m_j <= n
            cosamp_params.sparsity = s_j;
            try
                cosamp_params.Z = Z(1:m_j,:)/sqrt(m_j);
            catch ME
                disp(['Value of m_j is ', num2str(m_j)])
            end
            h_j = (theta*epsilon)/(n*tau*sigma_j);
            cosamp_params.delta = max(h_j,1e-4); % add a floor to prevent delta too small.
            cosamp_iters = ceil(log(theta/4)/log(rho));
            cosamp_params.n = cosamp_iters;
            cosamp_params.x = x;
            [f_k,grad_estimate] = CosampGradEstimate(f, fparam, cosamp_params);
            num_queries(k+1) = num_queries(k+1) + m_j + 1;
        % == Step 3
        else
            h_j = 2*theta*epsilon/(sigma_j*sqrt(n));
            fd_params.delta = h_j;
            fd_params.x = x;
            [f_k,grad_estimate] = FDGradEstimate(f, fparam, fd_params); 
            num_queries(k+1) = num_queries(k+1) + n + 1;
        end
        
        % == Step 4
        x_plus = x - (1/(sigma_j))*grad_estimate;
        if requires_params
            f_plus = f(x_plus, fparam);
        else
            f_plus = f(x_plus);
        end

        if (f_k - f_plus >= epsilon^2/(2*sigma_j))
            flag = true;
            disp(['Step size is ', num2str((1/(sigma_j)))])
            disp(['ZORO-FA. Obj val is ', num2str(f_plus) , 'num queries is ' num2str(num_queries(k+1))])
            disp(['jk for this iteration is ', num2str(j)])
            disp(['Dimension is ',num2str(n), ' and current target sparsity is ',num2str(s_j)])
            x = x_plus;
            objval_seq(k+1) = f_plus;
            sparse_seq(k+1) = s_j;
        else 
            j = j+1;
        end
        if (num_queries(k+1)>param.budget)
            objval_seq(k+1) = f_k;
            sparse_seq(k+1) = s_j;
            disp(['Max queries hit!'])
            % Package and return results
            num_iter = k;
            objval_seq = objval_seq(1:num_iter+1);
            % sol_seq = sol_seq(1:num_iter+1);
            sparse_seq = sparse_seq(1:num_iter+1);
            num_queries = num_queries(1:num_iter+1);

            % package
            Result.objval_seq = objval_seq;
            Result.sparse_seq = sparse_seq;
            Result.num_queries = num_queries;
            Result.sol = x;
            Result.converged = false;
            return
        end
    end
    if isfield(fparam, 'fmin') && early_stopping == true
        tolerance = target_reduction_frac * (f0 - fmin);
        if (f_plus < fmin + tolerance) 
            if verbose
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(f_k)]);
            end
            converged = true;
            break;
        end
    end
end

if (k>=maxit) || (num_queries(k+1)>param.maxit)
    if verbose>1
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

% Package and return results
num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
% sol_seq = sol_seq(1:num_iter+1);
sparse_seq = sparse_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);

% package
Result.objval_seq = objval_seq;
Result.sparse_seq = sparse_seq;
Result.num_queries = num_queries;
Result.sol = x;
Result.converged = converged;