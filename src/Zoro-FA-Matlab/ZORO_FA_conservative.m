function Result = ZORO_FA_conservative(fparam, param)
% Implementation of ZORO with Line Search designed to work with
% MGH benchmarking code.
% "Conservative" version, which means that by default it does not assume
% gradient vectors to be sparse.
% Daniel McKenzie and Geovani Nunes Grapiglia.
% April 2025.

algname = 'ZORO-FA-conservative';

%% INITIALIZATION
Result = struct;
Result.algname = algname;
n = param.n;
target_reduction_frac = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
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
if isfield(param, 'verbose')
    verbose = param.verbose;
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
if isfield(param, 'b')
    b = param.b;
else
    b = 1; % default value
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
sigma_0 = param.sigma0;
s0 = param.sparsity; %ToDo: Do I need this?

tau = 11; % Constant coming from the compressed sensing theory.
rho = 0.5; % From compressed sensing theory.
b = 1; % From compressed sensing theory.
% Step 0
Z = 2*((rand(n,n) > 0.5) - 1); % generate Rademacher vectors

% Initialize extra parameters for the conservative variant
m_k = n; % the default number of sampling directions to use.
sigma_k = sigma_0;
L_k = sigma_0; 
s_k = n;
count_successful_sparse_gradients = 0;

for k=1:maxit
    % Step 1
    num_queries(k+1) = num_queries(k);
    if m_k >= n
        flag = false;
        % == Step 2
        j = max(ceil(log2(sigma_0/sigma_k)) + 1, 0);
        while flag == false 
            % == Step 2.1
            % compute a non-sparse gradient
            sigma = (2^j)*sigma_k; 
            h = 2*theta*epsilon/(sigma*sqrt(n));
            % call the function which computes a gradient approximation.
            fd_params.delta = h;
            fd_params.x = x;
            [f_k,grad_estimate] = FDGradEstimate(f, fparam, fd_params); 
            num_queries(k+1) = num_queries(k+1) + n + 1;

            % == Step 2.2
            x_plus = x - (1/(sigma))*grad_estimate;
            if requires_params
                f_plus = f(x_plus, fparam);
            else
                f_plus = f(x_plus);
            end

            if (f_k - f_plus >= epsilon^2/(2*sigma))
                % == Step 3
                [s_flag, s] = check_sparsity(grad_estimate, epsilon);
                x = x_plus;
                sigma_k = sigma/2;
                L_k = max(L_k, sigma);
                s_k = s;
                m_k = ceil(b*s_k*log(n));

                % various flags and outputs. 
                flag = true;
                disp(['Step size is ', num2str((1/(sigma)))])
                disp(['Current function value is ', num2str(f_plus)])
                disp(['jk for this iteration is ', num2str(j)])
                disp(['Dimension is ',num2str(n), ' and current gradient sparsity is ',num2str(s_k)])

                objval_seq(k+1) = f_plus;
                sparse_seq(k+1) = s;
            else 
                j = j+1;
            end
            % check to see if maximum queries reached
            disp(['Alg is ZORO_FA_conservative and number of queries is ', num2str(num_queries(k+1))])
            if (num_queries(k+1)>param.budget)
                objval_seq(k+1) = f_k;
                sparse_seq(k+1) = s_k;
                disp(['Max queries hit!'])
                % Package and return results
                num_iter = k;
                objval_seq = objval_seq(1:num_iter+1);
                % sol_seq = sol_seq(1:num_iter+1);
                sparse_seq = sparse_seq(1:num_iter+1);
                num_queries = num_queries(1:num_iter+1);
                disp(['ZORO_FA_C number of successful sparse gradients is ', num2str(count_successful_sparse_gradients)])
    
                % package
                Result.objval_seq = objval_seq;
                Result.sparse_seq = sparse_seq;
                Result.num_queries = num_queries;
                Result.sol = x;
                Result.converged = false;
                return
            end
        end
    else % if m < n
        % == Step 4
        h = (theta*epsilon)/(n*tau*L_k);
        cosamp_params.Z = Z(1:m_k,:)/sqrt(m_k);
        cosamp_params.delta = max(h,1e-4); % add a floor to prevent delta too small.
        cosamp_iters = ceil(log(theta/4)/log(rho));
        cosamp_params.n = cosamp_iters;
        cosamp_params.x = x;
        cosamp_params.sparsity = s_k;

        % == Step 5
        [f_k, grad_estimate, min_query_point, perturbation_yielding_min_query, err] = CosampGradEstimate(f, fparam, cosamp_params);
        num_queries(k+1) = num_queries(k+1) + m_k + 1;

        % == Step 6
        x_plus = x - (1/L_k)*grad_estimate; % new trial point.
        if requires_params
            f_plus = f(x_plus, fparam);
        else
            f_plus = f(x_plus);
        end
        if (f_k - f_plus) >= epsilon^2/(2*L_k)
            disp(['Step size is ', num2str((1/(L_k)))])
            disp(['Current function value is ', num2str(f_plus)])
            disp(['ZORO_FA_C: Dimension is ',num2str(n), ' and current target sparsity is ',num2str(s_k) , ' and m_k is ', num2str(m_k)])
            x = x_plus;
            objval_seq(k+1) = f_plus;
            sparse_seq(k+1) = s_k;
            count_successful_sparse_gradients = count_successful_sparse_gradients + 1; 
        else
            [min_func_value, idx_min_func_value] = min([f_k, f_plus, min_query_point]);
            if idx_min_func_value == 2
                x = x_plus;
            elseif idx_min_func_value == 3
                x = x + perturbation_yielding_min_query;
            end
            disp('Failed sparse gradient step')
            m_k = n;
        end
         % check to see if maximum queries reached
        if (num_queries(k+1)>param.budget)
            objval_seq(k+1) = f_k;
            sparse_seq(k+1) = s;
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
            disp(['ZORO_FA_C number of successful sparse gradients is ', num2str(count_successful_sparse_gradients)])
            return
        end
    end
end


        
    