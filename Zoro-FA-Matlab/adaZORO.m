function Result = adaZORO(fparam, param)
%        output = ZORO(function_handle, function_params, params)
%
% Adaptive version of ZORO, as proposed in original ZORO paper.
% ============================ Inputs ================================== %
% function_handle ....... Handle of zeroth order oracle. Must be preceded
%                         by an "@"
% params ................ Struct containing the parameters: sparsity,
%                         sampling radius (delta), x0 (initial point), 
%                         num_iters and step_size                    
%
% ============================= Outputs ================================ %
% output ............... Struct containing all output data, specifically:
%
% Geovani Luis Grapiglia and Daniel McKenzie
% March 2022

%% INITIALIZATION
Result = struct;
Result.algname = 'adaZORO';
n = param.n;
eps = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
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
if requires_params
    objval_seq(1) = f(x, fparam); %initialization
else
    objval_seq(1) = f(x); %initialization
end

num_queries = zeros(maxit+1,1);
num_queries(1) = 1;  % evaluated f(x) at initialization

% Algorithm specific parameters.
sparsity = param.sparsity;
delta = param.delta;
step_size = param.step_size;
phi = 2.5e-1; %tolerance parameter. Hardcoding this for now.

for k=1:maxit
    % Try computing gradient using support of previous grad estimate
    % if k > 1
    %     support_gk = find(grad_estimate);
    %     column_mask = eye(n);
    %     column_mask = column_mask(:, support_gk);
    %     y = [];
    %     for j = 1:sparsity
    %         new_y_val = (f(x + delta*Z(:,j), fparam) - f(x, fparam))/delta;
    %         y = [y; new_y_val];
    %     end
    %     lsqr_grad_estimate = Z(support_gk, 1:s)\y;
    % end
    flag = false;
    j = 0;
    s_j = sparsity;
    num_queries_this_iter = 0;
    while flag == false
        num_samples = min(n,2*ceil(s_j*log2(n))); % taking b1 = 2 for now. USED TO BE log2(n/s_j) but this is prone to errors.
        Z =(2*(rand(num_samples,n) > 0.5) - 1)/sqrt(num_samples);  % Generate Rademacher sampling vecs
        cosamp_params.Z = Z;
        cosamp_params.sparsity = s_j;
        cosamp_params.delta = delta;
        cosamp_params.n = 20; % Hardcode number of iterations of CoSaMP.
        cosamp_params.x = x;
        try
            [f_k,grad_estimate,err] = CosampGradEstimate(f,fparam, cosamp_params);
        catch ME
            disp('An error has occurred. Terminating run of adaZORO')
            objval_seq(k+1) = f_k;
            % Package and return results
            num_iter = k;
            num_queries(k+1) = num_queries(k) +num_queries_this_iter;
            objval_seq = objval_seq(1:num_iter+1);
            num_queries = num_queries(1:num_iter+1);

            % package
            Result.objval_seq = objval_seq;
            Result.num_queries = num_queries;
            Result.sol = x;
            Result.converged = false;
            return
        end
        num_queries_this_iter = num_queries_this_iter + num_samples;
        disp(['Num queries thus far', num2str(num_queries(k) +num_queries_this_iter)])

        % == Check to see if we have hit feval budget.
        if num_queries(k) +num_queries_this_iter >param.budget
            objval_seq(k+1) = f_k;
            disp('Max queries hit!')
            % Package and return results
            num_iter = k;
            num_queries(k+1) = num_queries(k) +num_queries_this_iter;
            objval_seq = objval_seq(1:num_iter+1);
            num_queries = num_queries(1:num_iter+1);

            % package
            Result.objval_seq = objval_seq;
            Result.num_queries = num_queries;
            Result.sol = x;
            Result.converged = false;
            return
        end

        % Check residual to see if gradient estimate is accepted
        %disp(['Sampling data error is ', num2str(sampling_data.err)])
        if err < phi
            x = x - step_size*grad_estimate;
            objval_seq(k+1) = f_k;
            num_queries(k+1) = num_queries(k) + num_queries_this_iter;
            disp(['Current function value is ', num2str(f_k), ' and num samples is ', num2str(num_queries(k+1))])
            flag = true;
        end
        j = j+1;
        s_j = s_j + 10;
    end
end

if (k>=maxit)
    if verbose
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

% Package and return results
num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);

% package
Result.objval_seq = objval_seq;
Result.num_queries = num_queries;
Result.sol = x;
Result.converged = converged;
    
    
