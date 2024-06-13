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
phi = 1e-1; %tolerance parameter. Hardcoding this for now.

% Sundry variables
total_fevals = 0;

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
    num_samples_this_iter = 0;
    num_old_samples = 0;
    clear cosamp_params;
    j = 0;
    s_j = sparsity;
    while flag == false
        num_samples = min(n-1, ceil(2*s_j*log2(n/s_j)));
        num_new_samples = num_samples - num_old_samples;
        num_old_samples = num_samples;
        Z_new = (2*(rand(num_new_samples,n) > 0.5) - 1)/n;
        cosamp_params.sparsity = s_j;
        cosamp_params.Z_new = Z_new;
        cosamp_params.num_samples = num_samples;
        cosamp_params.delta = delta;
        cosamp_params.x = x;
        cosamp_params.n = 20; % Hardcode number of iterations of CoSaMP
        [f_k,grad_estimate,sampling_data] = CosampGradEstimate_query_recycling(f, fparam, cosamp_params);
        num_samples_this_iter = num_samples_this_iter + num_new_samples;
        cosamp_params.y = sampling_data.y;
        cosamp_params.Z = sampling_data.Z;

        % == Check to see if we have hit feval budget.
        total_fevals = total_fevals + num_new_samples;
        if total_fevals>param.budget
            objval_seq(k+1) = f_k;
            disp('Max queries hit!')
            % Package and return results
            num_iter = k;
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
        if sampling_data.err < phi
            x = x - step_size*grad_estimate;
            objval_seq(k+1) = f_k;
            num_queries(k+1) = num_queries(k) + num_samples_this_iter;
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
    
    
