function Result = ZORO(fparam, param)
%        output = ZORO(function_handle, function_params, params)
%
% Vanilla version of ZORO. For comparison purposes.
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
Result.algname = 'ZORO';
n = param.n;
eps = 1e-10; maxit = 100;

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

% ==== Initialize parameters
num_samples = 2*ceil(sparsity*log2(n/sparsity)); % taking b1 = 2 for now.
Z =(2*(rand(num_samples,n) > 0.5) - 1)/sqrt(num_samples);  % Generate Rademacher sampling vecs
cosamp_params.Z = Z;
cosamp_params.sparsity = sparsity;
cosamp_params.delta = delta;
cosamp_params.n = 20; % Hardcode number of iterations of CoSaMP.

for k=1:maxit
    cosamp_params.x = x;
    disp('This is an iteration of ZORO.')
    try
        [f_k,grad_estimate,err] = CosampGradEstimate(f,fparam, cosamp_params);
    catch ME
        disp('An error has occurred. Terminating run of ZORO')
        converged = false;
        break
    end
    x = x - step_size*grad_estimate;
    objval_seq(k) = f_k;
    num_queries(k+1) = num_queries(k) + num_samples;
    disp(['Current function value is ', num2str(f_k), ' and total number of queries is ', num2str(num_queries(k+1))])
    if (num_queries(k+1)>param.budget)
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
    if isfield(fparam, 'fmin')
        % eps = EPS_MORE * (f0 - fmin)
        if (f_k < fmin + eps) 
            if verbose
                disp(['ZORO Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
            end
            converged = true;
            break;
        end
    end
end

if (k>=maxit)
    if verbose
        disp(['ZORO did not converge in ', num2str(maxit) , ' steps.']);
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
    
    
