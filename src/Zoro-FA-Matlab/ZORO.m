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

algname = 'ZORO';

%% INITIALIZATION
initialization;

% Algorithm specific parameters.
sparsity = param.sparsity;
delta = param.delta;
step_size = param.step_size;

% ==== Initialize parameters
num_samples = ceil(b*sparsity*log2(n)); 
Z =(2*(rand(num_samples,n) > 0.5) - 1);  % Generate Rademacher sampling vecs
cosamp_params.Z = Z;
cosamp_params.sparsity = sparsity;
cosamp_params.delta = delta;
cosamp_params.n = 20; % Hardcode number of iterations of CoSaMP.

for k=1:maxit
    cosamp_params.x = x;
    try
        [f_k,grad_estimate,err] = CosampGradEstimate(f,fparam, cosamp_params);
    catch ME
        disp('An error has occurred. Terminating run of ZORO')
        converged = false;
        break
    end
    x = x - step_size*grad_estimate;
    objval_seq(k) = f_k;
    num_queries(k+1) = num_queries(k) + num_samples + 1;
    disp(['ZORO: Current function value is ', num2str(f_k), ' and total number of queries is ', num2str(num_queries(k+1))])
    if (num_queries(k+1)>param.budget)
        if requires_params
            objval_seq(k+1) = f(x, fparam); 
        else
            objval_seq(k+1) = f(x); 
        end
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
    if isfield(fparam, 'fmin') && early_stopping == true
        % eps = EPS_MORE * (f0 - fmin)
        if (f_k < fmin + eps) 
            if verbose
                disp(['ZORO Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(f_k)]);
            end
            converged = true;
            break;
        end
    end
end

% Package and return results
num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);

% package
Result.objval_seq = objval_seq;
Result.num_queries = num_queries;
Result.sol = x;
    
    
