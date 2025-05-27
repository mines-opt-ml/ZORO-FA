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
% Daniel McKenzie
% March 2022

algname = 'adaZORO';

%% INITIALIZATION
initialization;

% Algorithm specific parameters.
sparsity = param.sparsity;
delta = param.delta;
step_size = param.step_size;
phi = 2.5e-1; %tolerance parameter. Hardcoding this for now.
if isfield(param, 'b')
    b = param.b;
else
    b = 0.5; % Use different default to ZORO and ZORO-FA as this seems to work better.
end

for k=1:maxit
    flag = false;
    j = 0;
    s_j = sparsity;
    num_queries_this_iter = 0;
    while flag == false
        num_samples = min(n,ceil(b*s_j*log2(n)));
        Z =(2*(rand(num_samples,n) > 0.5) - 1);  % Generate Rademacher sampling vecs
        cosamp_params.Z = Z;
        cosamp_params.sparsity = s_j;
        cosamp_params.delta = delta;
        cosamp_params.n = 20; % Hardcode number of iterations of CoSaMP.
        cosamp_params.x = x;
        try
            [f_k,grad_estimate,~, ~,err] = CosampGradEstimate(f,fparam, cosamp_params);
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
        disp(['Num queries thus far ', num2str(num_queries(k) +num_queries_this_iter)])

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
        if err < phi
            x = x - step_size*grad_estimate;
            objval_seq(k+1) = f_k;
            num_queries(k+1) = num_queries(k) + num_queries_this_iter;
            disp(['adaZORO: Obj val is ', num2str(f_k), ' and num samples is ', num2str(num_queries(k+1))])
            flag = true;
        end
        j = j+1;
        s_j = s_j + 10;

        % if isfield(fparam, 'fmin') && early_stopping == true
        % % eps = EPS_MORE * (f0 - fmin)
        % if (f_k < fmin + eps) 
        %     if verbose
        %         disp(['adaZORO Converged in ', num2str(k),' steps. Exit the loop']);
        %         disp(['Function val = ' , num2str(f_k)]);
        %     end
        %     converged = true;
        %     break;
        % end
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
    
    
