function Result = StochasticSubspaceDescent(fparam, param)
%        Result = StochasticSubspaceDescent(fparam, param)
% An implementation of the stochastic subspace descent algorithm analyzed
% in "A stochastic subspace approach to gradient-free optimization in high
% dimensions" by Kozak et al (2021).

algname = 'StochasticSubspaceDescent';
%% INITIALIZATION
initialization;

% Algorithm specific parameters
if isfield(param, 'num_samples')
    num_samples = param.num_samples;
else
    num_samples = n;
end
if isfield(param, 'step_size')
    alpha = (num_samples/n)*param.step_size; % adjust for subsampling.
else
    alpha = 0.01;
end

if isfield(param, 'delta') %sampling radius
    mu = param.delta;
else
    mu = 10^-5;
end

for k=1:maxit
    num_queries(k+1) = num_queries(k);
    fx = objval_seq(k); % use saved value
    delta = zeros(length(x),1); % the gradient approximation
    % Sample num_samples orthogonal directions. We use the approach
    % suggested by Kozak et al on pg. 348.
    Z = randn(n, num_samples);
    [Q,~] = qr(Z);
    for i=1:num_samples
        u = sqrt(n/num_samples)*Q(:, i);
        if requires_params
            fp = f(x+mu*u, fparam);
        else
            fp = f(x+mu*u); 
        end
        num_queries(k+1) = num_queries(k+1) + 1;
        d = (fp - fx)/mu; % directional derivative
        delta = delta + d*u;
    end
    
    delta = -alpha*delta; % move to the next iterate
    
    % Update
    if requires_params
        fxnew = f(x+delta, fparam);
    else
        fxnew = f(x+delta);
    end
    num_queries(k+1) = num_queries(k+1) + 1; %single query here
    x = x + delta;
    objval_seq(k+1) = fxnew;

    disp(['SSD. Obj val is ', num2str(fxnew), ' and num queries is ', num2str(num_queries(k+1))])
    
    if (num_queries(k+1)>param.budget)
        if verbose
            disp('Max queries hit!')
        end
        num_iter = k;
        break;
    end
end


objval_seq = objval_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);

% put into a struct for output
Result.objval_seq = objval_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.sol = x;
end
