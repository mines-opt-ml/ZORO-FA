function Result = DirectSearchRR(fparam, param)
%        Result = DirectSearchRR(fparam, param)
% An implementation of the probabilistic direct search algorithm of 
% "Direct Search Based on Probabilistic Descent in Reduced Subspaces" by 
% Roberts and Royer (2023).
% We implement the most performant variant: r=1 and Gaussian P_k

algname = 'DirectSearchRR';
%% INITIALIZATION
initialization;

% algorithm specific parameters
if isfield(param, 'num_samples_RR')
    num_samples = param.num_samples_RR;
else
    num_samples = 1;
end
if num_samples > 1
    disp('Sorry, we have not implemented that functionality')
    return
end
if isfield(param, 'alpha0')
    alpha = param.alpha0;
else
    alpha = 1.0; %value used in experiments of Roberts and Royer
end
if isfield(param, 'alpha_max')
    alpha_max = param.alpha_max;
else
    alpha_max = 1000; % default value used in experiments of Roberts and Royer
end
if isfield(param, 'gamma_inc')
    gamma_inc = param.gamma_inc;
else
    gamma_inc = 2.0; %value used in experiments of Royer and Roberts was 2
end
if isfield(param, 'gamma_dec')
    gamma_dec = param.gamma_dec;
else
    gamma_dec = 0.5; % value used in experiments of Roberts and Royer.
end
if isfield(param, 'c')
    c = param.c;
else
    c = 0.001; %could not find default value in R&R's paper.
end

for k=1:maxit
    flag = false;
    counter = 1;
    fx = objval_seq(k);
    if isfield(fparam, 'fmin')
        if (fx < fmin + eps) 
            if verbose
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fx)]);
            end
            converged = true;
            break;
        end
    end
    if (num_queries(k)>param.budget)
        if verbose
            disp('Max queries hit!')
        end
        converged = false;
        break;
    end
    if alpha < 1e-10
        disp('Sampling radius is too small. Terminating.')
        disp(['Function val = ' , num2str(fx)]);
        converged = false;
        break;
    end

    num_queries(k+1) = num_queries(k);
    while flag == false
        u = randn(n,1);
        x_plus = x + alpha*u;
        target = c*alpha^2*norm(u,2)^2/2;
        if requires_params
            f_plus = f(x_plus, fparam);
        else
            f_plus = f(x_plus);
        end
    
        if f_plus < fx - target
            x = x_plus;
            flag = true;
            num_queries(k+1) = num_queries(k+1) + 1;
            objval_seq(k+1) = f_plus;
            alpha = min(gamma_inc*alpha, alpha_max);
            disp(['DS.  Obj val is ', num2str(f_plus), ' num queries is + ', num2str(counter)])
        end

        if flag == false
            x_minus = x - alpha*u;
            if requires_params
                f_minus = f(x_minus, fparam);
            else
                f_minus = f(x_minus);
            end

            if f_minus < fx - target
                x = x_minus;
                flag = true;
                num_queries(k+1) = num_queries(k+1) + 2; %also count the query at x_plus
                objval_seq(k+1) = f_minus;
                alpha = min(gamma_inc*alpha, alpha_max);
                disp(['DS.  Obj val is ' num2str(f_minus), ' num queries is ', num2str(counter)])
            end
        end
        alpha = gamma_dec*alpha; % decrease sampling radius if unsuccessful.
        counter = counter + 1;
    end
end

num_iter = k;
objval_seq = objval_seq(1:num_iter);
num_queries = num_queries(1:num_iter);

% put into a struct for output
Result.objval_seq = objval_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.sol = x;

end
