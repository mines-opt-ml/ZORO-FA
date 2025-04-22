function [func_val,grad_estimate,min_query_point, perturbation_yielding_min_query, err] = CosampGradEstimate(function_handle, fparam, cosamp_params)
%        grad_estimate = CosampGradEstimate(num_samples,delta)2(function_handle,x,num_samples,delta)
% Zeroth order gradient estimator using CoSaMP.
% =============================== Inputs =============================== %
% function_handle ....... Handle of zeroth order oracle. Must be preceded by an "@"
% function_params ....... Additional params for function evaluation.
% cosamp_params ......... Params required for CoSaMP. Includes x, Z,
%                         sparsity, delta, n 
%
% ============================== Outputs =============================== %
% func_val ............. returns f(x), as given by zeroth order oracle.
% grad_estimate ........ gradient estimator.
% 
%
% Daniel Mckenzie
% 2nd March 2022
% Modified April 2025 to find the smallest value of f(x + hZ[i,:])
% encountered, return this value, as well as the perturbation hZ[i,:].
%

% === Initialization and unpack various parameters
tol = 1e-3; % tolerance for CoSaMP. Do we still want this?
x = cosamp_params.x;
Z = cosamp_params.Z;
sparsity = cosamp_params.sparsity;
delta = cosamp_params.delta;
n = cosamp_params.n;
sizes = size(Z);
num_samples = sizes(1);
y = zeros(num_samples,1);
if isfield(fparam, 'requires_params')
    requires_params = fparam.requires_params;
else
    requires_params = false;
end


% ==== Evaluate and store the function value at x
if requires_params
    func_val = function_handle(x, fparam);
else
    func_val = function_handle(x);
end

query_points = zeros(num_samples,1); % vector for storing all function evaluations used.

for i = 1:num_samples
    if requires_params
        query_point = function_handle(x + delta*Z(i,:)', fparam);
    else
        query_point = function_handle(x + delta*Z(i,:)');
    end
    query_points(i) = query_point;
    y(i) = (query_point - func_val)/(delta*sqrt(num_samples));
end

% find best query point examined
[min_query_point, idx_of_min_query_point] = min(query_points);
perturbation_yielding_min_query = delta*Z(idx_of_min_query_point,:)';
%if num_samples < length(x) & sparsity < length(x)/2 - 1
[grad_estimate, err] = cosamp(Z,y,sparsity,tol,n);
% else
%     grad_estimate = Z\y;
%     err = 0;
    %disp(['Dimension is ', num2str(sizes(2)), ' and size of Z is ', num2str(num_samples)])
    %disp('Using least squares gradient estimator instead.')
%end

end
