function [func_val,grad_estimate,sampling_data] = CosampGradEstimate(function_handle, cosamp_params)
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
% Daniel Mckenzie
% 2nd March 2022
%

% === Initialization and unpack various parameters
tol = 1e-8; % tolerance for CoSaMP. Do we still want this?
x = cosamp_params.x;
Z_new = cosamp_params.Z_new;
sparsity = cosamp_params.sparsity;
delta = cosamp_params.delta;
n = cosamp_params.n;
num_samples = cosamp_params.num_samples;
sizes = size(Z_new);
num_new_samples = sizes(1);
y_new = zeros(num_new_samples,1);


% ==== Evaluate and store the function value at x
func_val = function_handle(x);

for i = 1:num_new_samples
    query_point = function_handle(x + (delta/sqrt(num_samples))*Z_new(i,:)');
    y_new(i) = (query_point - func_val)/(delta);
end

if isfield(cosamp_params, 'Z')
    Z = [cosamp_params.Z; Z_new];
else
    Z = Z_new;
end

if isfield(cosamp_params, 'y')
    y = [cosamp_params.y; y_new];
else
    y = y_new;
end

grad_estimate = cosamp(Z/sqrt(num_samples),y/sqrt(num_samples),sparsity,tol,n);
sampling_data.Z = Z;
sampling_data.y = y;

end