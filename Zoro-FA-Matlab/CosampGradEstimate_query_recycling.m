function [func_val,grad_estimate,sampling_data] = CosampGradEstimate_query_recycling(function_handle,fparam, cosamp_params)
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
tol = 1e-3; % tolerance for CoSaMP. Do we still want this?
x = cosamp_params.x;
Z_new = cosamp_params.Z_new;
sparsity = cosamp_params.sparsity;
delta = cosamp_params.delta;
n = cosamp_params.n;
num_samples = cosamp_params.num_samples;
sizes = size(Z_new);
num_new_samples = sizes(1);
y_new = zeros(num_new_samples,1);

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
disp(['function value is ', num2str(func_val)])

for i = 1:num_new_samples
    if requires_params
        query_point = function_handle(x + delta*Z_new(i,:)', fparam);
    else
        query_point = function_handle(x + delta*Z_new(i,:)');
    end
    y_new(i) = (query_point - func_val)/(delta);
    disp(['query point is ', num2str(query_point)])
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

disp(['several entries of y ', num2str(y(7)), ' and ' num2str(y(50))])
if sparsity < length(x)/2 -1 
    %disp(['sparsity is ', num2str(sparsity)])
    [grad_estimate,err] = cosamp(Z/sqrt(num_samples),y/sqrt(num_samples),sparsity,tol,n);
else
    grad_estimate = Z\y;
    err = 0; %FIX THIS LATER.
end
sampling_data.Z = Z;
sampling_data.y = y;
sampling_data.err = err;

end