function [func_val, grad_estimate] = FDGradEstimate(function_handle, fparam, fd_params)
% Simple zeroth order gradient estimator using forward finite difference
% =============================== Inputs =============================== %
% function_handle ....... Handle of zeroth order oracle. Must be preceded by an "@"
% function_params ....... Additional params for function evaluation.
% cosamp_params ......... Holds the length-scale for finite differencing, x
%
% ============================== Outputs =============================== %
% func_val ............. returns f(x), as given by zeroth order oracle.
% grad_estimate ........ gradient estimator.
%
% Daniel Mckenzie
% August 2024.

delta = fd_params.delta;
x = fd_params.x;
n = length(x);
grad_estimate = zeros(n,1);
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

% ==== Compute gradient estimate.
for i=1:n
    e_i = zeros(n,1);
    e_i(i) = 1;
    if requires_params
        query_point = function_handle(x + delta*e_i, fparam);
    else
        query_point = function_handle(x + delta*e_i);
    end
    grad_estimate(i) = (query_point - func_val)/delta;
end


