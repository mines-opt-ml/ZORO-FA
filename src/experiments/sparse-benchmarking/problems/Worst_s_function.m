function [val,grad] = Worst_s_function(x_in,function_params)
%        val = Worst_s_function(x_in, function_params)
% Implements the variant of Nesterov's worst function in the world, as
% described in "A stochastic subspace method ..." by Kozak et al.
% =========================== INPUTS ================================= %
% x_in ...................... Point at which to evaluate
% function_params ........... Contains sparsity parameter (s) and a noise
%                             magnitude parameter (noise_mag).
%
% ========================== OUTPUTS ================================== %
% 
% val ...................... noisy function evaluation at x_in
% grad ..................... exact (ie no noise) gradient evaluation at
% x_in
%
% Daniel Mckenzie
% May 2025
%
 
% =========== Unpack function_params 
s = function_params.s;
D = function_params.n;
if isfield('lambda', function_params)
    lambda = function_params.lambda;
else
    lambda = 8;
end
if isfield('noise_mag', function_params)
    noise_mag = function_params.noise_mag;
    noise = noise_mag*randn(1)./sqrt(D);
else
    noise = 0;
end

val = (lambda/4)*(sum(x_in(1:s).^2) - sum(x_in(1:s-1).*x_in(2:s))) - (lambda/4)*x_in(1) + noise;
grad = zeros(D,1); % grad not yet implemented.

end

