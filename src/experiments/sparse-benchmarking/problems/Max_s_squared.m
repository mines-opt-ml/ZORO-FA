function [val,grad] = Max_s_squared(x_in,function_params)
%        val = Max_s_squared(x_in, function_params)
% Implements the Max-s-squared function, as described in ... . Allows for
% noisy evaluations.
%
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
% March 2022
%
 
% =========== Unpack function_params 
noise_mag = function_params.noise_mag;
s = function_params.s;
D = function_params.n;

noise = noise_mag*randn(1)./sqrt(D);
[b, I] = maxk(x_in.^2, s);
val = sum(b) + noise;
grad = zeros(D,1);
grad(I) = 2*x_in(I);

end

