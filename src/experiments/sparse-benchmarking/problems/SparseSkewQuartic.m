function [val,grad] = SparseSkewQuartic(x_in,function_params)
% Provides noisy evaluations of a sparse deg 
%
% =========================== INPUTS ================================= %
% x_in ...................... Point at which to evaluate
% S ......................... Suppose set of sparse quadric. Keep this the
% same
% D ......................... Ambient dimension
% sigma ..................... sigma/sqrt(D) is per component Gaussian noise level
%
% ========================== OUTPUTS ================================== %
% 
% val ...................... noisy function evaluation at x_in
% grad ..................... exact (ie no noise) gradient evaluation at
% x_in
%
% Daniel Mckenzie
% 26th June 2019
% Modified by Yuchen Lou
% August 2020
% Modified again by Daniel McKenzie
% March 2022
%

% =========== Check if x_in has been transposed
%This happens with Nelder-Mead?
size_x_in = size(x_in);
if size_x_in(1) == 1
    x_in = x_in';
end

% =========== Unpack function_params 
noise_mag = function_params.noise_mag;
S = function_params.S;
D = length(x_in);
s = length(S);
A = triu(ones(s))/s;
Ax_S = A*x_in(S);

noise = noise_mag*randn(1)./sqrt(D);
val = Ax_S'*Ax_S + 0.1*sum(Ax_S.^3) + 0.01*sum(Ax_S.^4) + noise;
grad = zeros(D,1);
grad(S) = 2*A'*Ax_S + 0.3*A'*Ax_S.^2 + 0.04*A'*Ax_S.^3;

end

