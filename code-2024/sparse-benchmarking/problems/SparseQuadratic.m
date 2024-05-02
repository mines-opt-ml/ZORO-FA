function [val,grad] = SparseQuadratic(x_in,function_params)
%        val = SparseQuadric(x_in)
% Provides noisy evaluations of a sparse quadric of the form x^T_{S}Ax^T
% where A is a symmetric s-by-s matrix and s = |S|.
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
% Modified to be more ill-conditioned by Daniel McKenzie (March 2022)
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
D = function_params.n;
A = function_params.A;

noise = noise_mag*randn(1)./sqrt(D);
val = x_in(S)'*A*x_in(S) + noise;
grad = zeros(D,1);
grad(S) = 2*A*x_in(S);

end

