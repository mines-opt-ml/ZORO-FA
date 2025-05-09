% ===================================================================== %
% Benchmarking ZORO-FA against a variety of other algorithms.
% This script generates sample trajectories.
% Test function: sparse quadratic, sparse quartic, and max-s-squared
% functions
% Geovani Luis Grapiglia and Daniel McKenzie.
% March 2022--December 2024
% ===================================================================== %

clear, close all, clc

%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath(genpath('./problems/'))

algorithms = {@DFQRM_B, @ZORO, @adaZORO, @ZORO_FA, @Nelder_Mead, @StochasticSubspaceDescent};

% ==== Parameters determining the run
n = 1000;
s = 30; %true sparsity
budget = 500; %NB: the number of fevals allowed is budget*(problem dim + 1)
S = datasample(1:n,s,'Replace', false); % Sample s random indices in range 1:d
fparam.s = s;
fparam.S = S;
fparam.n = n;
B = rand(s);
fparam.A = B'*B;
fparam.noise_mag = 0; % no noise for now.
fparam.fmin = 0; % true minimum value.
%temp_fun = @SparseQuadratic;
temp_fun = @Max_s_squared;
%temp_fun = @SparseSkewQuartic;
fparam.requires_params = false;
fparam.f = @(x)temp_fun(x, fparam);

% ==== Common params
x0 = 10*randn(n,1);
fx0 = fparam.f(x0);
num_iters = 1e6;

% ==== Define the parameters for ZORO
param.sparsity = s; %feed ZORO the true sparsity.
param.maxit = num_iters;
param.delta = 0.001;
param.step_size = 0.5;
param.x0 = x0;
param.budget = (n+1)*budget;
param.n = n;
param.verbose = true;
param.num_samples = 30; % dimension of the subspace for Stochastic Subspace.

% ==== Plotting options
% NB: have switched 3rd and 4th options to accomodate ZORO. Remember to
% switch it back.
colors  = ['b' 'm' 'c' 'k' 'r' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' 'v' '^' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
labels{1} = 'DFQRM';
labels{2} = 'ZORO';
labels{3} = 'adaZORO';
labels{4} = 'ZORO-FA';
labels{5} = 'Nelder-Mead';
labels{6} = 'SSD'; %stochastinc subspace descent

hl = zeros(4,1);
for j=1:length(algorithms)
    if j == 4
        param.sparsity = ceil(0.01*n);
        param.epsilon = 0.01;
        param.sigma0 = 1;
        param.theta = 0.25;
    end
    temp_Results = feval(algorithms{j}, fparam, param);
    sl = mod(j-1,3) + 1;
    option1 = [char(lines(sl)) colors(j)];
    num_queries = temp_Results.num_queries;
    function_values = temp_Results.objval_seq;
    hl(j) = semilogy(num_queries/(n+1), function_values, option1, 'LineWidth', 3);
    hold on
end

legend(labels)
axis([0 505 0 1.1*fx0])
set(gca, 'FontSize', 18)
set(gca, 'LineWidth', 1)