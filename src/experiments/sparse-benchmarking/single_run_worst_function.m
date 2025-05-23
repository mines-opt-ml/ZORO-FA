% ===================================================================== %
% Benchmarking ZORO-FA against a variety of other algorithms.
% This script generates sample trajectories.
% Test function: Sparse version of Nesterov's worst function in the world.
% Geovani Luis Grapiglia and Daniel McKenzie.
% March 2022--May 2025
% ===================================================================== %

clear, close all, clc

%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath(genpath('./problems/'))

algorithms = {@DFQRM_B, @ZORO, @adaZORO, @ZORO_FA, @Nelder_Mead, @StochasticSubspaceDescent, @DirectSearchRR};

% ==== Parameters determining the run
n = 1000;
s = 30; %true sparsity
budget = 500; %NB: the number of fevals allowed is budget*(problem dim + 1)
fparam.s = s;
fparam.n = n;
fparam.noise_mag = 0; % no noise for now.
temp_fun = @Worst_s_function;
fparam.requires_params = false;
fparam.f = @(x)temp_fun(x, fparam);
lambda = 100;
fparam.lambda = lambda;
fparam.fmin = -lambda*s/(8*(s+1)); % true minimum value.

% ==== Common params
x0 = 100*randn(n,1);
fx0 = fparam.f(x0);
maxit = 1e6;

% ==== Define the parameters for ZORO
param.sparsity = s; %feed ZORO the true sparsity.
param.maxit = maxit;
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
labels{6} = 'SSD'; %stochastic subspace descent
labels{7} = 'DS-RS'; %Direct Search in Reduced Subspaces

hl = zeros(length(algorithms),1);
for j=1:length(algorithms)
    if j == 4 % ZORO-FA
        param.sparsity = ceil(0.01*n);
        param.epsilon = 0.01;
        param.sigma0 = 1;
        param.theta = 0.25;
    end
    temp_Results = feval(algorithms{j}, fparam, param);
    if j == 5
        save('worst_function_Nelder_Mead.mat',"temp_Results")
    end
    sl = mod(j-1,3) + 1;
    option1 = [char(lines(sl)) colors(j)];
    num_queries = temp_Results.num_queries;
    function_values = temp_Results.objval_seq;
    hl(j) = semilogy(num_queries/(n+1), function_values - fparam.fmin, option1, 'LineWidth', 3);
    hold on
end

%axis([0 505 0 1.1*fx0])
set(gca, 'FontSize', 18)
set(gca, 'LineWidth', 1)
legend(labels, 'FontSize', 12)