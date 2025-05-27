% ===================================================================== %
% Benchmarking ZORO-FA against a variety of other algorithms.
% This script generates sample trajectories.
% Test functions: sparse quadratic, sparse quartic, max-s-squared, worst
% function in the world.
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

algorithms = {@DFQRM_B, @ZORO, @adaZORO, @ZORO_FA, @Nelder_Mead, @StochasticSubspaceDescent, @DirectSearchRR};

% ==== Parameters determining the run
n = 1000;
s = 30; %true sparsity
budget = 100; %NB: the number of fevals allowed is budget*(problem dim + 1)
lambda = 8; % for worst function only.
% following are only relevant for sparse quartic
% S = datasample(1:n,s,'Replace', false); % Sample s random indices in range 1:d
% fparam.S = S;
% B = rand(s);
% fparam.A = B'*B;

fparam.s = s;
fparam.n = n;
fparam.noise_mag = 0; % no noise for now.
fparam.lambda = lambda;
fparam.fmin = 0; % true minimum value for max-s-squared.
fparam.fmin = 1.01*(-lambda*s/(8*(s+1)));  %for worst function. Multiplying by factor as noted numerical error.
%temp_fun = @SparseQuadratic;
%temp_fun = @Max_s_squared;
temp_fun = @Worst_s_function;
%temp_fun = @SparseSkewQuartic;
fparam.requires_params = false;


fparam.f = @(x)temp_fun(x, fparam);

% ==== Common params
x0 = 10*randn(n,1);
fx0 = fparam.f(x0);
maxit = 1e6;

% ==== Define the parameters for ZORO
param.sparsity = s; %feed ZORO the true sparsity.
param.maxit = maxit;
param.delta = 0.0001;
param.step_size = 1/(lambda);
param.x0 = x0;
param.budget = (n+1)*budget;
param.n = n;
param.verbose = true;
param.num_samples = s; % dimension of the subspace for Stochastic Subspace.
%param.b = 1.0; % Compressed sensing parameter for ZORO-type algorithms.
param.early_stopping = false;

% ==== Plotting options
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
fmin_computed = 100;
for j=1:length(algorithms)
    if j == 4
        param.sparsity = ceil(0.02*n);
        param.epsilon = 0.01;
        param.sigma0 = 2.5;
        param.theta = 0.25;
    end
    if j == 5
        % Nelder-Mead is time consuming, so we don't want to run it each
        % time. Next few lines are a cheap way to check what hte objective
        % function is, and then load the relevant NM results.
        if fparam.fmin == 0
            load('max_s_squared_Nelder_Mead.mat');
        else
            load('worst_function_Nelder_Mead.mat');
        end
    else
        temp_Results = feval(algorithms{j}, fparam, param);
    end
    sl = mod(j-1,3) + 1;
    option1 = [char(lines(sl)) colors(j)];
    num_queries = temp_Results.num_queries;
    function_values = temp_Results.objval_seq;
    fmin_temp = min(function_values);
    hl(j) = semilogy(num_queries/(n+1), function_values - fparam.fmin, option1, 'LineWidth', 3);
    hold on
end

legend(labels)
axis([0 budget + 5 0 1.1*fx0])
set(gca, 'FontSize', 18)
set(gca, 'LineWidth', 1)