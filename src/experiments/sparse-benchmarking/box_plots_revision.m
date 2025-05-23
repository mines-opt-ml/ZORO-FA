% ===================================================================== %
% Benchmarking ZORO-FA against a variety of other algorithms
% This script generates box plots
% Test function: sparse quadratic, sparse quartic, and max-s-squared
% functions
% Geovani Luis Grapiglia and Daniel McKenzie.
% March 2022--December 2024
% May 2025: Modified to include results for Direct Search and Stochastic
% Subspace Descent.
% ===================================================================== %

clear, close all, clc

%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath(genpath('./problems/'))

algorithms = {@StochasticSubspaceDescent, @DirectSearchRR};
num_algorithms = length(algorithms);
num_trials = 10;
final_values = zeros(num_trials, num_algorithms);
times = zeros(num_trials, num_algorithms);

% ==== Parameters determining the run
n = 1000;
s = 30; %true sparsity
budget = 500; %NB: the number of fevals allowed is budget*(problem dim + 1)
S = datasample(1:n,s,'Replace', false); % Sample s random indices in range 1:d
fparam.s = s;
fparam.S = S;
fparam.n = n;
fparam.noise_mag = 0; % no noise for now.
fparam.fmin = 0; % true minimum value.
%temp_fun = @SparseQuadratic;
temp_fun = @Max_s_squared;
%temp_fun = @SparseSkewQuartic;
num_iters = 1e6;

% ==== Define the parameters for ZORO
param.sparsity = s; %feed ZORO the true sparsity.
param.maxit = num_iters;
param.delta = 0.001;
param.step_size = 0.5;
param.budget = (n+1)*budget;
param.n = n;
param.verbose = true;

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

for k = 1:num_trials
    % random parameters that determine the run.
    A = rand(s);
    fparam.B = A'*A;
    x0 = 10*randn(n,1);
    param.x0 = x0;
    fparam.requires_params = false;
    fparam.f = @(x)temp_fun(x, fparam);
    fx0 = fparam.f(x0);
    %%%%
    for j=1:num_algorithms
        if j == 4
            param.sparsity = ceil(0.01*n);
            param.epsilon = 0.01;
            param.sigma0 = 1;
            param.theta = 0.25;
        end
        tic
        temp_Results = feval(algorithms{j}, fparam, param);
        times(k,j) = toc;
        param.sparsity = s; %reset sparsity for ZORO and adaZORO.
        final_values(k,j) = temp_Results.objval_seq(end);
    end
end

figure;
boxplot(final_values);
set(gca, 'XTickLabel', labels)
set(gca, 'FontSize', 18)
set(gca, 'YScale', 'log')
ylabel('$f(x_K) - f_{\star}$ (log scale)', 'Interpreter', 'latex');

figure;
boxplot(times);
set(gca, 'XTickLabel', labels)
set(gca, 'FontSize', 18)
set(gca, 'YScale', 'log')
ylabel('Run time in seconds (log scale)', 'Interpreter', 'latex');
