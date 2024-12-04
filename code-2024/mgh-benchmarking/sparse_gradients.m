%%%
% Testing gradients of functions in the MGH test set for sparsity.
% Daniel McKenzie and Geovani Luis Grapiglia, August 2022.
% Updated May 2024
%%%
clear, clc, close all
% Some test hyperparams
addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath("problems/")

fnames = {'band', 'bv', 'ie', 'lin', 'lin0', 'lin1', 'pen1', 'rosex', 'singx', 'trid', 'trig', 'vardim'};  
num_func = length(fnames);
num_trials = 1000;
maxit=1e5; % not relevant for this experiment, but needed by function_builder
budget=100; % not relevant for this experiment, but needed by function_builder

for i = 1:num_func
    n = 500; 
    s = 0; % use either 0 or 1 to toggle the initial point
    fname = fnames(i);
    function_builder;
    Sorted_grads = zeros(num_trials, n); 

    for j=1:num_trials
        scale = randsample([0.0001, 0.001, 0.01, 0.1, 1, 10],1);
        x = scale*randn(n,1);
        [f_vec, J] = f_M(x);
        grad = J'*f_vec; % gradient of nonlinear least squares.
        sorted_grad = flip(sort(abs(grad)));
        Sorted_grads(j,:) = sorted_grad;
    end
    figure
    grad_mean = mean(Sorted_grads);
    grad_min = min(Sorted_grads);
    grad_max = max(Sorted_grads);
    
    alpha(0.5)
    plot(grad_mean, 'LineWidth', 2)
    hold on
    x = 1:n;
    x2 = [x, fliplr(x)];
    inBetween = [grad_max, fliplr(grad_min)];
    h = fill(x2, inBetween, 'b');
    set(h, 'facealpha', 0.5);
    set(gca, 'FontSize',18)
    %title(fname)
end