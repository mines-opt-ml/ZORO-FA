% ===================================================================== %
% Small test for ZORO with Line Search (ZORO_LS). Compare to Vanilla ZORO
% Geovani Luis Grapiglia and Daniel McKenzie
% March 2022
% ===================================================================== %

clear, close all, clc

%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath(genpath('./problems/'))

algorithms = {@adaZORO}; %{@DFQRM_B, @Nelder_Mead, @ZORO, @ZORO_FA, @adaZORO};

% ==== Parameters determining the run
n = 500;
s = 30; %true sparsity
budget = 100; %NB: the number of fevals allowed is budget*(problem dim + 1)
S = datasample(1:n,s,'Replace', false); % Sample s random indices in range 1:d
fparam.s = s;
fparam.S = S;
fparam.n = n;
B = rand(s);
fparam.A = B'*B;
fparam.noise_mag = 0; % no noise for now.
fparam.fmin = 0; % true minimum value.
temp_fun = @SparseQuadratic;
%temp_fun = @Max_s_squared;
%temp_fun = @SparseSkewQuartic;
fparam.requires_params = false;
fparam.f = @(x)temp_fun(x, fparam);

% ==== Common params
x0 = 10*randn(n,1);
fx0 = fparam.f(x0);
num_iters = 1e6;

% ==== Define the parameters for ZORO
param.maxit = num_iters;
param.delta = 0.001;
param.step_size = 0.5;
param.x0 = x0;
param.budget = (n+1)*budget;
param.n = n;
param.verbose = true;

% ==== Plotting options
% NB: have switched 3rd and 4th options to accomodate ZORO. Remember to
% switch it back.
colors  = ['b' 'r' 'm' 'k' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' 'v' '^' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
labels{1} = 'adaZORO';
% labels{1} = 'DFQRM';
% labels{2} = 'Nelder-Mead';
% labels{3} = 'ZORO';
% labels{4} = 'ZORO-FA';

hl = zeros(4,1);
for j=1:length(algorithms)
    if j == 1
        param.sparsity = ceil(0.05*n);
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
axis([0 105 0 1.1*fx0])
set(gca, 'FontSize', 18)
set(gca, 'LineWidth', 1)

