% test prima

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath("problems/")

problems = {'band', 'bl', 'bv', 'ie', 'lin', 'lin0', 'lin1', 'pen1', 'pen2', 'rosex', 'singx', 'trid', 'trig', 'vardim'};

% Parameters determining the run
maxit=1e5; %so large it will never be reached.
budget=5000; %NB: the number of fevals allowed is budget*(problem dim + 1)
n = 500; % use same dimension for all problems
s = 0; % use either 0 or 1 to toggle the initial point
tolerance = 1e-3;

fname = problems(1);
function_builder; % script that creates param, fparam

% reshape params for prima
prima_params.objective = fparam.f;
prima_params.x0 = param.x0;
prima_opts.solver = 'newuoa';
prima_opts.iprint = 1;
prima_opts.fortran = false;
prima_opts.maxfun = budget;



[x, fx, exitflag, output] = prima(fparam.f, param.x0, prima_opts);