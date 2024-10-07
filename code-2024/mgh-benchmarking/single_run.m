%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../../Zoro-FA-Matlab'))
addpath(genpath('../../Benchmark-algorithms'))
addpath("problems/")

problems = {'band', 'bl', 'bv', 'ie', 'lin', 'lin0', 'lin1', 'pen1', 'pen2', 'rosex', 'singx', 'trid', 'trig', 'vardim'};
algorithms = {@DFQRM_B, @Nelder_Mead, @ZORO, @ZORO_FA, @adaZORO};%, 'prima'};
Results = cell(2*length(problems), length(algorithms));
function_evaluations = cell(2*length(problems), length(algorithms));
starting_value = zeros(2*length(problems),1);
ZORO_avg_sparsity = zeros(2*length(problems),1);
N = zeros(2*length(problems),1); % record size of the problems.

% Parameters determining the run
maxit=1e5; %so large it will never be reached.
budget= 100; %100; %NB: the number of fevals allowed is budget*(problem dim + 1)
n = 500; % use same dimension for all problems
%s = 0; % use either 0 or 1 to toggle the initial point
%tolerance = 1e-1;

for k=1:2*length(problems)
    if mod(k,2) == 0
        s = 1;
    else
        s = 0;
    end
    fname = problems(floor((k+1)/2));
    function_builder; % script that creates param, fparam
    starting_value(k) = fx0;
    % ==== additional parameters for ZORO-FA
    param.sparsity = ceil(0.05*n); % Let's be consistent with this initial sparsity choice
    param.epsilon = 0.01;
    param.sigma0 = 1.;
    param.theta = 0.25;
    N(k) = n;
    % ==== additional parameters for ZORO
    param.step_size = 0.001;
    param.verbose = true;
    param.delta = 0.005;

    for j=1:length(algorithms)
        if strcmp(algorithms{j}, 'prima')
            % reshape inputs for prima
            prima_params.objective = fparam.f;
            prima_params.x0 = param.x0;
            prima_opts.solver = 'newuoa';
            prima_opts.iprint = 1;
            prima_opts.fortran = false;
            prima_opts.maxfun = param.budget;

            % run prima
            [x, fx, exitflag, output] = prima(fparam.f, param.x0, prima_opts);

            % reshape output
            temp_Results.algname = 'newuoa';
            temp_Results.objval_seq = output.fhist;
            temp_Results.num_queries = 1:length(output.fhist);
        elseif j== 1 %j==3 || j == 5  % ZORO or adaZORO
            params.sparsity = ceil(0.1*n);
            temp_Results = feval(algorithms{j}, fparam, param);
            params.sparsity = ceil(0.05*n); % set sparsity back to default value.
        else
            temp_Results = feval(algorithms{j}, fparam, param);
            % min_feval_achieved = min(temp_Results.objval_seq);
        end
        Results{k, j} = temp_Results;
        % function_evaluations{k, j} = temp_Results.num_queries;
        if isfield(temp_Results, 'sparse_seq')
           ZORO_avg_sparsity(k) = mean(temp_Results.sparse_seq);
        end
    end
end

tolerances = [1e-1, 1e-2, 1e-3];
for i=1:3
    tolerance = tolerances(i)
    % Generate data profiles
    h1 = data_profile(Results,N,starting_value, tolerance);
    set(gca, 'FontSize',18)
    set(gca, 'LineWidth', 1)
    figure
end




