% ====================================================================== %
% Benchmarking ZORO with line search (ZORO-LS) against various classical
% ZOO methods. This code builds upon that released with "Curvature-aware 
% derivative-free optimization" by Kim et al
% ====================================================================== %


%%%%%
% Add the directories containing ZORO-LS and benchmark algorithms
%%%%%

addpath(genpath('../Zoro-LS-Matlab'))
addpath(genpath('../Benchmark-algorithms'))

%%%%%%%%
% Parameters
%%%%%%%%%
ftns = [3,9,14,17,18,19,23,24,26,28,29,30,31];  % all probs with variable dimension
rep = 1;
eps = 1e-3;
budget = 1e5;
verbose = 1;
showplot = 0;
showplot_final = 1;
noise_lvl = 0; %eps*5e-2; %1e-4;
str_to_add = 'VRNS';
rM = 1e20;
xlims = [0, 12];
[EVALS, Results] = more_test(ftns, rep, eps, noise_lvl, budget, verbose, showplot);
save(['EVALS_',num2str(eps),'_',date,'_',str_to_add]);
disp(' ');
algnames = cell(length(Results),1);
for i=1:length(Results)
    algnames{i} = Results{i}.name;
end
[tau, rho, r] = performance_profile(algnames, EVALS, budget, rM, xlims, eps, showplot_final);
save(['perf_prof_',num2str(eps),'_',num2str(noise_lvl),'_',date],'tau','rho');
% disp('Evaluation done. Plotting');
% figure(); hold on;
% Opts_to_disp = 1:length(Results);
% for s= Opts_to_disp
%     linewidth = 2;
%     if s==8
%         line_spec = '-.';
%     else
%         line_spec = '-'; % default;
%     end
%     if (s==3) % only provides the total number of queries
%         plot(1:2:Results{s}.num_queries, log10(Results{s}.objval_seq), line_spec, 'LineWidth', linewidth);
%     else 
%         plot(Results{s}.num_queries, log10(Results{s}.objval_seq), line_spec, 'LineWidth', linewidth);
%     end
% end
% fontsize=16;
% legend(algnames{Opts_to_disp}, 'Location', 'best','Orientation','vertical');
% xlim([0,budget]);
% xlabel('Function Queries', 'FontSize', fontsize);
% ylabel('log_10(f(x))', 'FontSize', fontsize);