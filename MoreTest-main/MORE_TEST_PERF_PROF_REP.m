REP = 10;

ftns = 1; %1:34;
rep = 1;
eps = 1e-9;
budget = 5e4;
verbose = 0;
showplot = 0;
showplot_final = 0;
noise_lvl = eps*10e-2; %1e-4;

str_to_add = 'Quartic_ALL_REP';
rM = 1e20;
xlims = [0, 12];

EVALS = cell(REP,1);
Results = cell(REP,1);
for iREP = 1:REP

[EVALS{iREP}, Results{iREP}] = more_test(ftns, rep, eps, noise_lvl, budget, verbose, showplot);
fprintf('.');
end
disp(' ');

save(['EVALS_',num2str(eps),'_',date,'_',str_to_add]);

disp(' ');
nOpts = length(Results{1});
algnames = cell(nOpts,1);
for i=1:nOpts
    algnames{i} = Results{1}{i}.name;
end

avg_nqueries = cell(nOpts,1);
avg_objval = cell(nOpts,1);

% for i=1:nOpts
%     if i~=3
%         avg_nqueries{i} = Results{1}{i}.nqueries;
%     else % STP
%         avg_nqueries{i} = 1:2:Results{1}{i}.nqueries;
%     end
%     avg_objval{i} = Results{1}{i}.objval_seq;
%     for iREP = 2:REP
%         if i~=3
%             avg_nqueires{i} = avg_nqueires{i} + Results{iREP}{i}.nqueries;
%         else % STP
%             avg_nqueires{i} = avg_nqueires{i} + 1:2:Results{iREP}{i}.nqueries;
%         end
%         avg_objval{i} = avg_objval{i} + Results{1}{i}.objval_seq;
%     end
%     avg_objval{i} = avg_objval{i}/REP;
%     avg_nqueries{i} = avg_nqueries{i}/REP;
% end
% avg_objval{i} = avg_objval{i}/REP;

% [tau, rho, r] = performance_profile(algnames, EVALS, budget, rM, xlims, eps, showplot_final);

% save(['perf_prof_',num2str(eps),'_',num2str(noise_lvl),'_',date],'tau','rho');

disp('Evaluation done. Plotting');
figure(); hold on;
Opts_to_disp = 1:length(Results);
for s= Opts_to_disp
    linewidth = 2;
    if s==1 || s==2 % CARS, CARS-NQ
        line_spec = '-.';
    elseif s ==3 || s ==4 % STP,SMTP
        line_spec = '-.';
    elseif s==6 || s==7
        line_spec = '--';
    else
        line_spec = '-';
    end
    if (s==3) % only provides the total number of queries
        plot(1:2:Results{s}.num_queries, log10(Results{s}.objval_seq), line_spec, 'LineWidth', linewidth);
    else 
        plot(Results{s}.num_queries, log10(Results{s}.objval_seq), line_spec, 'LineWidth', linewidth);
    end
end
fontsize=16;
legend(algnames{Opts_to_disp}, 'Location', 'best','Orientation','vertical');
xlim([0,budget]);
xlabel('Function Queries', 'FontSize', fontsize);
ylabel('$log_{10}(f(x))$', 'FontSize', fontsize, 'Interpreter','latex');