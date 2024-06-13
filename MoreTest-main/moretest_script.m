%---------------------------
% Run the same experiment again?
again = false;
%---------------------------

ProbType = 'MOR'; moretest1;

tic;
Results = cell(0,1);

etime2=zeros(0,1);

i = 0;
algname = cell(0,1);

param.eps_dep_mu = false;
param.adaptive_Lhat = false;
param.coord_change = false;
param.budget = budget;

%% CARS
i = i+1;
algname{i} = 'CARS';
if verbose>1
    disp(['Starting ', algname{i}, ' ...']);
end
tic;
Results{i} = CARS(fparam, param, 0);
etime2(i) = toc;
Results{i}.name = algname{i};



%% STP - variable step size
% i = i+1;
% 
% algname{i} = 'STP-vs';
% if verbose>1
%     disp(['Starting ', algname{i}, ' ...']);
% end
% Results{i} = SMTP(fparam, param, algname{i}, false);
% Results{i}.name = algname{i};
% 
% SMTP
param.fixed_step = false;
i = i+1;
algname{i} = 'SMTP';
if verbose>1
    disp(['Starting ', algname{i}, ' ...']);
end
Results{i} = SMTP(fparam, param, algname{i}, true);
Results{i}.name = algname{i};


%% Nesterov-Spokoiny for comparison (RGF)
%param.num_samples = 10; % number of directions per iteration
i = i+1;
algname{i} = 'Nesterov';
if verbose>1
    disp(['Starting ', algname{i}, ' ...']);
end
tic;
Results{i} = NesterovRS(fparam,param);
etime2(i) = toc;
Results{i}.name = algname{i};

% additional parameters for ZORO-LS
param.sparsity = ceil(0.05*n); % Let's be consistent with this initial sparsity choice
param.epsilon = 0.001;
param.sigma0 = 1.;
param.theta = 0.25;

i = i+1;
algname{i} = 'ZORO-LS';
if verbose>1
    disp(['Starting ', algname{i}, ' ...']);
end
tic;
Results{i} = ZORO_LS(fparam,param);
etime2(i) = toc;
Results{i}.name = algname{i};

% i = i+1;
% algname{i} = 'ZORO-ada-LS';
% if verbose>1
%     disp(['Starting ', algname{i}, ' ...']);
% end
% tic;
% Results{i} = ZORO_ada_LS(fparam,param);
% etime2(i) = toc;
% Results{i}.name = algname{i};


%% SPSA methods
% i = i+1;
% algname{i} = 'SPSA';
% if verbose>1
%     disp(['Starting ', algname{i}, ' ...']);
% end
% tic;
% Results{i} = RealSPSA(fparam,param);
% etime2(i) = toc;
% Results{i}.name = algname{i};

% i = i+1;
% algname{i} = '2-SPSA';
% if verbose>1
%     disp(['Starting ', algname{i}, ' ...']);
% end
% tic;
% Results{i} = Real2SPSA(fparam,param);
% etime2(i) = toc;
% Results{i}.name = algname{i};

%% AdaDGS
i = i+1;
algname{i} = 'AdaDGS';
if verbose>1
    disp(['Starting ', algname{i}, ' ...']);
end
tic;
Num_Quad_Pts = 5;
Results{i} = AdaDGS(fparam, param, Num_Quad_Pts);
etime2(i) = toc;
Results{i}.name = algname{i};

% % CARS DGS
% i = i+1;
% algname{i} = 'CARS-DGS';
% if verbose>1
%     disp(['Starting ', algname{i}, ' ...']);
% end
% tic;
% Results{i} = CARS_DGS(fparam, param, 5);
% etime2(i) = toc;
% Results{i}.name = algname{i};
%% Display Results
if showplot
    errors = display_result(Results(1:(end)), param.x0, fparam, ProbType, n);
end
% errors = display_result(Results(1:(end)),true_sol, fparam, ProbType, n);


sp = 1;
% plot(sol_sq(1:sp:end,1), sol_sq(1:sp:end,2),'-.*');
% plot(sol_sq_n(1:sp:end,1), sol_sq_n(1:sp:end,2),'-.*');
% plot(x2(minj),y2(mini),'r*');
% legend('', Results{idx_compare1}.name, Results{idx_compare2}.name, 'True min');
% xlim([min(x), max(x)]); ylim([min(y), max(y)]);
%%
% figure(); hold on;
% plot(Results{idx_compare1}.num_queries, log10(Results{idx_compare1}.mu_seq));
% plot(Results{idx_compare2}.num_queries, log10(Results{idx_compare2}.mu_seq));
% legend(Results{idx_compare1}.name, Results{idx_compare2}.name, 'location', 'best');
% title('\mu sequence');
% figure();
% plot(Results{idx_compare1}.num_queries(2:end), log10(Results{idx_compare1}.ndelta)); hold on;
% plot(Results{idx_compare2}.num_queries(2:end), log10(Results{idx_compare2}.ndelta));
% legend(Results{idx_compare1}.name, Results{idx_compare2}.name, 'location', 'best');
% title('delta norms');