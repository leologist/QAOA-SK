if ~exist('Q_fun', 'file')
    addpath([pwd, filesep, 'SK_utils'])
end

global p

p = 5;
gammas0 = [0.25, 0.5, 0.5, 0.6, 0.6];
betas0 = [0.6, 0.5, 0.4, 0.2, 0.1];

param0 = [gammas0, betas0];


fprintf('============== p = %d ===============\n', p);

SAVE_FILE_NAME = sprintf('SK_QAOA_p=%d_opt.mat', p);



SK_QAOA_p_objfun_helper;
myfun = @(param) SK_QAOA_p_objfun(param(1:p), param(p+1:end));


opt_start_time = tic;
myoptions = optimoptions('fminunc','GradObj','off','Display','iter',...
    'TolX',1e-4,'TolFun',1e-4, 'Algorithm', 'quasi-newton', 'FiniteDifferenceType','forward',...
    'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'PlotFcns',{@optimplotfval, @optimplotx});
[param, fval, exitflag, output] = fminunc(myfun, param0, myoptions);

fprintf('===== Optimization Done for p=%d after %0.2f sec =========\n', p, toc(opt_start_time));

%%
save(SAVE_FILE_NAME, 'param', 'fval', 'exitflag', 'output', 'param0')
