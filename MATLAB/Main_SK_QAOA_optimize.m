if ~exist('Q_fun', 'file')
    addpath([pwd, filesep, 'SK_utils'])
end

global p maxD

p = 8;

USE_MPS = 1;

if p <= 12
    load('SK_opts.mat','SK_inf');
    param0 = SK_inf(p).param;
else
    gammas0 = [0.25, linspace(0.4,0.66,p-1)];
    betas0 = linspace(0.6,0.1,p);
    param0 = [gammas0, betas0];
end

%%

fprintf('============== p = %d ===============\n', p);

SAVE_FILE_NAME = sprintf('SK_QAOA_p=%d_opt.mat', p);


if USE_MPS
    maxD = 4^(floor(p/2));
    [myUs, ~, myLs] = getLowRankBasis(p);
    myfun = @(param) SK_MPS_objfun4optim(p, param, myUs, myLs, 1);
else
    SK_QAOA_p_objfun_helper;
    myfun = @(param) SK_QAOA_p_objfun(param(1:p), param(p+1:end));
end

opt_start_time = tic;
myoptions = optimoptions('fminunc','GradObj','off','Display','iter',...
    'TolX',1e-6,'TolFun',1e-6, 'Algorithm', 'quasi-newton', 'FiniteDifferenceType','central',...
    'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'PlotFcns',{@optimplotfval, @optimplotx});
[param, fval, exitflag, output] = fminunc(myfun, param0, myoptions);

fprintf('===== Optimization Done for p=%d after %0.2f sec =========\n', p, toc(opt_start_time));

%%
save(SAVE_FILE_NAME, 'param', 'fval', 'exitflag', 'output', 'param0')
