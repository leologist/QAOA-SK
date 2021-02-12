global maxD
maxD = 1000;

p = 11;

        % BENCHMARK @ maxD = 500
        % p   time    error     fmap time   
record = [9,  3.1207, -1.49e-5, 0.28;...
         10, 43.7388, -5.29e-4, 3.80;...
         11, 125.064,  2.72e-3, 9.89;...
         12, 504.538, -1.14e-2, 37.6];

%   Compare to GPU evluation
%   On Harvard cluster: p=11 takes <= 6.5 hours,  6875 MB memory
%                       p=10 takes <= 44 minutes, 3000 MB memory

fprintf('==== V_p @ p = %d ===== MPS maxD=%d =====\n', p, maxD);
if p <= 12
    gammas = SK_inf(p).param(1:p);
    betas =  SK_inf(p).param(p+1:end);
    exact_val = -SK_inf(p).fval;
else
    [gammas, betas] = get_ParamGuess(p, SK_inf);
    exact_val = nan;
end

[obj_approx, fidSqX, fidSqW] = SK_MPS_objfun(p, [gammas, betas]);

%% printing information

fprintf('** Comparing obj: approx = %0.8f, exact = %0.8f\n', obj_approx, exact_val);
fprintf('**** absolute error = %0.3e\n', obj_approx - exact_val);
fprintf('**** relative error = %0.3e\n', (obj_approx - exact_val)/exact_val);
fprintf('**** a prediction of relative error = %0.3e\n', 1-fidSqX^(1+p+1 + (p-1)^2)*fidSqW);