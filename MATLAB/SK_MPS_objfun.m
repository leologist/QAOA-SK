function [obj_approx, fidSqX, fidSqW, W_answers] = SK_MPS_objfun(p, param)

global d

d = 4;
gammas = param(1:p);
betas = param(p+1:end);

%%
tstart = tic;

jk_pairs = nchoosek(1:p, 2);
Idelete = diff(jk_pairs,[],2) == 1;
jk_pairs = [jk_pairs; fliplr(jk_pairs(~Idelete,:))];

Ipairs = sub2ind([p,p],jk_pairs(:,1), jk_pairs(:,2));

%% get X_{b'}

[XMPS, fidSqX] = getXbprime(gammas, betas, p);

fprintf('X_MPS obtained after %0.4f s\n', toc(tstart));

%% get low-rank basis

[myUs, ~, myLs] = getLowRankBasis(p);


%% initialize W = X, only in the U space

W_lamb_t = zeros(p, p);
for ind = 1:size(jk_pairs,1)
    W_lamb_t(jk_pairs(ind,1), jk_pairs(ind,2)) = MPSoverlap(XMPS, myUs(ind));
end

W0time = toc(tstart);
fprintf('W_0 obtained after %0.4f s\n', W0time);

%% applying the nonlinear map f
fidSqW = 1;
for ind = 1:p
    [Wout, fidSq] = fmap(W_lamb_t, gammas, XMPS, myUs, p);
    fidSqW = fidSqW * fidSq;
    W_lamb_t(Ipairs) = Wout;
    
    fprintf('W_{t=%d} obtained after %0.4f s\n', ind, toc(tstart));
end

fprintf('-- avg fmap evaluation time = %0.4f s\n', (toc(tstart)-W0time)/p);

%% computing V_p

[W_answers, fidSq] = fmap(W_lamb_t, gammas, XMPS, myLs, p);
fidSqW = fidSqW * fidSq;


obj_approx = 2i*(W_answers(1:p).*W_answers(p+1:end))*gammas(:);

fprintf('V_p obtained after %0.4f s\n', toc(tstart));

