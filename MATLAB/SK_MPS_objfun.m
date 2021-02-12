function [obj_approx, fidSqX, fidSqW, XMPS, MPSout] = SK_MPS_objfun(p, param)

global d

d = 4;
gammas = param(1:p);
betas = param(p+1:end);

%%

jk_pairs = nchoosek(1:p, 2);
Idelete = diff(jk_pairs,[],2) == 1;
jk_pairs = [jk_pairs; fliplr(jk_pairs(~Idelete,:))];

Ipairs = sub2ind([p,p],jk_pairs(:,1), jk_pairs(:,2));

%% get X_{b'}

tstart = tic;

[XMPS, fidSqX] = getXbprime(gammas, betas, p);

fprintf('X_MPS obtained after %0.4f s\n', toc(tstart));

%% get low-rank basis

[myUs, ~, myLs] = getLowRankBasis(p);


%% initialize W = X, only in the U space

W_U_t = zeros(p, p);
for ind = 1:size(jk_pairs,1)
    W_U_t(jk_pairs(ind,1), jk_pairs(ind,2)) = MPSoverlap(XMPS, myUs(ind));
end

W0time = toc(tstart);
fprintf('W_0 obtained after %0.4f s\n', W0time);

%% applying the nonlinear map f
fidSqW = 1;
for ind = 1:p
    [Wout, fidSq, MPSout] = fmap(W_U_t, gammas, XMPS, myUs, p);
    fidSqW = fidSqW * fidSq;
    diffW = norm(W_U_t(Ipairs) - Wout.');
    W_U_t(Ipairs) = Wout;
    
    fprintf('W_{t=%d} obtained after %0.4f s, |Delta W|=%0.4e\n', ind, toc(tstart), diffW);
    if diffW < 1
        W_L = arrayfun(@(m) MPSoverlap(MPSout, m), myLs);
        fprintf('           ------ obj approx = %0.8f\n\n', 2i*(W_L(1:p).*W_L(p+1:end))*gammas(:));
    end
end

fprintf('-- avg fmap evaluation time = %0.4f s\n', (toc(tstart)-W0time)/p);

%% computing V_p

[W_L, fidSq, MPSout] = fmap(W_U_t, gammas, XMPS, myLs, p);

fidSqW = fidSqW * fidSq;


obj_approx = 2i*(W_L(1:p).*W_L(p+1:end))*gammas(:);

fprintf('V_p obtained after %0.4f s\n', toc(tstart));

