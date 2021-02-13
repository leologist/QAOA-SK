function [obj_approx, fidSqX, fidSqW]  = SK_MPS_objfun4optim(p, param, myUs, myLs, verboseLevel)
%SK_MPS_objfun4optim
%
%   [obj_approx, fidSqX, fidSqW] ...
%        = SK_MPS_objfun4optim(p, param, myUs, myLs, verboseLevel)

global d

d = 4;
gammas = param(1:p);
betas = param(p+1:end);

%%

if nargin <= 4
    verboseLevel = 1;
end

jk_pairs = nchoosek(1:p, 2);
Idelete = diff(jk_pairs,[],2) == 1;
jk_pairs = [jk_pairs; fliplr(jk_pairs(~Idelete,:))];

Ipairs = sub2ind([p,p],jk_pairs(:,1), jk_pairs(:,2));

%% get X_{b'}

tstart = tic;

[XMPS, fidSqX] = getXbprime(gammas, betas, p);

logMessage(2, 'X_MPS obtained after %0.4f s\n', toc(tstart));

%% initialize W = X, only in the U space

W_U_t = zeros(p, p);
for ind = 1:size(jk_pairs,1)
    W_U_t(jk_pairs(ind,1), jk_pairs(ind,2)) = MPSoverlap(XMPS, myUs(ind));
end

W0time = toc(tstart);
logMessage(2, 'W_0 obtained after %0.4f s\n', W0time);

%% applying the nonlinear map f
fidSqW = 1;
for ind = 1:p-1
    if ind == p-1
        [W_U_t(Ipairs), fidSq, MPSout] = fmap(W_U_t, gammas, XMPS, myUs, p);
    else
        [W_U_t(Ipairs), fidSq] = fmap(W_U_t, gammas, XMPS, myUs, p);
    end
    fidSqW = fidSqW * fidSq;
    
    logMessage(2, 'W_{t=%d} obtained after %0.4f s\n', ind, toc(tstart));
end

logMessage(1, 'avg fmap evaluation time = %0.4f s\n', (toc(tstart)-W0time)/p);

%% computing V_p

% [W_L, fidSq] = fmap(W_U_t, gammas, XMPS, myLs, p);
% fidSqW = fidSqW * fidSq;

W_L = arrayfun(@(mps) MPSoverlap(MPSout, mps), myLs);

obj_approx = 2i*(W_L(1:p).*W_L(p+1:end))*gammas(:);

if abs(imag(obj_approx)) > 1e-7
    error('obj has large imag part? %0.4e', imag(obj_approx));
end

obj_approx = -real(obj_approx);

% fprintf('gammas = [');
% fprintf('%0.6f, ', gammas);
% fprintf('\b\b];\n');
% fprintf('betas = [');
% fprintf('%0.6f, ', betas);
% fprintf('\b\b];\n');

logMessage(1, '** V_p=%0.12f obtained after %0.4f s\n', -obj_approx, toc(tstart));

    function logMessage(lvl, varargin)
        if lvl <= verboseLevel
            fprintf(varargin{:});
        end
    end
end
