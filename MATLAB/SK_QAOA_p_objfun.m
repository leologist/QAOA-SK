function obj = SK_QAOA_p_objfun(gammas, betas)
%%SK_QAOA_p_objfun evaluates the objective function -E[<C>]/n for
%   averaged over random instances of the Sherrington-Kirkpatrick problem
%     i.e., C = \sum_{i<j} J_{ij} with J_{ij} ~ Normal(0, 1)
%
%   Usage:
%     First call SK_QAOA_p_objfun_helper, then given parameters
%     gammas = [gamma1, gamma2, ...], betas = [beta1, beta2, ...]
%     this function outputs:
%
%           obj = SK_QAOA_p_objfun(gammas, betas)
%
%     where obj = - E[ <C> ] / n

global p aList aHashes
global uniqueBarInds_sorted this_aHashes this_aBarHashes
global quickPrimer MrPlus MrMinus

VERBOSE = false;

start_time = tic;

%% Calculate list of X_u

Qs = Q_fun(aList, betas);
Fs = F_fun(aList, gammas, betas);

Xs = Qs.*Fs; % X(u) = Q(u)*F(u)

% Make a sorted Ys vector so that
%   Ys_by_hash(hash+1) = Y_fun(hash2config(hash), gammas, betas)

temp = sortrows([aHashes, Xs]); % sort by first column
Xs_by_hash = temp(:,2);
clear temp


%% Calculate list of R_u = Y_u

numRelevantConfigs = (2^(2*p)-2^p)/2;
numNoPlayWellConfigs = (4^p - 2^p - 2^(2*p-1))/2;

% only needs the rows of K_ij where i does not play well with something

% IdoesNotPlayWell = uniqueBarInds_sorted(any(PlaysWellTruthTable(mod(aHashes(uniqueBarInds_sorted), 2^p)+1, :), 2));
% Note that mod(hash, 2^p)+1 gives the index of the partition that a
% hashed configuration belongs to, due to the way hash index is calculated
% and PlaysWellTruthTable(A, B) is 1 <=> partition A plays well with B

IdoesNotPlayWell = uniqueBarInds_sorted(1:numNoPlayWellConfigs);

RList = ones(1, numRelevantConfigs);
for ind = 1:length(IdoesNotPlayWell)
    Klist = Xs(IdoesNotPlayWell(ind)) * ...
        Delta_fun(aList(IdoesNotPlayWell(ind), :), ...
                  aList(uniqueBarInds_sorted, :), ...
                  gammas);
    RList = RList.*exp(RList(ind) * Klist);
end

FullRList_by_hash = ones(2^(2*p),1);

FullRList_by_hash(this_aHashes + 1) = RList;
FullRList_by_hash(this_aBarHashes + 1) = RList;

%% Calculate the objective function

Ws_by_hash = Xs_by_hash .* FullRList_by_hash;

obj = 0;
for r = 1:p
    u = hash2config(MrPlus{r});
    ur = u(:, r);
    v = hash2config(MrMinus{r});
    vr = v(:, r);

    up_hash = quickPrimer(MrPlus{r});
    vp_hash = quickPrimer(MrMinus{r});

    % use the fact the 
    %   \sum_{u,v} u_r v_r W_{u'} W_{v'} = (\sum_u u_r W_{u'}) (\sum_v v_r W_{v'})

    obj = obj + 2 * gammas(r) * ...
        sum(ur .* Ws_by_hash(up_hash+1)) * ...
        sum(vr .* Ws_by_hash(vp_hash+1));
end

% note the added minus sign for minimization
%   (that corresponds to actual maximization of E[<C>]/n)
obj = -1i*obj; 

if VERBOSE
    fprintf('*** SK_QAOA_p_objfun done after %0.6f sec\n', toc(start_time));
end

% fprintf('obj=%0.16f\n', obj) 
% fprintf('gammas = [');
% fprintf('%0.16f, ', gammas);
% fprintf('\b\b];\n');
% fprintf('betas = [');
% fprintf('%0.16f, ', betas);
% fprintf('\b\b];\n');

if imag(obj) < 1e-10
    obj = real(obj);
else
    error('objective function has large imaginary part!')
end
