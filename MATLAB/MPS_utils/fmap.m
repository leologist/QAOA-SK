function [Wout, fidSq, XMPS] = fmap(W_U, gammas, XMPS, LMPSes, p)
%fmap implements the non-linear map f as a tensor network contraction
%
%   The map f is defined in the configuration basis as
%      [f(W)]_b = exp(\sum_a W_a \Delta_{a,b}) X_b
%
%   Usage:
%       [Wout, fidSq] = fmap(W_U, gammas, XMPS, LMPSes, p)
%
%       0 <= fidSq <= 1 is the estimated relative squared fidelity remaining 
%       truncating Schmidt values

global d
assert(d == 4);

Uswap = swap_gate(d);

fidSq = 1;

current_order = 1:p;
for round = 1:p
    
    start_ind = mod(round-1,2) + 1; % alternate 1, 2, 1, 2 between rounds
    
    for jj = start_ind:2:p-1
        my_j = current_order(jj);
        my_k = current_order(jj+1);
        
        [j, k] = fix_order(my_j, my_k);

        if k-j == 1
            T_jk_op = Uswap  ...
                * T_jk_pm(my_j, my_k, 1, W_U(j,k), gammas);
        else
            T_jk_op = Uswap  ...
                * T_jk_pm(my_j, my_k, -1, W_U(k,j), gammas) ...
                * T_jk_pm(my_j, my_k, 1, W_U(j,k), gammas);
        end

        [XMPS, err] = MPSTwoSiteOp(XMPS, T_jk_op, jj, false);

        fidSq = fidSq *(1 -  err);

        % swap the qudits
        current_order(jj) = my_k;
        current_order(jj+1) = my_j;

    end
end

% fprintf('DEBUG: final order: ');
% fprintf('%d, ', current_order);
% fprintf('\b\b\n');

XMPS = MPSreverse(XMPS);

Wout = zeros(size(LMPSes));

for ind = 1:numel(LMPSes)
    Wout(ind) = MPSoverlap(XMPS, LMPSes(ind));
end


end

function [j, k] = fix_order(j,k)
    if j > k
        l = j; j = k; k = l;
    end
end

function op = T_jk_pm(j, k, sign, W_in, gammas)
    b_p = [1,1,-1,-1];
    b_m = [1,-1,1,-1];

    if j < k
        op = exp(- W_in * gammas(j) * gammas(k) * kron(b_p + sign * b_m, b_p - b_m) );
    else % j > k
        op = exp(- W_in * gammas(j) * gammas(k) * kron(b_p - b_m, b_p + sign * b_m) );
    end
    op = diag(op);
end