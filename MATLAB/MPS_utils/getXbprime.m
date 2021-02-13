function [XMPS, fidSq] = getXbprime(gammas, betas, p)
%getXbprime returns X_{b'} (approximately) as an MPS
%
%   [XMPS, fidSq] = getXbprime(gammas, betas, p)
%
%       0 <= fidSq <= 1 is the estimated relative squared fidelity remaining 
%       truncating Schmidt values

global d

d = 4;


%% Initialize X_b' = Q_b'

QMPS = getQbprime_as_MPS(p, betas);
XMPS = MPScanonicalize(QMPS);


%%
bs = [1     1;...
      1    -1;...
     -1     1;...
     -1    -1];

b_p = bs(:, 1);
b_m = bs(:, 2);


%% one-bit operation
for jj = 1:p
    G_j_op = diag(exp(-gammas(jj)^2/2 * (b_p - b_m).^2));
    XMPS = MPSOneSiteOp(XMPS, G_j_op, jj);
end


%% two-bit operation
Uswap = swap_gate(d);

fidSq = 1;

current_order = 1:p;
for round = 1:p
    
    start_ind = mod(round-1,2) + 1; % alternate 1, 2, 1, 2 between rounds

    leading_site_inds = start_ind:2:p-1;
    if start_ind == 2
        leading_site_inds = flip(leading_site_inds);
    end

    for jj = leading_site_inds

        my_j = current_order(jj);
        my_k = current_order(jj+1);
        G_jk_op = Uswap * G_jk(my_j, my_k, gammas, betas);

        [XMPS, err] = MPSTwoSiteOp(XMPS, G_jk_op, jj, false);

        fidSq  = fidSq*(1-err);

        % swap the qudits in the recorded ordering
        current_order(jj) = my_k;
        current_order(jj+1) = my_j;
    end
end


XMPS = MPSreverse(XMPS);


end

function op = G_jk(j, k, gammas, betas)
    if j > k
        l = j; j = k; k = l;
    end
    b_p = [1,1,-1,-1];
    b_m = [1,-1,1,-1];
    
    beta_factor = 1;
    for r = j:k-1
        beta_factor = beta_factor * cos(2*betas(r));
    end
    op = exp(-gammas(j) * gammas(k).*kron(b_p-b_m, b_p - b_m) * beta_factor);
    op = diag(op);
end