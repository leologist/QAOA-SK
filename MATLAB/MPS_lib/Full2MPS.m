function MPS = Full2MPS(psi, din)
%% Need to test more - doesn't quite work

global d
if nargin >= 2
    d = din;
end
svdecon = @(X) svd(X, 'econ');

N = reallog(numel(psi))/reallog(d);
if N ~= round(N)
    error('Incompatible dimension: dim(psi) = %d, d = %d', numel(psi), d);
end

G = cell(N,1);
L = cell(N-1, 1);

psi = full(psi);

temp = reshape(psi, [], d);


[U, S, V] = svdecon(temp);
Ikeep = abs(diag(S)) > 1e-10;
U = U(:, Ikeep); S = S(Ikeep, Ikeep); V = V(:, Ikeep);

L{1} = S;
DR = length(S);
G{1} = reshape(V', [DR, d]);


Phis = U;
DL = DR;

for site_ind = 2:N-1
    Phis = reshape(Phis, [], d* DL);
    if nnz(isnan(Phis) | isinf(Phis))
        pause
    end
    [U, S, V] = svdecon(Phis);
    
    Ikeep = abs(diag(S)) > 1e-10;
    U = U(:, Ikeep); S = S(Ikeep, Ikeep); V = V(:, Ikeep);
%     [size(U), size(S), size(V)]

    DR = length(S);
    L{site_ind} = S;
    G{site_ind} = permute(reshape(V, [d, DL, DR]), [2,3,1]);
    Phis = U;
    DL = DR;
end
G{N} = U';


MPS.G = G;
MPS.L = L;




end