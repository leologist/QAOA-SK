function MPS = initMPSstate(N, d, k)
%initMPSstate(N, d, k)
%   [Gammas Lambdas] = initMPSstate(N, d, k)
%   creates tensor product of |kkkkkkk> in N qudits if k is a integer
%   creates tensor product of |k1k2k3...> if k is a 1xN vector of integers
%   creates tensor product of |v1>|v2>... if k is a dxN matrix of the form
%       k = [v1 v2 v3 ... vN], vi \in \mathbb{C}^d

if length(k) == 1
    k = k*ones(1, N);
    temp = zeros(d, 1);
    temp(k) = 1;
    kvecs = repmat(temp, 1, N);
elseif size(k) == [1, N]
    kvecs = [];
    for ind = 1:N
        temp = zeros(d,1);
        temp(k(ind)) = 1;
        kvecs = [kvecs, temp];
    end
elseif all(size(k) == [d, N])
    kvecs = k;
else
    error('input k not valid');
end

%%
Gammas = cell(N, 1);
Lambdas = cell(N-1, 1);

for ind = 1:N-1
    Lambdas{ind} = 1;
end

for ind = 2:N-1
    Gammas{ind} = zeros(1,1,d);
    Gammas{ind}(1,1,:) = kvecs(:, ind);
end

Gammas{1} = zeros(1, d);
Gammas{1}(1, :) = kvecs(:, 1);

Gammas{N} = zeros(1, d);
Gammas{N}(1, :) = kvecs(:, N);

MPS.G = Gammas;
MPS.L = Lambdas;

end