function EV = MPSoverlap(MPS1, MPS2)
%MPSoverlap calculates overlap between two MPS
%   EV = MPSoverlap(MPS1, MPS2)
%   return <MPS2 | MPS1>

global d

Gammas1 = MPS1.G;
Lambdas1 = MPS1.L;

Gammas2 = MPS2.G;
Lambdas2 = MPS2.L;

N_sites = length(Gammas1);
if ~(length(Lambdas1) == N_sites-1 && length(Gammas2) == N_sites && length(Lambdas2) == N_sites-1)
    error('Invalid input - mismatched # sites');
end



%% contracting tensors from left to right

% A = Dt1 x Db1 matrix
A = Lambdas1{1}*Gammas1{1}*Gammas2{1}'*Lambdas2{1};

for site_ind = 2:N_sites-1
    B = 0; %zeros(Ds1(site_ind), Ds2(site_ind));
    for phys_ind = 1:d % seems faster than doing a bunch of reshape/permute
        B = B + Gammas1{site_ind}(:,:,phys_ind).' * A * conj(Gammas2{site_ind}(:,:,phys_ind));
    end
    A = Lambdas1{site_ind}*B*Lambdas2{site_ind};
end

% C = Dtn x Dbn matrx
C = Gammas1{N_sites}*Gammas2{N_sites}';

% EV = full((sum(sum(A.*C)));
EV = full(A(:).'*C(:));



end
