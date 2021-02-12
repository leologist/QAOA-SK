function psi = MPS2Full(MPS)
%MPStoFull Convert MPS representation to full state vector (d^N-dimension)
%
%   psi = MPStoFull(MPS)
%
%   Input:
%     MPS = struct containing MPS in canonical form
%           required to have two fields:
%             MPS.G = a cell array containing N Gamma_{ab}^s tensors
%             MPS.L = a cell array containing N-1 diagonal matrices of
%                     Schmidt values

Gammas = MPS.G;
Lambdas = MPS.L;

N_sites = length(Gammas);
d = size(Gammas{1}, 2);

psi = zeros(d^N_sites, 1);

for ind = 1:d^N_sites
    physical_inds = wrev(de2bi(ind-1, N_sites, d));
    physical_inds = physical_inds + 1; %convert 0-based index to 1-based index
    
    coeff = Gammas{1}(:, physical_inds(1)).'*Lambdas{1}; %1xD
    for site_ind = 2:N_sites-1
        coeff  = coeff*Gammas{site_ind}(:,:,physical_inds(site_ind)) * Lambdas{site_ind};
    end
    coeff = coeff * Gammas{N_sites}(:, physical_inds(N_sites));
    if length(coeff) ~= 1
        error('What happened?');
    end
    if abs(coeff) > eps
        psi(ind) = coeff;
    end
end

end

