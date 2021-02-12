function MPS = MPScanonicalize(MPS)
%MPScanonicalize canonicalize the input MPS matrices
%
%   MPS = MPScanonicalize(MPS)

global d

N_sites = length(MPS.G);

for ind = 1:N_sites-1
    MPS = MPSTwoSiteOp(MPS, speye(d^2), ind, 0);
end

for ind = (N_sites-2):-1:1
    MPS = MPSTwoSiteOp(MPS, speye(d^2), ind, 0);
end


end

