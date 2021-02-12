function MPS = MPSOneSiteOp(MPS, Op, k)
%MPSOneSiteOp apply a single site operator to MPS
%
%   MPS = MPSOneSiteOp(MPS, Op, k)

global d

G_in = MPS.G;
N_sites = length(G_in);

G_out = G_in;

if k == 1 || k == N_sites
    G_out{k} = G_in{k}* (Op.');
else
    DL = size(G_in{k},1);
    DR = size(G_in{k},2);
    G_out{k} = reshape(reshape(G_in{k},[DL*DR, d])* (Op.'), [DL, DR, d]);
end

MPS.G = G_out;

end