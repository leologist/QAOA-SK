function MPS = MPSreverse(MPS)
%MPSreverse reverses the site ordering in the MPS
%
%   MPS = MPSreverse(MPS)

    MPS.L = flip(MPS.L);
    Gout = flip(MPS.G);
    for ind = 2:length(Gout)-1
        Gout{ind} = permute(Gout{ind}, [2,1,3]);
    end
    MPS.G = Gout;
end