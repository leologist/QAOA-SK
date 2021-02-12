function MPS1 = getAllOnesMPS(N, d)
%getAllOnesMPS returns a simple MPS that evaluates to 1 for all indices
%
%   MPS1 = getAllOnesMPS(N, d)
%       N qudits

    G = cell(N,1);
    L = cell(N-1,1);
    
    for ind = 2:N-1
        G{ind} = ones(1,1,d);
    end
    G{1} = ones(1,d);
    G{N} = ones(1,d);
    
    for ind = 1:N-1
        L{ind} = 1;
    end
    MPS1.G = G;
    MPS1.L = L;
end