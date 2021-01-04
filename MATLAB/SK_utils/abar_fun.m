function abar = abar_fun(a)
%abar_fun bars the given bit-configuration(s)
%   abar = abar_fun(a)

    global abar_ops p
    partID = mod(config2hash(a), 2^p) + 1;
    abar = abar_ops(partID, :) .* a;
end