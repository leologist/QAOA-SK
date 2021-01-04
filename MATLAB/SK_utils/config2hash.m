function hash = config2hash(a)
%config2hash Returns a list of integers (hashes) that uniquely identify a
%    bit-configuration

    global p
    a(:, p+1:end) = a(:, p+1:end).*fliplr(a(:, 1:p));
    hasher = fliplr(2.^(0:(2*p)-1));
    hash = (1-a)/2 * hasher';
end