function a = hash2config(hash)
%hash2config Returns a list of bit-configurations given a list of hashes

    global p
    a = fliplr(1-2*de2bi(hash,2*p));
    a(:, p+1:end) =  a(:, p+1:end).*fliplr(a(:, 1:p));
end