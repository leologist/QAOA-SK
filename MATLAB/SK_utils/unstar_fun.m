function aprime = unstar_fun(a)
    p = size(a,2)/2;
    aprime = a;
    for r = 1:p-1
        aprime(:, r) = a(:, r) .* a(:, r+1);
        aprime(:, 2*p+1-r) = a(:, 2*p+1-r) .* a(:, 2*p+1-r-1);
    end
end