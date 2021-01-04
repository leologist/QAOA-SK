function aprime = unstar_fun(a)
%unstar_fun is the inverse of star_fun
%       in an earlier version `unstar' operation = `prime' operation
%
    p = size(a,2)/2;
    aprime = a;
    for r = 1:p-1
        aprime(:, r) = a(:, r) .* a(:, r+1);
        aprime(:, 2*p+1-r) = a(:, 2*p+1-r) .* a(:, 2*p+1-r-1);
    end
end