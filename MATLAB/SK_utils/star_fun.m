function astar = star_fun(a)
    p = size(a,2)/2;
    astar = a;
    for r = p-1:-1:1
        astar(:, r) = astar(:,r) .* astar(:,r+1);
        astar(:, 2*p+1-r) = astar(:, 2*p+1-r) .* astar(:, 2*p+1-r-1);
    end
end