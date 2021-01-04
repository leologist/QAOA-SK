function sigma = sigma_fun(a)
    global p
    sigma = 0;
    for r = 1:p
        sigma = sigma - a(:,r) + a(:, 2*p+1-r);
    end
    sigma = sigma/2;
end