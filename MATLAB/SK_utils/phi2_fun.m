function phi2 = phi2_fun(c, gammas)
%phi2_fun computes Phi_{c}^2 for configuration c
%
    global p
    phi2 = 0;
    for r = 1:p
        phi2 = phi2 + gammas(r)*(c(:,r) - c(:, 2*p+1-r));
    end
    phi2 = phi2.^2;
end