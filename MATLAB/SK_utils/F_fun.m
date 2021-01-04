function F = F_fun(a, gammas, betas)

    global bigset
    Q = Q_fun(bigset, betas);
    
    F = 0;
    for ind = 1:size(bigset,1)
        F = F - phi2_fun(star_fun(a.*bigset(ind,:)), gammas)*Q(ind);
    end
    F = exp(F/2);

end