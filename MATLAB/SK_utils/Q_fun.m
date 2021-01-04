function Q = Q_fun(a, betas)
%qQ_fun computes the product of Q(a) for configuration a at given betas
%
%   input a can be L*2p matrix, wheree each row is a configuration
%
	global p
    Q = (1i).^sigma_fun(a);
    for r = 1:p
        powers = (a(:,r) + a(:,2*p+1-r))/2;
        Q = Q .* sin(betas(r)).^(1-powers) .* cos(betas(r)).^(1+powers);
    end
end