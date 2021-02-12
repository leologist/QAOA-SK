function [gammas, betas] = get_ParamGuess(p, SK_inf)
%
p_in = 8;

interp_param = @(param, p_in, p_out) [interp1(linspace(0,1,p_in), param(1:p_in), linspace(0,1,p_out)), ...
     interp1(linspace(0,1,p_in), param(p_in+1:end), linspace(0,1,p_out))];


if nargin < 2
    load('SK_opts.mat', 'SK_inf')
end
param0 = SK_inf(p_in).param;

for p_out = p_in+1:p
    param0 = interp_param(param0, p_out-1, p_out);
end



gammas = param0(1:p);
betas = param0(p+1:end);



gammas0 = gammas;

load('param_guess_from_powfit.mat', 'fitres_bet', 'fitres_bet_f', 'fitres_b1')
qcutoff_g = 6;
qcutoff_b = 5;
for q = 1:p
    if q <= qcutoff_g
        gammas(q) = gamma_fun(q, p);
    else
%         gammas(q) = file.param(p_in-(p-q)); % retro
        gammas(q) = gammas(q) - (gammas0(qcutoff_g)-gammas(qcutoff_g))/(p-qcutoff_g)*(p-q); % lingap
    end
    
    if q <= qcutoff_b
        betas(p+1-q) = fitres_bet{q}(p+1-q);
    elseif p+1-q <= 6
        betas(p+1-q) = fitres_bet_f{p+1-q}(p);
    end
end
betas(1) = fitres_b1(p);
end