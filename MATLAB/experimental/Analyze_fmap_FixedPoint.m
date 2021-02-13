global p aList

p = 5;

load('../SK_opts.mat');
param0 = SK_inf(p).param;
gammas = param0(1:p);
betas = param0(p+1:end);

% gammas = linspace(0.22,0.64,p);
% betas = linspace(0.65,0.1,p);

Drank = (p-1)^2;

fprintf('**** p=%d ****\n', p);

%% setting up

SK_QAOA_p_objfun_helper;

% Note: aList is 4^p x 2p matrix, where each row is a configuration
%       organized from top to bottom as configurations from
%       A_{p+1}, A_{p}, ..., A_2, A_1 [see Section 6.1 of arXiv:1910.08187]

%%
aPrime = unstar_fun(aList);
Dmat = Delta_fun(aPrime, aPrime, gammas);


%% prepare Xs(aPrime)

Qs = Q_fun(aPrime, betas);
Fs = F_fun(aPrime, gammas, betas);
Xs = Qs.*Fs;


%% get Ws(aPrime) - for exact calculation of objective function
tic;

Ws = Xs;
current_ind = 2^p;
mysize = 2^p;

% iterate over A_{p}, ... A_2, A_1
for ind = 1:p 
    new_inds = current_ind + (1:mysize);
    Ws(new_inds) = exp(Dmat(1:current_ind, new_inds).' * Ws(1:current_ind) / 2) .* Xs(new_inds);
    current_ind = current_ind + mysize;
    mysize = 2*mysize;
end

fprintf('Ws are made after %0.6f s\n', toc);


%% construct L_\pm basis vectors

tp = @(ind) (aList(:, ind) + aList(:, 2*p+1-ind))/2;
tm = @(ind) (aList(:, ind) - aList(:, 2*p+1-ind))/2;

Ls = nan(4^p, 2*p);
for ind = 1:p
    Ls(:, ind) = tp(ind);
    Ls(:, p+ind) = tm(ind);
end

%% obtaining W as fixed point of evolution

% Wt = normrnd(0, 1, 4^p, 1)/2^p; % random roughly unit vector in 4^p dim
% Wt = Qs;
Wt = Xs;

%%% exact
myfun = @(w) exp(Dmat.'*w/2) .* Xs;

%%% approximate
% myfun = @(w) Xs + (Dmat.'*w/2) .* Xs ...
%     + (Dmat.'*w/2).^2/2 .* Xs ...
%     + (Dmat.'*w/2).^3 / 6 .* Xs + (Dmat.'*w/2).^4 / 24 .* Xs;

W_L = Ls.'*Wt;
myobj = 2i*gammas*(W_L(1:p).*W_L(p+1:end));
fprintf('====== fmap iter begins ======\n')
fprintf('-- init obj_guess = %0.6f\n', myobj);
    
for ind = 1:p
    temp = myfun(Wt);
    
    W_L = Ls.'*temp;
    temp2 = 2i*gammas*(W_L(1:p).*W_L(p+1:end));
    
    fprintf('-- after iter %d, |?W| = %0.3e, ?V= %0.3e\n----------- obj_guess = %0.6f\n', ...
        ind, norm(temp-Wt), temp2-myobj, temp2); 
    Wt = temp;
    myobj = temp2;
end
fprintf('  ||W - W_t||  = %0.6e\n', norm(Ws - Wt));


%% closer look at fixed point evolution - error various order
temp = Dmat.'*Ws/2;

figure(5);

maxorder = 3;
errs = nan(maxorder+1,1);
yy = - Ws;
for ind = 0:maxorder
    yy = yy + temp.^ind / factorial(ind) .*Xs;
    plot(sort(abs(yy)),'.','markersize',10)
    hold on
    errs(ind+1) = norm(yy);
end
hold off, grid on
legend(arrayfun(@(n, err) sprintf('%d-order err = %0.2e', n, err), (0:maxorder)', errs, 'UniformOutput',0))
legend('location', 'northwest')
set(gca,'yscale','log')
xlabel('index b')
ylabel('$W_b - [\hat{f}^{(k)}(W)]_b X_b$', 'Interpreter','latex','fontsize',24)
