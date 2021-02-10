%% needs to have run initialize.m so that ../ and ../SK_utils are in MATLAB path
if ~exist('Delta_fun','file')
    cd ..
    initialize;
    cd experimental
end

global p aList

p = 4;


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
%
%       uniqueBarInds_sorted are indices of rows of aList correspondning
%       to configurations in D, where A - A_{p+1} = D \cup D^c

% short hand for generating vectors of a_{j} +/- a_{-j}

fp = @(ind) (aList(:, ind) + aList(:, 2*p+1-ind))/2;
fm = @(ind) (aList(:, ind) - aList(:, 2*p+1-ind))/2;

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


%% get the essential Low-Rank Basis vectors
[myU, myV, myS] = get_Delta_SVD(aList, gammas);

Ls = nan(4^p, 2*p);
for ind = 1:p
    Ls(:, ind) = fp(ind);
    Ls(:, p+ind) = fm(ind);
end

Ls = sparse(Ls);
myU = sparse(myU);
myV = sparse(myV);
Lbasis = [myU, myV, Ls];
dimL = size(Lbasis,2);

%% verify obj function value

temp = Ls.'*Ws;
obj = 2i*gammas*(temp(1:p).*temp(p+1:end));

fprintf('p=%d, obj = %0.6f, compared to standard %0.6f\n', p, real(obj), -SK_inf(p).fval)



%% FAST APPROXIMATE objective function evaluation in Low-Rank basis

Xtau = Lbasis.'*Xs;

tic;
J1 = nan(dimL, Drank);
for ind = 1:Drank
    J1(:, ind) = (myV(:,ind).*Lbasis)'*Xs;
end
fprintf('J1 done after %0.2f s\n', toc);

tic;
J2 = nan(dimL, Drank, Drank);
for ind = 1:Drank
    for ind2 = 1:Drank
        J2(:, ind, ind2) = (myV(:,ind).*myV(:,ind2).*Lbasis)'*Xs;
    end
end
J2 = reshape(J2, [dimL, Drank^2]);
J2 = sparse(J2);
fprintf('J2 done after %0.2f s\n', toc);

tic;
J3 = nan(dimL, Drank, Drank, Drank);
for ind = 1:Drank
    for ind2 = ind:Drank
        for ind3 = ind2:Drank
            temp = (myV(:,ind).*myV(:,ind2).*myV(:,ind3).*Lbasis)'*Xs;
            temp(abs(temp) < 1e-14) = 0;
            for v = perms([ind, ind2, ind3])'
                J3(:, v(1), v(2), v(3)) = temp;
            end
        end
    end
end
J3 = reshape(J3, [dimL, Drank^3]);
J3 = sparse(J3);
fprintf('J3 done after %0.2f s\n', toc);


% tic;
% J4 = nan(dimL, Drank, Drank, Drank, Drank);
% for ind = 1:Drank
%     for ind2 = ind:Drank
%         for ind3 = ind2:Drank
%             for ind4 = ind2:Drank
%                 temp = (myV(:,ind).*myV(:,ind2).*myV(:,ind3).*myV(:,ind4).*Lbasis)'*Xs;
%                 for v = perms([ind,ind2, ind3, ind4]).'
%                     J4(:, v(1), v(2), v(3), v(4)) = temp;
%                 end
%             end
%         end
%     end
% end
% J4 = reshape(J4, [dimL, Drank^4]);
% J4 = sparse(J4);
% fprintf('J4 done after %0.2f s\n', toc);
% 
% tic;
% J5 = nan(dimL, Drank, Drank, Drank, Drank, Drank);
% for ind = 1:Drank
%     for ind2 = ind:Drank
%         for ind3 = ind2:Drank
%             for ind4 = ind3:Drank
%                 for ind5 = ind4:Drank
%                     temp = (myV(:,ind).*myV(:,ind2).*myV(:,ind3).*myV(:,ind4).*myV(:,ind5).*Lbasis)'*Xs;
%                     for v = perms([ind,ind2, ind3, ind4, ind5]).'
%                         J5(:, v(1), v(2), v(3), v(4), v(5)) = temp;
%                     end
%                 end
%             end
%         end
%     end
% end
% J5 = reshape(J5, [dimL, Drank^5]);
% J5 = sparse(J5);
% fprintf('J5 done after %0.2f s\n', toc);



%%

tic;
Wtau = Xtau + normrnd(0,1,size(Xtau));
for ind = 1:p
    W2 = Wtau(1:Drank).*myS/2;
    Wtau = Xtau + J1*W2 + J2*kron(W2,W2)/2 + J3*krons(W2, W2, W2)/6;
    
    % (optional) higher order corrections:
    % Wtau = Wtau + J4*krons(W2, W2, W2, W2)/24  + J5*krons(W2, W2, W2, W2, W2)/120;
end

myW = Wtau;
temp = myW(2*Drank + (1:2*p)); % extract the Ws that correspond to Ls
approx_val = 2i*gammas*(temp(1:p).*temp(p+1:end));

fprintf('Low-Rank Basis evaluation time = %0.6f s\n', toc);
fprintf('  Result  = %0.6f %+0.6fj\n', real(approx_val), imag(approx_val));
fprintf(' True val = %0.6f\n',  obj);
fprintf('     Diff = %0.2e\n',   abs(approx_val - obj));


return

%% ============ EXTRAS BELOW ============

%% plotting W in Low-Rank basis

set(0,'DefaultAxesFontSize', 16);

W_Ls = nan(2*p, 1);
for ind = 1:p
    W_Ls(ind) = fp(ind)'*Ws;
    W_Ls(p+ind) = fm(ind)'*Ws;
end

W_Us = myU.'*Ws;
W_Vs = myV.'*Ws;

figure(4);
subplot(2,1,1);
plot(abs(Ws),'o')
grid on
set(gca, 'yscale','log','xlim',[0,4^p])
ylabel('$|W_a|$','Interpreter','latex');
xlabel('$a$ (Configuration Basis)','Interpreter','latex');

subplot(2,1,2);
plot(1:2*p, abs(W_Ls),'o')
hold on
plot(1:(p-1)^2, abs(W_Us), 'x')
% plot(1:(p-1)^2, abs(W_Vs), 's')
hold off, grid on
xlabel('$\lambda$ (Low-Rank Basis)','Interpreter','latex');
ylabel('$|\tilde{W}^\lambda|$','Interpreter','latex');


%% obtaining W as fixed point of evolution

Wt = normrnd(0, 1, 4^p, 1)/2^p;
% Wt = Xs;

%%% exact
myfun = @(w) exp(Dmat.'*w/2) .* Xs;

%%% approximate
% myfun = @(w) Xs + (Dmat.'*w/2) .* Xs ...
%     + (Dmat.'*w/2).^2/2 .* Xs ...
%     + (Dmat.'*w/2).^3 / 6 .* Xs + (Dmat.'*w/2).^4 / 24 .* Xs;

for ind = 1:p-1
    Wt = myfun(Wt);
end
fprintf('  ||W - W_t||  = %0.6f\n', norm(Ws - Wt));
fprintf('||W - fun(W)|| = %0.6f\n', norm(Ws - myfun(Ws)));


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

%% tensor factorization: Schmidt decomposition of X_b = X_{b1, b2, ...}

subs = cell(1,2*p);
for ind = 1:p
    subs{2*ind-1} = (1-aList(:,ind))/2+1;
    subs{2*ind} = (1-aList(:,2*p+1-ind))/2+1;
%     subs{2*ind-1} = (1-aPrime(:,ind))/2+1;
%     subs{2*ind} = (1-aPrime(:,2*p+1-ind))/2+1;
end
inds = sub2ind(ones(1,2*p)*2, subs{:});

%%% tensor index order:
% [a_1, a_{-1}, a_2, a_{-2}, a_3, a_{-3}, ..., a_p, a_{-p}]
Xtensor = nan(ones(1,2*p)*2);

Xtensor(inds) = Xs;
Xtensor(inds) = Qs; % for low rank example, set X = Q

Wtensor = Xtensor;
Wtensor(inds) = Ws;

%%% alternate reordering
% Wtensor = permute(Wtensor, [2*(1:p)-1, 2*(p:-1:1)]);
% --> [1,2, ..., p, -p, ..., -2, -1]

% Schmidt decomposition (SVD)

k = 2;
Xsvd = svd(reshape(Xtensor, [4^k, 4^(p-k)]));
Xsvd2 = svd(reshape(Wtensor, [4^k, 4^(p-k)]));

Xsvd(Xsvd < 1e-14) = 0; % get rid of numerical artifacts
Xsvd2(Xsvd2 < 1e-14) = 0;

figure(9);
semilogy(Xsvd,'o')
hold on
semilogy(Xsvd2,'x')
hold off, grid on
xlabel('index')
ylabel('Schmidt values')
legend('X tensor', 'W tensor')
title(sprintf('Cut [1...%d], [%d..%d] -- rank <= %d', k, k+1, p, min(4^k, 4^(p-k))));


%% finding new basis vectors

LLs = nan(4^p, p*(p-1)*(p-2)/6);
myind = 1;
for ind = 1:p-2
    for ind2 = ind+1:p-1
        for ind3 = ind2+1:p
            LLs(:, myind) = fp(ind).*fm(ind2).*fm(ind3);
            myind = myind + 1;
        end
    end
end
LLs'*[myU, myV, Ls]
