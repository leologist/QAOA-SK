if ~exist('Delta_fun', 'file')
    addpath([pwd, filesep, 'SK_utils'])
end

global p

p = 3;

gammas = rand(1,p);
betas = rand(1,p);

%% setting up

SK_QAOA_p_objfun_helper;

global aList bigset uniqueBarInds_sorted

% Note: aList is 4^p x 2p matrix, where each row is a configuration
%       organized from top to bottom as configurations from
%       A_{p+1}, A_{p}, ..., A_2, A_1 [see Section 6.1 of arXiv:1910.08187]
%
%       uniqueBarInds_sorted are indices of rows of aList correspondning
%       to configurations in D, where A - A_{p+1} = D \cup D^c

%% Efficiently construct the X_u

tic;

Qs_D = Q_fun(aList(uniqueBarInds_sorted,:), betas);
Fs_D = F_fun(aList(uniqueBarInds_sorted,:), gammas, betas);

Xs_D = Qs_D.*Fs_D; % X(u) = q*Q(u)*F(u)

Xs = nan(4^p,1);
Xs(1:2^p) = Q_fun(bigset, betas);
Xs(uniqueBarInds_sorted) = Xs_D;
Xs(uniqueBarBarInds_sorted) = -Xs_D;
toc;

fprintf('Done making Xs after %0.6f s\n', toc);

%% Construct matrix of Delta_{a,b} and K_{a,b} = X_a Delta_{a,b}
%   Note Delta_{a,b} = (Phi_{(a bar)* b}^2 - Phi_{a*b}^2)/2

tic;
Dmat = Delta_fun(aList, aList, gammas);

Kmat = Xs.*Dmat;
toc;

% reorder the matrices so the increasing indices goes through A_1, A_2 ... A_{p+1}
Dmat = rot90(Dmat, 2);
Kmat = rot90(Kmat, 2);

Drank = rank(Dmat);

fprintf('p=%d --- rank(D) = rank(K) = %d\n', p, Drank);

%% singular decompositions
[Ud, Sd, Vd] = svd(Dmat);
Sd = diag(Sd);

[Uk, Sk, Vk] = svd(Kmat);
Sk = diag(Sk);

%% plot singular values
set(0,'DefaultAxesFontSize', 14);

figure(11);
subplot(2,1,1);
plot(Sd,'ob')
grid on
ylabel('singular values of \Delta')

subplot(2,1,2);
plot(Sk,'xr');
grid on
ylabel('singular values of K')

xlabel('1 \leq index \leq 4^p')

%% plot singluar vectors of Delta matrix

figure(12);
x = (1:4^p)/(2^p);
for ind = 1:Drank
    vec = Ud(:,ind);
    vec = vec/max(abs(vec));
    
%     vec(abs(vec)<1e-10) = nan; % uncomment if plot non-zero entries only

    subplot(p-1, p-1, ind)
    plot(x, vec, '.','markersize',10);
    hold on
    plot([1;1]* (1:2^p), [-1,1]*2,'--r','LineWidth',1);
    hold off, grid on
    xlabel('index / 2^p ')
    title(sprintf('singular vector # %d', ind));
    set(gca,'xlim',[0,2^p],'ylim',[-1,1]*1.4);
end

