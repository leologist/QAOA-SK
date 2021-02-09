if ~exist('Delta_fun', 'file')
    addpath([pwd, filesep, 'SK_utils'])
    addpath(pwd);
end

global p abar_ops aList aHashes partitions bigset
global uniqueBarInds_sorted this_aHashes this_aBarHashes
global quickPrimer quickBarrer MrPlus MrMinus PlaysWellTruthTable

tic;

halfList = 1-2*de2bi((0:2^p-1));
halfList = fliplr(halfList);

partitions = cell(2^p, 1);
abar_ops = ones(2^p, 2*p); 
for ind = 1:2^p
    partitions{ind} = [halfList, fliplr(halfList).*halfList(ind, :)];
    temp = find(halfList(ind, :) == -1, 1, 'first');
    abar_ops(ind, p+1-temp) =  -1;
    abar_ops(ind, p+temp) =  -1;
end

bigset = partitions{1};

%% quick helper functions

% list of configurations - 2^(2p) x 2p matrix
aList = cat(1, partitions{:});

% hashes of the configurations, which are unique integer representation of configurations
aHashes = config2hash(aList);

% create quick functions for manipulating configurations
aBarHashes = config2hash(abar_fun(aList));
aStarHashes = config2hash(star_fun(aList));
aPrimeHashes = config2hash(unstar_fun(aList));

quickPrime_helper = sortrows([aHashes, aPrimeHashes]);
quickPrime_helper = quickPrime_helper(:,2);
quickPrimer = @(hash) quickPrime_helper(hash + 1);

quickBar_helper = sortrows([aHashes, aBarHashes]);
quickBar_helper = quickBar_helper(:,2);
quickBarrer = @(hash) quickBar_helper(hash + 1);

aHashToInd = sortrows([aHashes, (1:2^(2*p))']);
quickHashToInder = @(hash) aHashToInd(hash+1, 2);

% bigsetHash = config2hash(bigset);
bigsetHash = (0:2^p-1)*2^p;

uniqueBarHashes = setdiff(unique(min(aHashes, aBarHashes)), bigsetHash);
uniqueBarInds = quickHashToInder(uniqueBarHashes);
[uniqueBarInds_sorted, IsortUnique] = sort(uniqueBarInds);
this_aHashes = aHashes(uniqueBarInds_sorted);
this_aBarHashes = quickBarrer(this_aHashes);
uniqueBarBarInds_sorted = quickHashToInder(this_aBarHashes);

fprintf('List of configuration, hashes, and helper functions made after %0.6f s\n', toc);

%% Constructing sets

tic;
MrPlus = cell(p,1);
MrMinus = cell(p, 1);
for r = 1:p
    MrPlus{r} = config2hash(aList(aList(:, r) == aList(:, 2*p+1-r), :));
    MrMinus{r} = config2hash(aList(aList(:, r) == -aList(:, 2*p+1-r), :));
end
fprintf('The sets Mr + and - constructed after %0.6f s\n', toc);


%% Get truth table of who plays well with whom

tic;
% get representative config from each partition
As = ones(2^p, 2*p);
for ind = 1:2^p
    As(ind, :) = partitions{ind}(1,:);
end

PlaysWellTruthTable = Delta_fun(As, As, rand(1,p)) ~= 0;
fprintf('Plays Well Truth Table made after %0.6f s\n', toc);

%%
fprintf('*** SK QAOA precomputing done for p=%d\n', p);