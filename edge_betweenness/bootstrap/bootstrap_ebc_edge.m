function  pvalue = bootstrap_ebc_edge(normalizeddir, edgeno, threshold)

% Goto the location where both the classes have been merged

% Read all the files into an array
malemats = dir(strcat(normalizeddir, '/males/*.mat'));

mergeddir = strcat(normalizeddir, '/merged/');
files = dir(strcat(mergeddir, '*.mat')); 

for i=1:numel(files)
    fullfiles{i} = strcat(mergeddir, files(i).name);
    fullfilename = char(fullfiles{i});
    M = load(fullfilename);
    D = M.fibergraph;
    E = find(D);
    D(E) = 1./D(E);
    EBC = edge_betweenness_wei(D);
    EBCRS = reshape(EBC, prod(size(EBC)), 1);
    EBCVEC(i) = EBCRS(edgeno);
end

assert(numel(EBCVEC) == numel(files))
% Bootstrap for 3000 iteration. Each iteration will generate a sample of size 114 (total number of subjects) through repeated sampling
BT = bootstrp(3000, @bootstrap_brain_func, EBCVEC, numel(malemats), 'Option', statset('UseParallel', 'always'));
% Plot the histogram of the test statistic
hold on;

hist(BT);
% Find out number of observations more extreme then the alternative hypothesis
C = (abs(BT) > threshold);
num = prod(size(find(C)));
num
% Calculate its probablity. This is the p-value
pvalue = num/3000;

xlabel('Test statistic', 'FontSize',14);
ylabel('Frequency', 'FontSize',14);
title('Test statistic distribution', 'FontSize',16);

hold off;
