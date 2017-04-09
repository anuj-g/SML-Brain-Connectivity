function  pvalue = bootstrap_ccf_node(normalizeddir, nodeno, threshold)

% Goto the location where both the classes have been merged

% Read all the files into an array
malemats = dir(strcat(normalizeddir, '/high/*.mat'));

mergeddir = strcat(normalizeddir, '/merged/');
files = dir(strcat(mergeddir, '*.mat')); 

%fullfiles = cellstr();

for i=1:numel(files)
    fullfiles{i} = strcat(mergeddir, files(i).name);
    fullfilename = char(fullfiles{i});
    M = load(fullfilename);
    ccfvec = clustering_coef_wd(M.fibergraph);
    NodeCCF(i) = ccfvec(nodeno);
end

assert(numel(NodeCCF) == numel(files))
% Bootstrap for 3000 iteration. Each iteration will generate a sample of size 114 (total number of subjects) through repeated sampling
BT = bootstrp(3000, @bootstrap_brain_func, NodeCCF, numel(malemats), 'Option', statset('UseParallel', 'always'));
% Plot the histogram of the test statistic
hold on;

hist(BT);
% Find out number of observations more extreme then the alternative hypothesis
C = (abs(BT) > threshold);
num = prod(size(find(C)));

% Calculate its probablity. This is the p-value
pvalue = num/3000;
xlabel('Test statistic', 'FontSize',14);
ylabel('Frequency', 'FontSize',14);
title('Test statistic distribution', 'FontSize',16);

hold off;
