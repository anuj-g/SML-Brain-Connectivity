function  pvalue = bootstrap_partcf(normalizeddir, nodeno, threshold)
%normalizeddir: Directory which contains normalized data
%nodeno: Node no whose participation coefficient needs to be analyzed if
%it is statistically significant
%thresh-hold: The observed difference.

% Goto the location where both the classes have been merged

% Read all the files into an array
malemats = dir(strcat(normalizeddir, '/males/*.mat'));
mergeddir = strcat(normalizeddir, '/merged/');
files = dir(strcat(mergeddir, '*.mat')); 

for i=1:numel(files)
    fullfiles{i} = strcat(mergeddir, files(i).name);
    fullfilename = char(fullfiles{i});
    M  = load(fullfilename);
    CI = modularity_dir(M.fibergraph);
    P  = participation_coef(M.fibergraph, CI);
    PBOOT(i) = P(nodeno);
end

assert(numel(PBOOT) == numel(files))
% Bootstrap for 3000 iteration. Each iteration will generate a sample of size 114 (total number of subjects) through repeated sampling
BT = bootstrp(3000, @bootstrap_brain_func, PBOOT, numel(malemats), 'Option', statset('UseParallel', 'always'));
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
