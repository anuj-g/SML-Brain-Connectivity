function  ccf_normalized_median(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/males/');
mfiles = dir(strcat(maledir, '*.mat'));

for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    %Calculate clustering coeffecient
    ccf = clustering_coef_wd(M.fibergraph);
    %Put the row(clustering coff of each node) into a row of a matrix 
    OM(k,:) = ccf;
end

%Calaculate the Median CCF of each node(Mean of each column)
MedianCCFNodeMales = median(OM);

%Plot this using a Blue color
plot(MedianCCFNodeMales,'color', 'blue');

%Use the same figure for further plots
hold on;

% Do the exact same thing for other sex as well
femdir = strcat(normalizeddir, '/females/');
femfiles = dir(strcat(femdir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femdir, femfiles(k).name));
    %Calculate clustering coefficient
    ccf = clustering_coef_wd(M.fibergraph);
    %Put the returned array of CCF for each node in a matrix
    OF(k,:) = ccf;
end

%Find the median of each node
MedianCCFNodeFemales = median(OF);
%Plot it using red
plot(MedianCCFNodeFemales,'color', 'red');

xlabel('Brain region', 'FontSize',14);
ylabel('Clustering Coefficient', 'FontSize',14);
title('Median of clustering coefficient for each brain region', 'FontSize',16);
legend('Males', 'Females');

hold off;
clearvars;