function [MeanCCFNodeMales, errorMales, MeanCCFNodeFemales, errorFemales] = ccf_normalized_errorbars(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/high/');
mfiles = dir(strcat(maledir, '*.mat'));
for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    %Calculate the Clustering Coeff of each node
    ccf = clustering_coef_wd(M.fibergraph);
    %Put it as a row in a matrix
    OM(k,:) = ccf;
end

%Find the mean ccf of each node(mean over each column)
MeanCCFNodeMales = mean(OM);

%Std error = stddev/sqrt(Number of points)
% Calculate std dev of ccf of each node and divide it by sqrt of number of measurements 
errorMales = std(OM)/sqrt(size(OM,1));

%Plot the mean ccf of each node with error bars in blue
%errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');

%Use the same figure for further plots
%hold on;

% Do the exact same thing for other sex as well
femdir = strcat(normalizeddir, '/low/');
femfiles = dir(strcat(femdir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femdir, femfiles(k).name));
    %Get the clustering coeffecient
    ccf = clustering_coef_wd(M.fibergraph);
    %Append the obtained row to a file
    OF(k,:) = ccf;
end
%Calculate mean ccf of each node(mean of each column)
MeanCCFNodeFemales = mean(OF);
%Calculate standard error(for each node)
errorFemales = std(OF)/sqrt(size(OF,1));
%Plot this with error bars in red
%{
errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales, errorFemales,'.k', 'color', 'red');

xlabel('Brain region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Mean of clustering coefficient for each brain region',  'FontSize',16);
legend('high', 'low');

hold off;
clearvars;
%}