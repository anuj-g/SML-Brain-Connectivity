function [MeanPCFNodeMales, errorMales, MeanPCFNodeFemales, errorFemales] = partcf_normalized_errorbars(normdir)
%Process Sex:0

%Clear variables
maledir = strcat(normdir, '/high/');
mfiles = dir(strcat(maledir, '*.mat'));

for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    CI = modularity_dir(M.fibergraph);
    %Calculate the Clustering Coeff of each node
    P = participation_coef(M.fibergraph, CI);
    %Put it as a row in a matrix
    OM(k,:) = P;
end

%Find the mean ccf of each node(mean over each column)
MeanPCFNodeMales = mean(OM);

%Std error = stddev/sqrt(Number of points)
% Calculate std dev of ccf of each node and divide it by sqrt of number of measurements 
errorMales = std(OM)/sqrt(size(OM,1));

%Plot the mean ccf of each node with error bars in blue
%errorbarMale = errorbar(1:prod(size(MeanPCFNodeMales)), MeanPCFNodeMales,errorMales,'.k', 'color', 'blue');

%Use the same figure for further plots
%hold on;

femaledir = strcat(normdir, '/low/');
femfiles = dir(strcat(femaledir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femaledir, femfiles(k).name));

    CI = modularity_dir(M.fibergraph);
    %Get the clustering coeffecient
    P = participation_coef(M.fibergraph, CI);
    %Append the obtained row to a file
    OF(k,:) = P;
end
%Calculate mean ccf of each node(mean of each column)
MeanPCFNodeFemales = mean(OF);
%Calculate standard error(for each node)
errorFemales = std(OF)/sqrt(size(OF,1));
%Plot this with error bars in red
%{
errorbarFemale = errorbar(1:prod(size(MeanPCFNodeFemales)), MeanPCFNodeFemales, errorFemales,'.k', 'color', 'red');
xlabel('Brain region',  'FontSize',14);
ylabel('Participation Coefficient',  'FontSize',14);
title('Mean of participation coefficient for each brain region',  'FontSize',16);
legend('Males', 'Females');

hold off;
clearvars;
%}