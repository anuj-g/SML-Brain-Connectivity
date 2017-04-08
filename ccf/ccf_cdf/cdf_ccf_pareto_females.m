function  cdf_ccf_pareto_females(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/females/');
femmfiles = dir(strcat(maledir, '*.mat'));
for k = 1:numel(femmfiles)
    M = load(strcat(maledir, femmfiles(k).name));
    %Calculate the Clustering Coeff of each node
    ccf = clustering_coef_wd(M.fibergraph);
    %Put it as a row in a matrix
    OM(k,:) = ccf;
end

%Find the mean ccf of each node(mean over each column)
MeanCCFNodeFemales = mean(OM);

hold on;
pareto(MeanCCFNodeFemales);

xlabel('Brain Region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Ranking of brain region according to mean clustering coefficient (for females)',  'FontSize',16);
hold off;

clearvars;