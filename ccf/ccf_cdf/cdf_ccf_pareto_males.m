function  cdf_ccf_pareto(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/males/');
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

hold on;
pareto(MeanCCFNodeMales);

xlabel('Brain Region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Ranking of brain region according to mean clustering coefficient (for males)',  'FontSize',16);
hold off;

clearvars;