function [MeanEBCMales, errorEBCMales, MeanEBCFemales, errorEBCFemales] = edge_betweenness_normalized_errorbars_final(normdir, edgelow, edgehigh)
%Process Sex:0

%Clear variables
maledir = strcat(normdir, '/high/');
mfiles = dir(strcat(maledir, '*.mat'));

for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    D = M.fibergraph;
    E= find(D);
    D(E) = 1./D(E);
    EBC = edge_betweenness_wei(D);
    RS = reshape(EBC, prod(size(EBC)),1);
    OM(k,:) = RS;
end
MeanEBCMales = mean(OM);
errorEBCMales = std(OM)/sqrt(size(OM,1));

femaledir = strcat(normdir, '/low/');
femfiles = dir(strcat(femaledir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femaledir, femfiles(k).name));
    D = M.fibergraph;
    E = find(D);
    D(E) = 1./D(E);
    EBC = edge_betweenness_wei(D);
    RS = reshape(EBC, prod(size(EBC)),1);
    OF(k,:) = RS;
end
MeanEBCFemales = mean(OF);
errorEBCFemales = std(OF)/ sqrt(size(OF,1));
%{
figure;
hold on;
errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red')
xlabel('Inter-region Edge',  'FontSize',14);
ylabel('Edge Betweenness',  'FontSize',14);
title('Mean of edge betweenness for each inter-region edge',  'FontSize',14);
legend('Males', 'Females');

hold off;
clearvars;
%}