clear all;
clc;
datasetdir='C:\Users\minug\Desktop\SML-Project\brainnetworks\smallgraphs\';
maledatasetDir='C:\Users\minug\Desktop\SML-Project\brainnetworks\males\';
femaledatasetDir='C:\Users\minug\Desktop\SML-Project\brainnetworks\females\';
path=strcat(datasetdir,'*.mat');
malePath=strcat(maledatasetDir,'*.mat');
femalePath=strcat(femaledatasetDir,'*.mat');
outdir = strcat(datasetdir, '/normalized');
mkdir(outdir);
numberOfFiles = dir(path);
numofMaleFiles=dir(malePath);
numofFemaleiles=dir(femalePath);
ZM = [];
ZF=[];
for k = 1:numel(numofMaleFiles)
    file=strcat(strcat(maledatasetDir, numofMaleFiles(k).name));
    M = load(file);
    F = full(M.fibergraph);  
    F=F+F';
    normalizedRow = diag(1./sum(F,2))*F;
    maxWeight = max(normalizedRow(:));    
    MalenormalizedMatrix = normalizedRow./maxWeight;    
    ccfMales(k,:) = clustering_coef_wd(MalenormalizedMatrix);    
    %mean edge connectivity
    ZM = cat(3, ZM, MalenormalizedMatrix);    
    %participation coeff   
    %finds non zero values
    E= find(MalenormalizedMatrix);
    MalenormalizedMatrix(E) = 1./MalenormalizedMatrix(E);
    EBC = edge_betweenness_wei(MalenormalizedMatrix);
    RS = reshape(EBC, prod(size(EBC)),1);
    OM(k,:) = RS;     
    
end
MeanZM = mean(ZM, 3);
MeanCCFNodeMales = mean(ccfMales);
errorMales = std(ccfMales)/sqrt(size(ccfMales,1));
errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');
hold on;
MeanEBCMales = mean(OM);
errorEBCMales = std(OM)/sqrt(size(OM,1));

for l = 1:numel(numofFemaleiles)
    file=strcat(strcat(femaledatasetDir, numofFemaleiles(l).name));
    M1 = load(file);
    F1 = full(M1.fibergraph);  
    F1=F1+F1';
    normalizedRow = diag(1./sum(F1,2))*F1;
    maxWeight = max(normalizedRow(:));    
    FemalenormalizedMatrix = normalizedRow./maxWeight;
    ccfFemale(l,:) = clustering_coef_wd(FemalenormalizedMatrix);      
    ZF = cat(3, ZF, FemalenormalizedMatrix);
    E= find(FemalenormalizedMatrix);
    FemalenormalizedMatrix(E) = 1./FemalenormalizedMatrix(E);
    EBC = edge_betweenness_wei(FemalenormalizedMatrix);
    RS = reshape(EBC, prod(size(EBC)),1);
    OF(l,:) = RS;       
end
MeanZF = mean(ZF, 3);
MeanCCFNodeFemales = mean(ccfFemale);
%ccf bootstrap computation
MeanDiff=MeanCCFNodeFemales-MeanCCFNodeMales;
[g,t]=sort(MeanDiff,'descend');
for i=1:70    
    pval(i)=bootstrap_ccf_node(t(i),g(i));
end
consider=t(pval<0.05);
figure(1);
errorFemales = std(ccfFemale)/sqrt(size(ccfFemale,1));
errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales,errorFemales,'.k', 'color', 'red');
xlabel('Brain region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Mean of clustering coefficient for each brain region',  'FontSize',16);
legend('Males', 'Females');
hold off;
MeanEBCFemales = mean(OF);
errorEBCFemales = std(OF)/ sqrt(size(OF,1));
edgelow=1;
edgehigh=4900;
figure(2);
writetoPAJ(MeanZM,'MeanSex0',1);
function get_ccf_nodes()
    
end

writetoPAJ(MeanZF,'MeanSex1',1);
heatmap((MeanZF - MeanZM), 1:70, 1:70);
hold on;
xlabel('Brain region', 'FontSize',14);
ylabel('Brain region', 'FontSize',14);
title('Differences in mean connectivity for each edge between the connectomes', 'FontSize',16);
colorbar;
%colormap bone;
hold off;
figure(3);
hold on;
errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red')
xlabel('Inter-region Edge',  'FontSize',14);
ylabel('Edge Betweenness',  'FontSize',14);
title('Mean of edge betweenness for each inter-region edge',  'FontSize',14);
legend('Males', 'Females');
hold off;
