clear all;
clc;
tic;
maledatasetDir='C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\males\';
femaledatasetDir='C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\females\';
malePath=strcat(maledatasetDir,'*.mat');
femalePath=strcat(femaledatasetDir,'*.mat');
numofMaleFiles=dir(malePath);
numofFemaleiles=dir(femalePath);
ZM = [];
ZF=[];
for k = 1:numel(numofMaleFiles)
    file=strcat(strcat(maledatasetDir, numofMaleFiles(k).name));
    M = load(file);
    F = full(M.fibergraph);
    MalenormalizedMatrix = F;
    %clustering coeff for male calculated
    ccfMales(k,:) = clustering_coef_wd(MalenormalizedMatrix);
    
    %mean edge connectivity for male calculated
    ZM = cat(3, ZM, MalenormalizedMatrix);
    
    %participation coeff    
    CI = modularity_dir(MalenormalizedMatrix);
    %Calculate the Clustering Coeff of each node
    P = participation_coef(MalenormalizedMatrix, CI);
    %Put it as a row in a matrix
    OM_PPF(k,:) = P;
    
    %Edgebetweenness connectivity
    E= find(MalenormalizedMatrix);
    MalenormalizedMatrix(E) = 1./MalenormalizedMatrix(E);    
    EBC = edge_betweenness_wei(MalenormalizedMatrix);
    RS = reshape(EBC, prod(size(EBC)),1);
    OM(k,:) = RS;
    
end
MeanZM = mean(ZM, 3);
MeanCCFNodeMales = mean(ccfMales);
MeanEBCMales = mean(OM);
errorEBCMales = std(OM)/sqrt(size(OM,1));

%PPF
MeanPPFNodeMales = mean(OM_PPF);
errorPPFMales = std(OM_PPF)/sqrt(size(OM_PPF,1));
%PPF over

for l = 1:numel(numofFemaleiles)
    file=strcat(strcat(femaledatasetDir, numofFemaleiles(l).name));
    M1 = load(file);
    F1 = full(M1.fibergraph);
    FemalenormalizedMatrix = F1;
    ccfFemale(l,:) = clustering_coef_wd(FemalenormalizedMatrix);
    ZF = cat(3, ZF, FemalenormalizedMatrix);
    E= find(FemalenormalizedMatrix);
    FemalenormalizedMatrix(E) = 1./FemalenormalizedMatrix(E);   
    EBC = edge_betweenness_wei(FemalenormalizedMatrix);
    RS = reshape(EBC, prod(size(EBC)),1);
    OF(l,:) = RS;
    
    %PPF
     FemaleCI = modularity_dir(FemalenormalizedMatrix);
    %Get the clustering coeffecient
    P = participation_coef(FemalenormalizedMatrix, FemaleCI);
    %Append the obtained row to a file
    OF_PPF(l,:) = P;
   
end
MeanZF = mean(ZF, 3);
MeanCCFNodeFemales = mean(ccfFemale);
%ccf bootstrap computation
MeanCCFDiff=MeanCCFNodeFemales-MeanCCFNodeMales;
[g,t]=sort(MeanCCFDiff,'descend');

for i=1:70
    pval(i)=bootstrap_ccf_node('C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\', t(i),g(i));
end
consider=t(pval<0.05);

errorMales = std(ccfMales)/sqrt(size(ccfMales,1));
errorBarCCFMale = errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');
maxCCFMale = MeanCCFNodeMales + errorBarCCFMale.LData;
minCCFMale = MeanCCFNodeMales - errorBarCCFMale.UData;

errorFemales = std(ccfFemale)/sqrt(size(ccfFemale,1));
errorBarCCFFemale = errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales,errorFemales,'.k', 'color', 'red');
maxCCFFemale = MeanCCFNodeFemales + errorBarCCFFemale.LData;
minCCFFemale = MeanCCFNodeFemales - errorBarCCFFemale.UData;

for i=1:70
    disp(minCCFMale(i));
    if(minCCFMale(i)>maxCCFFemale(i))
        CCFdiff(i)=minCCFMale(i)-maxCCFFemale(i);
    elseif(minCCFFemale(i)>maxCCFMale(i))
        CCFdiff(i)=minCCFFemale(i)-maxCCFMale(i);
    else
        CCFdiff(i) = -1;
    end
end
[CCFval,CCFindex]=sort(CCFdiff, 'descend');
CCFnodes = CCFindex(CCFval>0);
for i=1:size(CCFnodes,2)
        newPval(i)=bootstrap_ccf_node('C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\', CCFnodes(i),MeanCCFDiff(CCFnodes(i)));
end
considerCCF=CCFnodes(newPval<0.05);
%plotting clustering coeff for visualization
figure(2);
hold on
errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');
errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales,errorFemales,'.k', 'color', 'red');
hold off
xlabel('Brain region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Mean of clustering coefficient for each brain region',  'FontSize',16);
legend('Females', 'Males');

MeanEBCFemales = mean(OF);
MeanEBCDiff=MeanEBCMales-MeanEBCFemales;
errorEBCFemales = std(OF)/ sqrt(size(OF,1));
edgelow=1;
edgehigh=4900;

figure(3);
writetoPAJ(MeanZM,'MeanSex0',1);
writetoPAJ(MeanZF,'MeanSex1',1);
heatmap((MeanZF - MeanZM), 1:70, 1:70);
hold on;
xlabel('Brain region', 'FontSize',14);
ylabel('Brain region', 'FontSize',14);
title('Differences in mean connectivity for each edge between the connectomes', 'FontSize',16);
colorbar;
%colormap bone;
hold off;

% significance computation for EBC
errorBarEBCMale=errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
maxEBCMale = MeanEBCMales + errorBarEBCMale.LData;
minEBCMale = MeanEBCMales - errorBarEBCMale.UData;

errorBarEBCFemale=errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red');
maxEBCFemale = MeanEBCFemales + errorBarEBCFemale.LData;
minEBCFemale = MeanEBCFemales - errorBarEBCFemale.UData;


for i=1:size(minEBCFemale,2)
    disp(minEBCMale(i));
    if(minEBCMale(i)>maxEBCFemale(i))
        EBCdiff(i)=minEBCMale(i)-maxEBCFemale(i);
    elseif(minEBCFemale(i)>maxEBCMale(i))
        EBCdiff(i)=minEBCFemale(i)-maxEBCMale(i);
    else
        EBCdiff(i) = -1;
    end
end

[EBCval,EBCindex]=sort(EBCdiff, 'descend');
EBCedges = EBCindex(EBCval>0);
for i=1:size(EBCedges,2)
        EBCPval(i)=bootstrap_connectivity_edge('C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\', EBCedges(i),MeanEBCDiff(EBCedges(i)));
end
considerEBC=EBCedges(EBCPval<0.05);

%PPF bootstrapping + high low visualization..
MeanPPFNodeFemales = mean(OF_PPF);
errorPPFFemales = std(OF_PPF)/sqrt(size(OF_PPF,1));
MeanPPFDiff=MeanPPFNodeFemales-MeanPPFNodeMales;

errorBarPPFMale = errorbar(1:prod(size(MeanPPFNodeMales)), MeanPPFNodeMales,errorPPFMales,'.k', 'color', 'blue');
maxPPFMale = MeanPPFNodeMales + errorBarPPFMale.LData;
minPPFMale = MeanPPFNodeMales - errorBarPPFMale.UData;

errorPPFFemales = std(ccfFemale)/sqrt(size(ccfFemale,1));
errorBarPPFFemale = errorbar(1:prod(size(MeanPPFNodeFemales)), MeanPPFNodeFemales,errorPPFFemales,'.k', 'color', 'red');
maxPPFFemale = MeanPPFNodeFemales + errorBarPPFFemale.LData;
minPPFFemale = MeanPPFNodeFemales - errorBarPPFFemale.UData;

for i=1:size(minPPFFemale,2)
    disp(minPPFMale(i));
    if(minPPFMale(i)>maxPPFFemale(i))
        PPFdiff(i)=minPPFMale(i)-maxPPFFemale(i);
    elseif(minPPFFemale(i)>maxPPFMale(i))
        PPFdiff(i)=minPPFFemale(i)-maxPPFMale(i);
    else
        PPFdiff(i) = -1;
    end
end

[PPFval,PPFindex]=sort(PPFdiff, 'descend');
PPFnodes = PPFindex(PPFval>0);
% PPFval=PPFval(PPFval>0);
for i=1:size(PPFnodes,2)
        PPFPval(i)=bootstrap_partcf('C:\Users\minug\Pictures\jagatsastry-brain-analysis-for-gender-classification-b4e99bb2b4c5\datasets\set1\normalized\', PPFnodes(i),MeanPPFDiff(PPFnodes(i)));
end
considerPPF=PPFnodes(PPFPval<0.05);

%plotting for EBC
figure(4);
hold on;
errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red')
xlabel('Inter-region Edge',  'FontSize',14);
ylabel('Edge Betweenness',  'FontSize',14);
title('Mean of edge betweenness for each inter-region edge',  'FontSize',14);
legend('Males', 'Females');
hold off;
toc;
