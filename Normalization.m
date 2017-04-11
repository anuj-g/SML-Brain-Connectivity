clear all;
clc;
tic;
highiqDatasetDir='datasets\fsiq\high\';
lowiqDatasetDir='datasets\fsiq\low\';

%CCF
[MeanCCFNodeMales,errorMales,MeanCCFNodeFemales,errorFemales] = ccf_normalized_errorbars('datasets\fsiq\');  

errorBarCCFMale = errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');
maxCCFMale = MeanCCFNodeMales + errorBarCCFMale.YPositiveDelta;
minCCFMale = MeanCCFNodeMales - errorBarCCFMale.YNegativeDelta;

errorBarCCFFemale = errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales,errorFemales,'.k', 'color', 'red');
maxCCFFemale = MeanCCFNodeFemales + errorBarCCFFemale.YPositiveDelta;
minCCFFemale = MeanCCFNodeFemales - errorBarCCFFemale.YPositiveDelta;

MeanCCFDiff=MeanCCFNodeFemales-MeanCCFNodeMales;

for i=1:70
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
        newPval(i)=bootstrap_ccf_node('datasets\fsiq\', CCFnodes(i),MeanCCFDiff(CCFnodes(i)));
end
considerCCF=CCFnodes(newPval<0.05);

%EBC
edgelow=1;
edgehigh=4900;

[MeanEBCMales, errorEBCMales, MeanEBCFemales, errorEBCFemales] = edge_betweenness_normalized_errorbars_final('datasets/fsiq/', edgelow, edgehigh);
errorBarEBCMale=errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
maxEBCMale = MeanEBCMales + errorBarEBCMale.YPositiveDelta;
minEBCMale = MeanEBCMales - errorBarEBCMale.YNegativeDelta;

errorBarEBCFemale=errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red');
maxEBCFemale = MeanEBCFemales + errorBarEBCFemale.YPositiveDelta;
minEBCFemale = MeanEBCFemales - errorBarEBCFemale.YNegativeDelta;

MeanEBCDiff=MeanEBCMales-MeanEBCFemales;

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
for i=1:10
        EBCPval(i)=bootstrap_ebc_edge('datasets\fsiq\', EBCedges(i),MeanEBCDiff(EBCedges(i)));
end
considerEBC=EBCedges(EBCPval<0.05);

%PPF
[MeanPPFNodeMales, errorPPFMales, MeanPPFNodeFemales, errorPPFFemales] = partcf_normalized_errorbars('datasets/fsiq/');
errorBarPPFMale = errorbar(1:prod(size(MeanPPFNodeMales)), MeanPPFNodeMales,errorPPFMales,'.k', 'color', 'blue');
maxPPFMale = MeanPPFNodeMales + errorBarPPFMale.YPositiveDelta;
minPPFMale = MeanPPFNodeMales - errorBarPPFMale.YNegativeDelta;
errorBarPPFFemale = errorbar(1:prod(size(MeanPPFNodeFemales)), MeanPPFNodeFemales,errorPPFFemales,'.k', 'color', 'red');
maxPPFFemale = MeanPPFNodeFemales + errorBarPPFFemale.YPositiveDelta;
minPPFFemale = MeanPPFNodeFemales - errorBarPPFFemale.YNegativeDelta;

for i=1:size(minPPFFemale,2)
    if(minPPFMale(i)>maxPPFFemale(i))
        PPFdiff(i)=minPPFMale(i)-maxPPFFemale(i);
    elseif(minPPFFemale(i)>maxPPFMale(i))
        PPFdiff(i)=minPPFFemale(i)-maxPPFMale(i);
    else
        PPFdiff(i) = -1;
    end
end
MeanPPFDiff=MeanPPFNodeFemales-MeanPPFNodeMales;
[PPFval,PPFindex]=sort(PPFdiff, 'descend');
PPFnodes = PPFindex(PPFval>0);

for i=1:size(PPFnodes,2)
        PPFPval(i)=bootstrap_partcf('datasets\fsiq\', PPFnodes(i),MeanPPFDiff(PPFnodes(i)));
end
considerPPF=PPFnodes(PPFPval<0.05);

%plotting clustering coeff for visualization
figure(1);
hold on
errorbar(1:prod(size(MeanCCFNodeMales)), MeanCCFNodeMales,errorMales,'.k', 'color', 'blue');
errorbar(1:prod(size(MeanCCFNodeFemales)), MeanCCFNodeFemales,errorFemales,'.k', 'color', 'red');
xlabel('Brain region',  'FontSize',14);
ylabel('Clustering Coefficient',  'FontSize',14);
title('Mean of clustering coefficient for each brain region',  'FontSize',16);
legend('Females', 'Males');
hold off

%plotting for EBC
figure(2);
hold on;
errorbar(edgelow:edgehigh, MeanEBCMales(edgelow:edgehigh), errorEBCMales(edgelow:edgehigh), '.k','color','blue');
errorbar(edgelow:edgehigh, MeanEBCFemales(edgelow:edgehigh), errorEBCFemales(edgelow:edgehigh),'.k','color', 'red')
xlabel('Inter-region Edge',  'FontSize',14);
ylabel('Edge Betweenness',  'FontSize',14);
title('Mean of edge betweenness for each inter-region edge',  'FontSize',14);
legend('Males', 'Females');
hold off;

%plotting for ppf
figure(3);
hold on;
errorbar(1:prod(size(MeanPPFNodeMales)), MeanPPFNodeMales,errorPPFMales,'.k', 'color', 'blue');
errorbar(1:prod(size(MeanPPFNodeFemales)), MeanPPFNodeFemales,errorPPFFemales,'.k', 'color', 'red');
xlabel('Node',  'FontSize',14);
ylabel('PPF',  'FontSize',14);
title('PPF',  'FontSize',14);
legend('High', 'Low');
hold off;

toc;