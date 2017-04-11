clc;
clear all;
directory='datasets/fsiq/merged/';
dataPath=strcat(directory,'*.mat');
numofFiles=dir(dataPath);

highiqDatasetDir='datasets\fsiq\high\';
lowiqDatasetDir='datasets\fsiq\low\';
dataPath=strcat(directory,'*.mat');
mergedFiles = dir('datasets/fsiq/merged/');
highiqfiles = dir('datasets/fsiq/high');
lowiqfiles = dir('datasets/fsiq/low');
highfiles = zeros(65);
highindex = 1;
lowindex = 1;
strHigh = strings(65);
strLow = strings(51);
labels = zeros(114,1);
for i = 1:63
    strHigh(i) = highiqfiles(i+2).name;
end
for i = 1:51
    strLow(i) = lowiqfiles(i+2).name;
end
labelIndex = 1;
for i=1:116
    disp(mergedFiles(i,1).name);
    if ismember( mergedFiles(i,1).name,strHigh)    
        labels(labelIndex,1) = 1;
        labelIndex = labelIndex + 1;
    elseif ismember( mergedFiles(i,1).name,strLow)
        labels(labelIndex,1) = 0;
        labelIndex = labelIndex + 1;
    end
end
ZM=[];
for k = 1:numel(numofFiles)
    file=strcat(strcat(directory, numofFiles(k).name));
    M = load(file);
    F = full(M.fibergraph);
    ccf(k,:) = clustering_coef_wd(F);     
    CI = modularity_dir(F);
    %Calculate the Clustering Coeff of each node   
    PPF(k,:) = participation_coef(F, CI);    
    %Edgebetweenness connectivity
    E= find(F);
    F(E) = 1./F(E);    
    EBC = edge_betweenness_wei(F);
    EdgeBetween(k,:) = reshape(EBC, prod(size(EBC)),1);        
    %connectivity
    ZM = cat(3, ZM, F);      
end