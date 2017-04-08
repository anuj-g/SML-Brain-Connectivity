%Process Sex:0
clearvars;
format long;
cd /home/jagat/brain_analysis/datasets/set1/normalized/males
files = dir('*.mat');

for k = 1:numel(files)
    M = load(files(k).name);
    D = M.fibergraph;
    E = find(D);
    D(E) = 1./D(E);
    EBC = edge_betweenness_wei(D);
    RS  = reshape(EBC, 1, prod(size(EBC)));
    GM(k,:) = [0 RS];
end
k = k + 1
size(GM)

cd /home/jagat/brain_analysis/datasets/set1/normalized/females
files = dir('*.mat');
for j = 1:numel(files)
    M = load(files(j).name);
    D = M.fibergraph;
    E = find(D);
    D(E) = 1./D(E);
    EBC = edge_betweenness_wei(D);
    RS  = reshape(EBC, 1, prod(size(EBC)));
    GM(k,:) = [1 RS];
    k = k + 1;
end
size(GM)
cd /home/jagat/brain_analysis/
dlmwrite('EBCCombined_SVM.csv',GM,'precision','%10.5f');
