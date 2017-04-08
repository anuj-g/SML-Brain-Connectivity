%Process Sex:0
clearvars;
format long;
cd /home/jagat/brain_analysis/datasets/set1/normalized/males
files = dir('*.mat');
for k = 1:numel(files)
    M = load(files(k).name);
    CI = modularity_dir(M.fibergraph);
    P =  participation_coef(M.fibergraph, CI);
    recppf = reshape(P, 1, prod(size(P)));
    GM(k,:) = [0 recppf];
end
k = k + 1
size(GM)
cd /home/jagat/brain_analysis/datasets/set1/normalized/females
files = dir('*.mat');
for j = 1:numel(files)
    M = load(files(j).name);
    CI = modularity_dir(M.fibergraph);
    P = participation_coef(M.fibergraph, CI);
    reppf = reshape(P, 1, prod(size(P)));
    GM(k,:) = [1 reppf];
    k = k + 1;
end
size(GM)
cd /home/jagat/brain_analysis
dlmwrite('PPFCombined_SVM.csv',GM,'precision','%10.5f');
