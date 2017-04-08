%Process Sex:0
clearvars;
format long;
cd /home/jagat/brain_analysis/datasets/set1/normalized/males
files = dir('*.mat');

for k = 1:numel(files)
    M = load(files(k).name);
    ccfm = clustering_coef_wd(M.fibergraph);
    reccfm = reshape(ccfm, 1, prod(size(ccfm)));
    GM(k,:) = [0 reccfm];
end
k = k + 1
size(GM)

cd /home/jagat/brain_analysis/datasets/set1/normalized/females
files = dir('*.mat');
for j = 1:numel(files)
    M = load(files(j).name);
    ccff = clustering_coef_wd(M.fibergraph);
    reccff = reshape(ccff, 1, prod(size(ccff)));
    GM(k,:) = [1 reccff];
    k = k + 1;
end
size(GM)
cd /home/jagat/brain_analysis/
dlmwrite('CCFCombined_SVM.csv',GM,'precision','%10.5f');
