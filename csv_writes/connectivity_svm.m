%Process Sex:0
clearvars;
format long;
cd /home/jagat/brain_analysis/datasets/set1/normalized/males
files = dir('*.mat');

for k = 1:numel(files)
    M = load(files(k).name);
    R = reshape(M.fibergraph', 1, prod(size(M.fibergraph)));
    GM(k,:) = [0 R];
end
k = k + 1
size(GM)
cd /home/jagat/brain_analysis/datasets/set1/normalized/females
files = dir('*.mat');
for j = 1:numel(files)
    M = load(files(j).name);
    R = reshape(M.fibergraph', 1, prod(size(M.fibergraph)));
    GM(k,:) = [1 R];
    k = k + 1;
end
size(GM)
cd /home/jagat/brain_analysis
dlmwrite('SexCombined.csv',GM,'precision','%10.5f');
%dlmwrite('Sex1.csv',GF,'precision','%10.5f');
