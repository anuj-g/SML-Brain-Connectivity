function  ccf_node55_hist(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/males/');
mfiles = dir(strcat(maledir, '*.mat'));
for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));    %Calculate the mean clustering coeff in that connectome
    ccf = clustering_coef_wd(M.fibergraph);
    %Store that mean in an array
    GM(k) = ccf(55);
end
figure(1)
hist(GM)

xlabel('Weighted Clustering Coefficient',  'FontSize',14);
ylabel('Count',  'FontSize',14);
title('Histogram of clustering coefficient of node 55 across male subjects',  'FontSize',16);


% Do the exact same thing for other sex as well
femdir = strcat(normalizeddir, '/females/');
femfiles = dir(strcat(femdir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femdir, femfiles(k).name));
    %Calculate the mean clustering coeff for that connectome
    ccf = clustering_coef_wd(M.fibergraph);
    %Store that mean in an array
    GF(k) = ccf(55);
end
figure(2)
hist(GF)
xlabel('Weighted Clustering Coefficient',  'FontSize',14);
ylabel('Count',  'FontSize',14);
title('Histogram of clustering coefficient of node 55 across female subjects',  'FontSize',16);

end