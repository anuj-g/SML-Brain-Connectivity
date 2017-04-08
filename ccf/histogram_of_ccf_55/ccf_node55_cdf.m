function  ccf_node55_cdf(normalizeddir)

%Process Sex 0
maledir = strcat(normalizeddir, '/males/');
mfiles = dir(strcat(maledir, '*.mat'));
for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));    %Calculate the mean clustering coeff in that connectome
    ccf = clustering_coef_wd(M.fibergraph);
    %Store that mean in an array
    GM(k) = ccf(55);
end


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

figure(1);
hm = cdfplot(GM);
set(hm, 'color', 'b');

hold on;

hf = cdfplot(GF);
set(hf, 'color', 'r');

xlabel('Clustering Coefficient of Node 55',  'FontSize',14);
ylabel('F(x)',  'FontSize',14);
title('Empirical CDF of Node 55',  'FontSize',16);
legend('Male', 'Female');
hold off;

end