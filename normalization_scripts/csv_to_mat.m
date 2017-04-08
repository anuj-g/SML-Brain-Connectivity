function csv_to_mat(csvdir)

maledir = strcat(csvdir, '/males/');
mfiles = dir(strcat(maledir, '*.csv'));
for k = 1:numel(mfiles)
    fname = strcat(maledir, mfiles(k).name);
    fibergraph = dlmread(fname,'',[1,1,70,70]);
    normalizedfilename = strcat(fname, '.mat');
    save(normalizedfilename, 'fibergraph');
end

femaledir = strcat(csvdir, '/females/');
femfiles = dir(strcat(femaledir, '*.csv'));
for k = 1:numel(femfiles)
    fname = strcat(femaledir, femfiles(k).name);
    fibergraph = dlmread(fname,'',[1,1,70,70]);
    normalizedfilename = strcat(fname, '.mat');
    save(normalizedfilename, 'fibergraph');
end

end