%Process Sex:0
function  normalizematrix(datasetdir)


outdir = strcat(datasetdir, '/normalized');
mkdir(outdir);

mkdir(strcat(outdir, '/merged'));

mfiles = dir(strcat(datasetdir, '/males/*.mat'));
mkdir(strcat(outdir, '/males'));
  
for k = 1:numel(mfiles)
    M = load(strcat(strcat(datasetdir, '/males/'), mfiles(k).name));
    % M is a sparse matrix and it is symmetric
    % The bottom half is 0 becuase it is sparse, we cant pass this 
    % directly for computing weights as the 0's will be over counted
    % Hnce we convert it to full matrix , make it symmetric and pass that
    % to get the weight distribution
    F = full(M.fibergraph);
    % Adding F to F transpose
    S = F + F';
    isequal(S, S');
    fibergraph = normalizematrixfunction(S);
    normalizedfilename = sprintf(strcat(outdir, '/males/normalized_%s'), mfiles(k).name);
    save(normalizedfilename, 'fibergraph');
    
    mergedfilename = sprintf(strcat(outdir, '/merged/normalized_%s'), mfiles(k).name);
    save(mergedfilename, 'fibergraph');

end

femfiles = dir(strcat(datasetdir, '/females/*.mat'));

mkdir(strcat(outdir, '/females'));

for k = 1:numel(femfiles)
    M = load(strcat(strcat(datasetdir, '/females/'), femfiles(k).name));
    % M is a sparse matrix and it is symmetric
    % The bottom half is 0 becuase it is sparse, we cant pass this 
    % directly for computing weights as the 0's will be over counted
    % Hnce we convert it to full matrix , make it symmetric and pass that
    % to get the weight distribution
    F = full(M.fibergraph);
    % Adding F to F transpose
    S = F + F';
    isequal(S, S');
    fibergraph = normalizematrixfunction(S);    
    normalizedfilename = sprintf(strcat(outdir, '/females/normalized_%s'), femfiles(k).name);
    save(normalizedfilename, 'fibergraph');
    
    mergedfilename = sprintf(strcat(outdir, '/merged/normalized_%s'), femfiles(k).name);
    save(mergedfilename, 'fibergraph');

end

