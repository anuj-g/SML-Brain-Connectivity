function efficiency_normalized_errorbars(normdir)
%Process Sex:0

%Clear variables

maledir = strcat(normdir, '/males/');
mfiles = dir(strcat(maledir, '*.mat'));
for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    %Calculate the effeciency of each node(local effeciency)
    lcm = efficiency_wei(M.fibergraph,1);
    %Store the local eff of each node as a row in a matrix
    OM(k,:) = lcm;
end

% Find the mean local eff of each node (mean of each column)
MeanLCENodeMales = mean(OM);
% Compute the standard error as well for each node
errorMales = std(OM)/sqrt(size(OM,1))
% Plot the mean local eff of each node with error bars in blue
errorbar(1:prod(size(MeanLCENodeMales)), MeanLCENodeMales,errorMales,'.k', 'color', 'blue');
% Use the same figure for further plots
hold on;

% Process the next sex
femaledir = strcat(normdir, '/females/');
femfiles = dir(strcat(femaledir, '*.mat'));
for k = 1:numel(femfiles)
    M = load(strcat(femaledir, femfiles(k).name));
    %Compute the local efficiency for each node
    lcf = efficiency_wei(M.fibergraph,1);
    OF(k,:) = lcf;
end

%Find the mean local eff for each node(mean for each column)
MeanLCFNodeFemales = mean(OF);

%Find the standard error as well
errorFemales = std(OF)/sqrt(size(OF,1))

%Plot the mean local eff of with error bars in red
errorbar(1:prod(size(MeanLCFNodeFemales)), MeanLCFNodeFemales, errorFemales,'.k', 'color', 'red');
xlabel('Brain region',  'FontSize',14);
ylabel('Local Efficiency',  'FontSize',14);
title('Mean of local efficiency for each brain region',  'FontSize',16);
legend('Males', 'Females');

hold off;
clearvars;
end