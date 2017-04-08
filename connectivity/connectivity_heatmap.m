function connectivity_heatmap(normdir)

%Process Sex:0
% Clear variables
%Main Logic: http://www.mathworks.com/matlabcentral/answers/34395

maledir = strcat(normdir, '/males/');
mfiles = dir(strcat(maledir, '*.mat'));
ZM = [];
%Initalize an empty matrix
for k = 1:numel(mfiles)
    M = load(strcat(maledir, mfiles(k).name));
    % Add this matrix to the 3rd dimension of the running matrix
    ZM = cat(3, ZM, M.fibergraph);
end
% Calculate mean of each element of matrix
MeanZM = mean(ZM, 3);

% Process the next sex
femdir = strcat(normdir, '/females/');
ZF = [];
femfiles = dir(strcat(femdir, '*.mat'));
for k = 1:numel(femfiles)
    F = load(strcat(femdir, femfiles(k).name));
    %Append each connectivity matrix to the running matrix
    ZF = cat(3, ZF, F.fibergraph);
end
% Find the mean for each element
MeanZF = mean(ZF, 3);

writetoPAJ(MeanZM,'MeanSex0',1);
writetoPAJ(MeanZF,'MeanSex1',1)

% Draw a heat map between of the difference between the Means
heatmap((MeanZF - MeanZM), 1:70, 1:70);
hold on;
xlabel('Brain region', 'FontSize',14);
ylabel('Brain region', 'FontSize',14);
title('Differences in mean connectivity for each edge between the connectomes', 'FontSize',16);
colorbar;
%colormap bone;
hold off;
end