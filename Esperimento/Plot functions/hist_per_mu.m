function hist_per_mu(parentFolder, idx)

%PROCESS_MAT_CELLS Load the i-th .mat file from each subfolder and extract its first variable (a cell).
%   parentFolder: path to the main folder containing subfolders
%   idx: index of the file to load in each subfolder
%   results: cell array of outputs from each subfolder

% Get list of subfolders
d = dir(parentFolder);
isub = [d(:).isdir];
names = {d(isub).name};
% Remove . and ..
names = names(~ismember(names, {'.', '..'}));
nSub = numel(names);

for k = 1:nSub

    subPath = fullfile(parentFolder, names{k});
    % List .mat files in subfolder
    matFiles = dir(fullfile(subPath, '*.mat'));
    if numel(matFiles) < idx
        warning('Subfolder %s has fewer than %d .mat files. Skipping.', names{k}, idx);
        continue;
    end
    % Select the i-th .mat file
    fileName = matFiles(idx).name;
    filePath = fullfile(subPath, fileName);
    % Load the .mat file
    S = load(filePath);
    % Extract the first variable (assumed to be a cell)
    varNames = fieldnames(S);
    firstVar = S.(varNames{1});
    if ~iscell(firstVar)
        error('First variable in %s is not a cell.', filePath);
    end
    C = firstVar;

    N = 1e4;
    op_end = zeros(N,1);

    for i = 1:200
        x = C{2,i};
        x = double(x)/1000;
        op_end = op_end + sort(x);
    end

    op_end = op_end/200;
    
    figure;
    histogram(op_end,51);
    ylim([0 10000]); 
    
end
end

