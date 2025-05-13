function plot_means_mu(main_folder)

% PLOT_MEANS_MU  Reads .mat files from 10 subfolders, computes mean of the
% third row across each file, and plots the results for mu = 0.5:0.05:0.95.
%
%   plot_means_mu(main_folder)
%       main_folder: path to the parent directory containing 10 subfolders.

    % Check input folder
    if ~isfolder(main_folder)
        error('Directory "%s" does not exist.', main_folder);
    end

    % Define x-axis values and mu list
    x = 0.01:0.01:0.25;                % 25 points
    mu_values = 0.5:0.05:0.95;         % 10 mu values
    legend_entries = arrayfun(@(m) sprintf('\\mu = %.2f', m), mu_values, 'UniformOutput', false);
    
    % Get subfolder list
    D = dir(main_folder);
    subdirs = D([D.isdir] & ~ismember({D.name}, {'.','..'}));
    if numel(subdirs) ~= numel(mu_values)
        warning('Expected %d subfolders but found %d.', numel(mu_values), numel(subdirs));
    end

    figure;
    hold on;
    % Loop over subfolders
    for i = 1:min(numel(subdirs), numel(mu_values))
        sub_path = fullfile(main_folder, subdirs(i).name);
        mat_files = dir(fullfile(sub_path, '*.mat'));
        % Preallocate
        means_i = nan(1, numel(mat_files));
        
        for j = 1:numel(mat_files)
            % Load .mat file (assume it contains one cell variable)
            data = load(fullfile(sub_path, mat_files(j).name));
            fn = fieldnames(data);
            C = data.(fn{1});             % the cell array 3x200
            % Extract third row cells and convert to numeric
            row3 = C(3, :);
            v = cell2mat(row3);          % 1x200 numeric
            means_i(j) = mean(v);
        end
        % Plot this mu curve
        plot(x, means_i, 'LineWidth', 1.5);
    end
    
    % Finalize plot
    hold off;
    xlabel('x');
    ylabel('ITERAZIONI MEDIE');
    title('ITERAZIONI MEDIE RISPETTO A \mu');
    legend(legend_entries, 'Location', 'best');
    grid on;
end