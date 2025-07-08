function [H, BC, CI] = community_plots_artigianale(folderName)

% caricamento dati
mus   = 0.5:-0.05:0.05;      % 10 valori di mu
epsilon = linspace(0.01,0.25,25);
nMus       = numel(mus);

% --- Elenco e ordinamento file ---
files = dir(fullfile(folderName,'*.mat'));
%files = files(~strncmp({files.name}, '._', 2));
if numel(files)~=nMus
    error('Trovati %d file, ma ci aspettavamo %d (numero di mu)', ...
          numel(files), nMus);
end
names = sort({files.name});
fprintf('--- files caricati ---\n')

for i = 1:nMus
    % fullfile(folderName, names{i})
    S = load(fullfile(folderName, names{i}));
    fn = fieldnames(S);
    results = S.(fn{1});
    indexes = fieldnames(results);

    for j = 1:numel(epsilon)
        h(:,j) = results.(indexes{1}){j};
        ci(:,j) = results.(indexes{2}){j};
        bc(:,j) = results.(indexes{3}){j};
    end
    
    H(i,:) = mean(h,1);
    CI(i,:) = mean(ci,1);
    BC(i,:) = mean(bc,1);
end

%% ----- Script: plot3D_surface_professionale_fixed.m -----
figure('Color','w')%,'Position',[150 150 900 600]);

% Griglia 2D di punti (epsilon, mu)
[EPS, MU] = meshgrid(epsilon, mus);

% Surface plot
% hSurf = surf(EPS, MU, H, ...
%     'EdgeColor','k', ...
%     'FaceColor','interp');

% Surface plot con bordi neri
hSurf = surf( EPS, MU, H, ...
    'FaceColor','interp', ...   % colore interpolato sui vertici
    'EdgeColor','k', ...        % bordi neri lungo la griglia dati
    'LineWidth',0.5 );          % spessore linea sottile


% Miglioramenti grafici
colormap(jet);               % scala colori vivida
c = colorbar;                  % barra dei colori

view(45,30);

% Label e stile assi
xlabel('$\varepsilon$','Interpreter','latex','FontSize',14);
ylabel('$\mu$','Interpreter','latex','FontSize',14);
zlabel('Shannon Entropy','Interpreter','latex','FontSize',14);
title('Shannon Entropy in funzione di $\varepsilon$ e $\mu$',...
    'Interpreter','latex','FontSize',16);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1);
grid on; box on;
axis tight;

%% ----- Script: plot3D_surface_professionale_fixed.m -----
figure('Color','w')%,'Position',[150 150 900 600]);

% Griglia 2D di punti (epsilon, mu)
[EPS, MU] = meshgrid(epsilon, mus);

% Surface plot
% hSurf = surf(EPS, MU, H, ...
%     'EdgeColor','k', ...
%     'FaceColor','interp');

% Surface plot con bordi neri
hSurf = surf( EPS, MU, BC, ...
    'FaceColor','interp', ...   % colore interpolato sui vertici
    'EdgeColor','k', ...        % bordi neri lungo la griglia dati
    'LineWidth',0.5 );          % spessore linea sottile


% Miglioramenti grafici
colormap(jet);               % scala colori vivida
c = colorbar;                  % barra dei colori

view(45,30);

% Label e stile assi
xlabel('$\varepsilon$','Interpreter','latex','FontSize',14);
ylabel('$\mu$','Interpreter','latex','FontSize',14);
zlabel('Bimodality Coefficient','Interpreter','latex','FontSize',14);
title('Bimodality Coefficient in funzione di $\varepsilon$ e $\mu$',...
    'Interpreter','latex','FontSize',16);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1);
grid on; box on;
axis tight;

%% ----- Script: plot3D_surface_professionale_fixed.m -----
figure('Color','w')%,'Position',[150 150 900 600]);

% Griglia 2D di punti (epsilon, mu)
[EPS, MU] = meshgrid(epsilon, mus);

% Surface plot
% hSurf = surf(EPS, MU, H, ...
%     'EdgeColor','k', ...
%     'FaceColor','interp');

% Surface plot con bordi neri
hSurf = surf( EPS, MU, CI, ...
    'FaceColor','interp', ...   % colore interpolato sui vertici
    'EdgeColor','k', ...        % bordi neri lungo la griglia dati
    'LineWidth',0.5 );          % spessore linea sottile


% Miglioramenti grafici
colormap(jet);               % scala colori vivida
c = colorbar;                  % barra dei colori

view(45,30);

% Label e stile assi
xlabel('$\varepsilon$','Interpreter','latex','FontSize',14);
ylabel('$\mu$','Interpreter','latex','FontSize',14);
zlabel('Centrality Index','Interpreter','latex','FontSize',14);
title('Centrality Index in funzione di $\varepsilon$ e $\mu$',...
    'Interpreter','latex','FontSize',16);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1);
grid on; box on;
axis tight;
end

