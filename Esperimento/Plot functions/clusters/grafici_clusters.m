%%

% APRIRE CLUSTERS.mat
% dentro ci sono 5 variabili, tutte matrici 10x25 (10 valori di mu, 25 di 
% epsilon):
% - n_clusters : numero medio di clusters
% - var_clusters : varianza del numero di clusters
% - outliers : numero medio di nodi fuori dai clusters
% - meanL e meanH: media dei valori attorno a cui gli outlier si
%                  distribuiscono. Sono due valori, uno tra [0,0.5] e uno 
%                  tra [0.5,1], perchè sono le medie per i valori nei due 
%                  sotto intervalli, altrimenti il valore medio era sempre 
%                  0.5 e non era significativo. Mostra come gli outlier si
%                  allontanano sempre di più all'aumentare di epsilon

% selezionare M come matrice che si vuole plottare, ad esempio n_clusters

M1 = meanH; %
M2 = meanL;

%%  GRAFICO 2D PER CLUSTERS

mu_values = linspace(0.5,0.05,10);
epsilon = linspace(0.01,0.25,25);
% ----- Script: plot2D_mat_professionale_fixed_legend.m -----
figure('Color','w','Position',[100 100 800 500]);
hold on;

% Prepara la mappa colori
nMU   = size(M,1);
cmap  = parula(nMU);
xlim([0, 0.25]);

% Disegna ciascuna curva con legend in math-mode
for i = 1:nMU
    plot(epsilon, M1(i,:), ...
        'LineWidth',1.8, ...
        'Color',cmap(i,:), ...
        'DisplayName',sprintf('$\\mu = %.2f$' , mu_values(i)));
end

for i = 1:nMU
    plot(epsilon, M2(i,:), ...
        'LineWidth',1.8, ...
        'Color',cmap(i,:), ...
        'DisplayName',sprintf('$\\mu = %.2f$' , mu_values(i)),'HandleVisibility','off');
end

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',18, ...
    'LineWidth',1);
grid on; box on;

% Legend esterna con interpreter LaTeX
lh = legend('Location','northeastoutside');
set(lh, ...
    'Interpreter','latex', ...
    'FontSize',16);
xlabel('$\varepsilon$','Interpreter','latex','FontSize',24);
ylabel('opinione media degli \emph{outlier}','Interpreter','latex','FontSize',20);

% Colorbar che rimappa i mu_values
%{
c = colorbar('eastoutside');
colormap(parula);
c.Ticks      = linspace(0,1,nMU);
c.TickLabels = arrayfun(@(x) sprintf('%.2f',x), mu_values,'UniformOutput',false);
c.Label.String      = '$\mu$';
c.Label.Interpreter = 'latex';
c.Label.FontSize    = 14;
%}


%%

mu_values = linspace(0.5,0.05,10);
epsilon = linspace(0.01,0.25,25);
% ----- Script: plot2D_mat_professionale_fixed_legend.m -----
figure('Color','w','Position',[100 100 800 500]);
hold on;

% Prepara la mappa colori
nMU   = size(M,1);
cmap  = cool(nMU);

% Disegna ciascuna curva con legend in math-mode
for i = 1:nMU
    plot(epsilon, M(i,:), ...
        'LineWidth',1.8, ...
        'Color',cmap(i,:), ...
        'DisplayName',sprintf('$\\mu = %.2f$' , mu_values(i)));
end

% Styling accademico

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',18, ...
    'LineWidth',1);
grid on; box on;
xlim([0, 0.25]);

% Legend esterna con interpreter LaTeX
lh = legend('Location','northeast');
set(lh, ...
    'Interpreter','latex', ...
    'FontSize',16);
xlabel('$\varepsilon$','Interpreter','latex','FontSize',24);
ylabel('varianza del \# cluster','Interpreter','latex','FontSize',20);

%%

% ----- Script: plot3D_surface_con_grid.m -----
figure('Color','w','Position',[150 150 900 600]);
hold on;

% Griglia 2D di punti (epsilon, mu)
[EPS, MU] = meshgrid(epsilon, mu_values);

% Surface plot con bordi neri
hSurf = surf( EPS, MU, M, ...
    'FaceColor','interp', ...   % colore interpolato sui vertici
    'EdgeColor','k', ...        % bordi neri lungo la griglia dati
    'LineWidth',0.5 );          % spessore linea sottile


% Colormap e colorbar
colormap(parula);
c = colorbar;                  


view(45,30);

% Etichette con LaTeX
xlabel('$\varepsilon$','Interpreter','latex','FontSize',14);
ylabel('$\mu$','Interpreter','latex','FontSize',14);

% Styling assi
set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',18, ...
    'LineWidth',1);
grid on; box on;  % mantiene anche la griglia di riferimento dell’asse
axis tight;

xlabel('$\varepsilon$','Interpreter','latex','FontSize',24);

