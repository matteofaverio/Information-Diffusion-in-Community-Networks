function plot_H_BC(H, BC)
%PLOT_H_BC_CI  Visualizza la nuvola (H, BC, CI) in uno scatter 3D.
%
%   plot_H_BC_CI(H, BC, CI)
%
%   INPUT
%     H  : vettore entropie normalizzate   (N×1)
%     BC : vettore Bimodality Coefficient  (N×1)
%     CI : vettore Concentration Index     (N×1) – opzionale ma consigliato
%
%   Se ci sono NaN, vengono ignorati automaticamente.

    if nargin < 3 || isempty(CI)
        CI = NaN(size(H));    % permetti chiamata con due soli vettori
    end

    % Pulizia: rimuovi righe con NaN simultanei
    valid = ~isnan(H) & ~isnan(BC);
    H  = H(valid);  BC = BC(valid);  CI = CI(valid);

    % Normalizza CI per usarlo come colore (se presente)
    if all(isnan(CI))
        cdata = 0.5 * ones(size(H));   % grigio fisso
    else
        cdata = (CI - min(CI)) ./ (max(CI) - min(CI) + eps);
    end

    % Dimensione marker legata a BC (scala log per equilibrio visivo)
    sz = 30 + 120 * (log(1 + BC - min(BC)) ./ log(1 + max(BC) - min(BC)));

    % Grafico
    figure('Color','w');
    scatter(H, BC, sz, cdata, 'filled', 'MarkerFaceAlpha', 0.65);
    % colormap('turbo');        % o 'parula' se preferisci
    if all(isnan(CI))
        cb.Label.String = 'No CI (fixed colour)';
    else
        cb.Label.String = 'Concentration Index (scaled 0-1)';
    end

    grid on
    
    xlabel('Shannon entropy H (0 = consenso, 1 = frammentazione)','interpreter','latex');
    ylabel('Bimodality coefficient BC','interpreter','latex');
    title('Distribuzione congiunta di H, BC','interpreter','latex');

    view(2)
    % view(45, 25);             % angolo isometrico gradevole
end
