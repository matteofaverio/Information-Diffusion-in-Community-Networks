function finalVec = mean_comm_var(processedFolder)
% aggregateProcessedResults – Calcola il valore medio complessivo per ciascun file nella cartella di risultati processati.
%   finalVec = aggregateProcessedResults(processedFolder) legge tutti i file .mat in
%   processedFolder, estrae la cella resultCell (1×200 cell) da ciascun file,
%   calcola la media di ogni vettore all'interno di resultCell, quindi calcola la
%   media dei 200 valori medi ottenuti, restituendo un vettore finalVec di lunghezza
%   pari al numero di file.

    % Trova e ordina i .mat nella cartella di input
    files = dir(fullfile(processedFolder, '*.mat'));
    fileNames = sort({files.name});
    nFiles = numel(fileNames);

    % Preallocazione del vettore risultato
    finalVec = zeros(1, nFiles);

    % Ciclo su ciascun file
    for k = 1:nFiles
        % Carica resultCell dal file .mat
        data = load(fullfile(processedFolder, fileNames{k}));
        resultCell = data.resultCell;    % 1×200 cell array

        % Calcola la media di ciascun vettore in resultCell
        nCells = numel(resultCell);
        perMeans = zeros(1, nCells);
        for i = 1:nCells
            v = resultCell{i};
            perMeans(i) = mean(v);
        end

        % Media complessiva dei 200 valori medi
        finalVec(k) = mean(perMeans);
    end
end
