function community_variance(testFolder, dataFile)
% processTestFiles – Elabora i .mat di test estraendo la 2ª riga di outputCell
%   processTestFiles(testFolder, dataFile) carica la cella C da dataFile,
%   scansiona in ordine alfabetico tutti i .mat in testFolder, estrae
%   la seconda riga di outputCell (1×200 cell), esegue le tue operazioni
%   su ciascun vettore e salva i risultati in una nuova cartella
%   “processed_results” sotto lo stesso path.

    % 1. Carica C dal file dati
    sData = load(dataFile);
    C = sData.C;  % usala come ti serve

    % 2. Trova tutti i .mat in testFolder e ordina i nomi
    files = dir(fullfile(testFolder, '*.mat'));
    fileNames = sort({files.name});

    % 3. Crea cartella di output sotto il root comune
    [rootPath, ~, ~] = fileparts(testFolder);
    if isempty(rootPath)
        rootPath = pwd;
    end
    outputFolder = fullfile(rootPath, 'processed_results');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % 4. Cicla su ciascun file
    for k = 1:numel(fileNames)
        % 4.1 Carica il file di test
        testFile = fullfile(testFolder, fileNames{k});
        sTest = load(testFile);
        inputCell = sTest.outputCell;         % 3×200 cell
        rowData   = inputCell(2, :);          % 1×200 cell

        % 4.2 Prepara la cell per i risultati
        resultCell = cell(1, 200);

        % 4.3 Ciclo sui 200 vettori
        for i = 1:200

            opinion_end = double(rowData{i})/1000;
            comm = C{i};
            n_comm = max(comm);
            comm_var = zeros(n_comm,1);

            for j = 1:n_comm
                comm_agents = find(comm == j);
                comm_var(j) = var(opinion_end(comm_agents));  
            end 
        
            resultCell{i} = comm_var;
            i

        end

        % 4.4 Salva il risultato
        [~, nameNoExt, ~] = fileparts(fileNames{k});
        outFile = fullfile(outputFolder, [nameNoExt '_processed.mat']);
        save(outFile, 'resultCell');
    end
end
