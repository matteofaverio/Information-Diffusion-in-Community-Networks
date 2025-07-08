 %% ISTOGRAMMA3D CON OPINIONI FINALI AL VARIARE DI EPSILON

histOpsVSEps('test_mu050_processed');
histOpsVSEps('test_mu055_processed');
histOpsVSEps('test_mu060_processed');
histOpsVSEps('test_mu065_processed');
histOpsVSEps('test_mu070_processed');
histOpsVSEps('test_mu075_processed');
histOpsVSEps('test_mu080_processed');
histOpsVSEps('test_mu085_processed');
histOpsVSEps('test_mu090_processed');
%%
histOpsVSEps('test_mu095_processed');



%% PLOT 25 ISTOGRAMMI (PIANI DI TAGLIO DI OPSvsEPS)

N = 1e4;
test = 'test_mu095_processed';
files = dir(fullfile(test,'*.mat'));

names = sort({files.name});

for i = 1:numel(files)

    S = load(fullfile(test, names{i}));
    fn = fieldnames(S);
    C  = S.(fn{1});

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

%% PLOT ITERAZIONI DEI 10 TEST

iterazioni('Results');

%% 10 ISTOGRAMMI PER LO STESSO EPS

epsilon = 0.09; % epsilon da analizzare
hist_per_mu('Results',epsilon*100);

%% IDENTIFICAZIONE CLUSTERS

epsilon = linspace(0.01,0.25,25);
n_clusters = zeros(10,25); 
var_clusters = zeros(10,25);
outliers = zeros(10,25);
meanH = zeros(10,25);
meanL = zeros(10,25);
mus = 50:5:95;

for k = 1:10 %:numel(mus)
    tic
    % costruiamo il nome della cartella, con padding a 3 cifre
    test = sprintf('test_mu%03d_processed', mus(k));
    
    % controlliamo che esista davvero
    if ~isfolder(test)
        warning('Cartella non trovata: %s', test);
        continue;
    end

    files = dir(fullfile(test,'*.mat'));
    names = sort({files.name});
    
    for i = 1:25
        
        S = load(fullfile(test, names{i}));
        fn = fieldnames(S);
        C  = S.(fn{1});
        epsilon_i = epsilon(i)/4; % distanza massima tra punti nello stesso cluster
        minpts = 700;    % minimo numero di punti in un cluster
        isolated = zeros(1,200);
        meanHigh = zeros(1,200); meanLow = zeros(1,200);
        n = zeros(1,200);
    
        for j = 1:200
            x = C{2,j};
            data = double(x)/1000;
            idx = dbscan(data, epsilon_i, minpts);
            if max(idx)>0 
                v = data((idx==-1));
                meanHigh(j) = mean( v(v>0.5) );   % media dei valori > 0.5
                meanLow(j)  = mean( v(v<0.5) );   % media dei valori < 0.5
                isolated(j) = length(v);
            end
            n(j) = max(0,max(idx));
        end
    
        % media e varianza numero di cluster
        n_clusters(k,i) = mean(n);
        var_clusters(k,i) = var(n);
        % tolgo casi senza cluster
        isolated = isolated(isolated~=0);
        meanHigh = meanHigh(meanHigh~=0);
        meanLow = meanLow(meanLow~=0);
        % media numero di outlier e i valori attorno a cui stanno
        outliers(k,i) = mean(isolated);
        meanH(k,i) = mean(meanHigh);
        meanL(k,i) = mean(meanLow);
        fprintf('\n\nepsilon = %d\n',epsilon(i));
        
    end
    
    fprintf('Processata cartella %s (%d di %d)\n', test, k, numel(mus));
    toc
end

save('clusters4.mat','n_clusters','var_clusters','outliers','meanL',"meanH");


%% 

x = outputCell{2,4};
data = double(x)/1000;
minpts = 1000; 
epsilon = 0.1;
 
idx = dbscan(data, epsilon, minpts);
n = max(idx);
clusters = zeros(n,2);

for i = 1:n
    idx_i = find(idx==i);
    clusters(i,:) = [length(idx_i) , mean(data(idx_i))];
end
format shortG
clusters
histogram(data)


%%

community_variance('Results/test_mu090_processed','Data/data_mu01.mat');


%%

testFolder = 'Results/test_mu090_processed';
files = dir(fullfile(testFolder, '*.mat'));
fileNames = sort({files.name});

variances_050_tot = zeros(1,25);

for k = 1:numel(fileNames)
    
    testFile = fullfile(testFolder, fileNames{k});
    sTest = load(testFile);
    inputCell = sTest.outputCell;         % 3×200 cell
    rowData   = inputCell(2, :);          % 1×200 cell

    % 4.2 Prepara la cell per i risultati
    var_per_test = zeros(1, 200);

    % 4.3 Ciclo sui 200 vettori
    for i = 1:200
        i

        opinion_end = double(rowData{i})/1000;
        var_per_test(i) = var(opinion_end);

    end

    variances_050_tot(k) = mean(var_per_test);

end

plot(1:25,variances_050_tot,1:25,variances_050_comm);

%%

opinion_cell = outputCell(2,:);
idx = 5;
op_end = double(opinion_cell{idx})/1000;
histogram(op_end)


c = C{idx};
for i = 1:max(c)
    comm = find(c==i);
    mean_comm(i) = mean(op_end(comm));
    var_comm(i) = var(op_end(comm));
end

figure;
histogram(mean_comm);

