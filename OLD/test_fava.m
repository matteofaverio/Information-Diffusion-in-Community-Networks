%% PARAMETRI E GENERAZIONE RETE

clear 

n = 3000;
gamma = 3;
gamma_c = 3;
d = 10;
d_min = 7; 
mu = 0.5;

tic
[A,~,c,dd] = LFR2(n,d,mu,gamma, gamma_c, d_min); 
toc

M = (c == c');
sameCommCounts1 = sum(A .* M, 2);
degrees1 = sum(A, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
mean(fractions1)
%%
tic
AA = rewiring(A,c,mu,10000);
toc
%%

M = (c == c');
sameCommCounts1 = sum(AA .* M, 2);
degrees1 = sum(AA, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
mean(fractions1)


%%

clear
clc

% Parametri
n       = 10000;
gamma   = 3;
gamma_c = 2;
d       = 10;
d_min   = 7;


tests   = 200;
mu      =  [0.51 0.56 0.62 0.64 0.70 0.76 0.78 0.84 0.88 0.92];

% Cartella di output
outDir = 'output';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for i = 1:length(mu)
    % Nome del file per questo valore di mu
    outFile = fullfile(outDir, sprintf('data_mu%02d.mat', i));
    if exist(outFile, 'file')
        delete(outFile);
    end

    % 1) Preallocazione su disco
    M       = cell(1, tests);   % per le matrici A
    C       = cell(1, tests);   % per i vettori c
    AVG_mu = zeros(1, tests);  % per i risultati AVG_MU(i,:)
    save(outFile, 'M', 'C', 'AVG_mu', '-v7.3');
    clear M C AVG_mu;

    % 2) Apre il MAT in scrittura parziale
    matObj = matfile(outFile, 'Writable', true);

    % 3) Loop interno: genera A e c, calcola, salva e libera RAM
    for j = 1:tests
        [A, ~, c, ~] = LFR2(n, d, mu(i), gamma, gamma_c, d_min);
      
        % calcolo AVG_MU(i,j)
        Mcomm           = (c == c');
        sameCommCounts  = sum(A .* Mcomm, 2);
        degs            = sum(A, 2);
        fracs           = sameCommCounts ./ degs;
        avg_val         = mean(fracs);

        % riduco la matrice a una matrice sparsa
        A = sparse(A);

        % scrive parzialmente su disco
        matObj.M(1, j)       = {A};
        matObj.C(1, j)       = {c};
        matObj.AVG_mu(1, j) = avg_val;

        clear A c
        fprintf('  [%2d,%3d] salvato in %s\n', i, j, outFile);
    end

    fprintf('Completato file %s (mu=%.2f)\n\n', outFile, mu(i));
end

%% 
clear
clc

n = 10000;
gamma = 3;
gamma_c = 2;
d = 10;
d_min = 7; 

t = 50;
tests = 20;

mu = linspace(0.47,0.97,t+1);
AVG_MU = zeros(t,tests);

for i = 1:t+1

    for j = 1:tests

        %tic
        [A,~,c,dd] = LFR2(n,d,mu(i),gamma, gamma_c, d_min); 
        %toc
        M = (c == c');
        sameCommCounts1 = sum(A .* M, 2);
        degrees1 = sum(A, 2);
        fractions1 = sameCommCounts1 ./ degrees1;
        fractions1(degrees1 == 0) = 0;
        AVG_MU(i,j) = mean(fractions1);

        j

    end

    fprintf('mu = %d\n\n', mu(i));

end
%%
nonzero_count = sum(AVG_MU ~= 0, 2); 
row_sum = sum(AVG_MU, 2);  
media_righe = row_sum ./ nonzero_count;

x = mu';
y = media_righe - x;

coeff = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);  % crea una griglia fine
y_fit = polyval(coeff, x_fit);
plot(x, y)%, x_fit, y_fit, '-')


%%

err = mu'-media_righe;
err_mediato = zeros(t,1);
h = 6;

for i = 1:t
        idx_start = max(1, i - h);
        idx_end   = min(t, i + h);
        err_mediato(i) = mean(err(idx_start:idx_end));
end
plot(mu,err_mediato);


%% ANALISI GRADO NODI

subplot(2,1,1)
histogram(dd,'BinLimits',[0,50])
subplot(2,1,2)
histogram(sum(A,2),'BinLimits',[0,50])

%% INDIVIDUAZIONE COMUNITÀ

Q_LFR = community_louvain(A);

NMI_LFR = nmi(c,Q_LFR);

fprintf('Numero di comunità rilevate: %d\n', max(Q_LFR));
fprintf('Normalized Mutual Information: %4f\n',NMI_LFR)
