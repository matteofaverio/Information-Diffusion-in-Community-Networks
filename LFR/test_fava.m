%% PARAMETRI E GENERAZIONE RETE

clear 

n = 3000;
gamma = 3;
gamma_c = 3;
d = 10;
d_min = 7; 
mu = 0.1;

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

M = (c == c');
sameCommCounts1 = sum(AA .* M, 2);
degrees1 = sum(AA, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
mean(fractions1)

%%

n = 1000;
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7; 

t = 50;
tests = 5;

mu = linspace(0.01,0.99,t);
AVG_MU = zeros(t,tests);


for i = 1:t

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

    end

    fprintf('mu = %d\n\n', mu(i));

end



%%

nonzero_count = sum(AVG_MU ~= 0, 2); 
row_sum = sum(AVG_MU, 2);  
media_righe = row_sum ./ nonzero_count;

x = mu;
y = media_righe - mu';

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

%%

clc; clear; close all;

n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.85;
%%
tic
[A,~,c,dd] = LFR2(n,d,mu,gamma, gamma_c, d_min);
toc

%%

www = [816
   939
   963
   967];

%% 

t = [ 2 3 ; 4 3 ; 3 2]

t(t(:,1) == 3,:) = [];

t


