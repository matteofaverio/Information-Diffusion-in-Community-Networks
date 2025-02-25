%% Generazione di una matrice di adiacenza test
% number of nodes
clc, clear
n = 15;
% Adjacency matrix
A = rand(n)*2-1;
A(A <= 0.7) = 0;
A = A./A;
A(isnan(A)) = 0;
A = A.*rand(n);
A = round(A,3);
G = digraph(A,'omitselfloops');
plot(G,'layout','circle','EdgeLabel',G.Edges.Weight)

%% Test Independent Cascade Model
clc,clear
v0 = [1, 0, 0, 0, 0, 0]';
A = [ 0, 0, 0.7, 0.8, 0, 0;
    0.3, 0, 0, 0.5, 0, 0;
    0, 0, 0, 0, 0.6, 0.1;
    0, 0, 0.2, 0, 0, 0.3;
    0.4, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0.9, 0];
[v,t] = opinionProp_IC(A,v0);
v

%% Test Independent Cascade con matrici più grandi

clc, clear, close all
n = 50;

% Matrice di incidenza pesata
A = rand(n)*2-1;
A(A <= 0.8) = 0;
A = A./A;
A(isnan(A)) = 0;
A = A.*(0.8*rand(n));
A = round(A,3);

%% Matrice small World
clc,clear
n = 100;
k = 4;
p = 0.3;
A = network_SM(n,k,p);
W = rand(n).*A;

for i=1:n
    W(n,:) = W(n,:)/sum(W(n,:));
end

%% insieme dei nodi di partenza
v = 1:n;
seed1 = randi([1,n]);
seed2 = randi([1,n]);
v0 = zeros(n,1);
v0(seed1) = 1;
if seed1 ~= seed2
    v0(seed2) = 1;
end
tr = 0.1*rand(n,1);

% G = digraph(A)

[opLT,tLT] = opinionProp_LT(W,v0,tr);
[opIC,tIC] = opinionProp_IC(W,v0);

%%
G = digraph(A,'omitselfloops');

for i = t-1:-1:1
    figure(i)
    p = plot(G,ArrowPosition=0.95,EdgeAlpha=0.3);
    selectedNodes = v(op(:,i) == 1);
    highlight(p, selectedNodes, 'NodeColor','r');
    p.NodeLabel = [];
    
end

%% test Barabàsi-Albert Network
clc, clear, close all
n = 10000;
m = 20;

A = network_BA(n,m);
G = graph(A);
figure(1)
% subplot(1,2,1)
p = plot(G,EdgeAlpha=0.3);
for i = 1:n
    d = sum(A(i,:));
    dim = 0.2;
    highlight(p,i,'MarkerSize',dim*G.degree(i),'NodeColor','b')
    p.NodeLabel = [];
end

D = degree(G);
mean(D)
%%
W = rand(n).*A;
for i=1:n
    W(n,:) = W(n,:)/sum(W(n,:));
end

v = 1:n;
v0 = binornd(1,0.01,n,1);


% G = digraph(A)

[op,t] = opinionProp_IC(W,v0);

%% test Barabàsi-Albert LCD
 clear
n =1000;
m = 20;

A = network_BA_LCD(n,m);
G = graph(A);
figure(1)
% subplot(1,2,2)
p = plot(G,EdgeAlpha=0.3);
for i = 1:n
    dim = 0.2;
    highlight(p,i,'MarkerSize',dim*G.degree(i),'NodeColor','b')
    p.NodeLabel = [];
end

D = degree(G);
mean(D)