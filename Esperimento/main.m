%% ISTOGRAMMA3D CON OPINIONI FINALI AL VARIARE DI EPSILON

histOpsVSEps('test_mu085_processed');

%% PLOT 25 ISTOGRAMMI (PIANI DI TAGLIO DI OPSvsEPS)

N = 1e4;
test = 'test_mu080_processed';
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

iterazioni('Data');


