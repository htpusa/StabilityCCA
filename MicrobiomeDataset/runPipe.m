function runPipe(Xvar,Yvar)

path = "res/" + Xvar + "_" + Yvar;
mkdir(path)

[X,Y] = prepareViews(Xvar,Yvar);

%% Run StabilityCCA
if ~isfile(path +'/stableCCA.mat')
    [A,B] = stableCCA(X.data,Y.data);
    save(path +'/stableCCA.mat','A','B');
end

%% Run cross-validation
if ~isfile(path +'/ranking.mat')
    [stabCCA,scca,sccaRank] = rankingCV(X,Y);
    save(path +'/ranking.mat','stabCCA','scca','sccaRank');
end