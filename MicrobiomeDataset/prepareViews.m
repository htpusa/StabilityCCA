function [X,Y] = prepareViews(Xvar,Yvar)

load('data.mat')

%% X view
if Xvar=="species"
    zeroVar = var(species.data)==0;
    X.ind = find(~zeroVar);
    X.data = species.data(:,~zeroVar);
    X.names = species.speciesID(~zeroVar);
    X.data = normalize(X.data);
    X.status = species.metaData.status;
elseif Xvar=="genus"
    zeroVar = var(genus.data)==0;
    X.ind = find(~zeroVar);
    X.data = genus.data(:,~zeroVar);
    X.names = genus.genusID(~zeroVar);
    X.data = normalize(X.data);
    X.status = genus.metaData.status;
elseif Xvar=="enzymes"
    zeroVar = var(enzymes.data)==0;
    X.ind = find(~zeroVar);
    X.data = enzymes.data(:,~zeroVar);
    X.names = enzymes.enzymeID(~zeroVar);
    X.data = normalize(X.data);
    X.status = enzymes.metaData.status;
    X.shortnames = X.names;
    for i=1:numel(X.names)
        tmp = split(X.names(i),': ');
        X.shortnames(i) = tmp(1);
    end
end

%% Y view
if Yvar=="idMetabolites"
    exact = ~ismissing(metabolites.exactMatch);
    zeroVar = var(metabolites.data)==0;
    Y.ind = find(exact & (~zeroVar)');
    Y.data = metabolites.data(:,exact & (~zeroVar)');
    Y.names = metabolites.best(exact & (~zeroVar)');
    Y.data = normalize(Y.data);
elseif Yvar=="metabolites"
    zeroVar = var(metabolites.data)==0;
    Y.ind = find(~zeroVar);
    Y.data = metabolites.data(:,~zeroVar);
    Y.names = metabolites.best(~zeroVar);
    Y.data = normalize(Y.data);
end