mkdir("tables")
load data

[X_sim,Y_sim] = getRankScore("species","idMetabolites");
[X_sm,Y_sm] = getRankScore("species","metabolites");
[X_eim,Y_eim] = getRankScore("enzymes","idMetabolites");
[X_em,Y_em] = getRankScore("enzymes","metabolites");

%% species
t = NaN(numel(species.speciesID),4);

t(X_sim.ind,1) = X_sim.score;
t(X_sim.ind,2) = X_sim.rank;

t(X_sm.ind,3) = X_sm.score;
t(X_sm.ind,4) = X_sm.rank;

colNames = ["ID","Id-Mets-Score","Id-Mets-Rank", ...
    "All-Mets-Score","All-Mets-Rank"];
% all
writematrix([colNames; [species.speciesID t]], "tables/species.csv")
% top 10
minRank = min(t(:,[2 4]),[],2);
inds = find(minRank<11);
t = t(inds,:);
ID = replace(species.speciesID(inds),"_","-");
[~,order] = sort(minRank(inds));
t = t(order,:);
ID = ID(order);
writematrix([colNames;[ID round(t,2)]],'tables/species_top10.csv')

%% enzymes
t = NaN(numel(enzymes.enzymeID),4);

t(X_eim.ind,1) = X_eim.score;
t(X_eim.ind,2) = X_eim.rank;

t(X_em.ind,3) = X_em.score;
t(X_em.ind,4) = X_em.rank;

colNames = ["ID","Id-Mets-Score","Id-Mets-Rank", ...
    "All-Mets-Score","All-Mets-Rank"];
% all
writematrix([colNames; [enzymes.enzymeID t]], "tables/enzymes.csv", ...
    'QuoteStrings','none')
% top 15
minRank = min(t(:,[2 4]),[],2);
inds = find(minRank<16);
t = t(inds,:);
ID = replace(enzymes.enzymeID(inds),"_","-");
[~,order] = sort(minRank(inds));
t = t(order,:);
ID = ID(order);
writematrix([colNames;[ID round(t,2)]],'tables/enzymes_top15.csv', ...
    'QuoteStrings','none')

%% metabolites
t = NaN(numel(metabolites.best),8);

t(Y_sim.ind,1) = Y_sim.score;
t(Y_sim.ind,2) = Y_sim.rank;

t(Y_sm.ind,3) = Y_sm.score;
t(Y_sm.ind,4) = Y_sm.rank;

t(Y_eim.ind,5) = Y_eim.score;
t(Y_eim.ind,6) = Y_eim.rank;

t(Y_em.ind,7) = Y_em.score;
t(Y_em.ind,8) = Y_em.rank;

colNames = ["ID","Species-Id-Mets-Score","Species-Id-Mets-Rank", ...
    "Species-All-Mets-Score","Species-All-Mets-Rank", ...
    "Enzymes-Id-Mets-Score","Enzymes-Id-Mets-Rank", ...
    "Enzymes-All-Mets-Score","Enzymes-All-Mets-Rank"];
% all
writematrix([colNames; [metabolites.best t]], "tables/metabolites.csv", ...
    'QuoteStrings','none')
% top 15
minRank = min(t(:,2:2:8),[],2);
inds = find(minRank<16);
t = t(inds,:);
ID = replace(metabolites.best(inds),"_","-");
[~,order] = sort(minRank(inds));
t = t(order,:);
ID = ID(order);
writematrix([colNames;[ID round(t,2)]],'tables/metabolites_top15.csv')
% id only
t = t(:,[1 2 5 6]);
keep = ~isnan(t(:,1));
t = t(keep,:);
ID = ID(keep);
minRank = min(t(:,[2 4]),[],2);
[~,order] = sort(minRank);
t = t(order,:);
ID = ID(order);
writematrix([colNames([1 2 3 6 7]);[ID round(t,2)]], ...
    'tables/metabolites_top15_IdOnly.csv', ...
    'QuoteStrings','none')
