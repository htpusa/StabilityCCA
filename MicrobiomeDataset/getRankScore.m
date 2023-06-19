function [X,Y] = getRankScore(Xvar,Yvar)

load("res/"+Xvar+"_"+Yvar+"/stableCCA")
[X,Y] = prepareViews(Xvar,Yvar);
[~,ranking] = sort([A.auc; B.auc], 'descend');

X.score = A.auc;
[~,X.rank] = ismember((1:numel(X.score))',ranking);

Y.score = B.auc;
[~,Y.rank] = ismember(((1:numel(Y.score))+numel(X.score))',ranking);