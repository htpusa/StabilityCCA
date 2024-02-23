function [thr Xinds Yinds] = tuneThr(X,Y,A,B,SCCA)

% Find an optimal threshold value via cross-validation

gridSize = 100;
rounds = 10;
K = 3;

auc = [A.auc;B.auc];
grid = linspace(min(auc),max(auc),gridSize);
score = zeros(gridSize,rounds);

for r=1:rounds
    part = crossvalind("Kfold",size(X,1),K);
    test = part==1;
    Xtrain = X(~test,:);
    Ytrain = Y(~test,:);
    Xtest = X(test,:);
    Ytest = Y(test,:);
    parfor i=1:gridSize
        [Xinds,Yinds] = stabilityScoreSelection(A,B,grid(i));
        if numel(Xinds)>0 & (numel(Yinds)>0)
            if SCCA=="SCCAec"
                [a,b] = SCCAec(Xtrain(:,Xinds),Ytrain(:,Yinds),...
                    'cxy',[sqrt(numel(Xinds)) sqrt(numel(Yinds))]);
            else
                [a,b] = PMDCCA(Xtrain(:,Xinds),Ytrain(:,Yinds),...
                    'cxy',[sqrt(numel(Xinds)) sqrt(numel(Yinds))]);
            end
            score(i,r) = corr(Xtest(:,Xinds)*a,Ytest(:,Yinds)*b);
        end
    end
end
score = mean(score,2);
[~,I] = max(score);
thr = grid(I);
[Xinds,Yinds] = stabilityScoreSelection(A,B,thr);