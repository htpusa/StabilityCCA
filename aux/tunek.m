function [k Xinds Yinds] = tunek(X,Y,A,B,SCCA)

% Binary search for optimal Top-k

maxIter = 100;
stopCrit = 1e-8;
rounds = 10;
K = 3;

px = numel(A.auc);
[~,rank] = sort([A.auc;B.auc],'descend');
range = [1 numel(A.auc)+numel(B.auc)];

iter = 0;
cont = 1;

while cont
    grid = linspace(range(1),range(2),5);
    grid = round([grid(2) grid(4)]);
    score = -ones(2,rounds);
    parfor r=1:rounds
        part = crossvalind("Kfold",size(X,1),K);
        test = part==1;
        Xtrain = X(~test,:); 
        Ytrain = Y(~test,:);
        Xtest = X(test,:);
        Ytest = Y(test,:); 
        for i=1:2
            [Xinds,Yinds] = stabilityScoreSelection(A,B,grid(i));
            if (numel(Xinds)>0) & (numel(Yinds)>0)
                if SCCA=="SCCAec"
                    [a,b] = SCCAec(Xtrain(:,Xinds),Ytrain(:,Yinds),...
                        'cxy',[sqrt(numel(Xinds)) sqrt(numel(Yinds))]);
                else
                    [a,b] = PMDCCA(Xtrain(:,Xinds),Ytrain(:,Yinds),...
                        'cxy',[sqrt(numel(Xinds)) sqrt(numel(Yinds))]);
                end
                score(i,r) = corr(Xtest(:,Xinds)*a(:,1),Ytest(:,Yinds)*b(:,1));
            end
        end
    end
    score = mean(score,2);
    [M,I] = max(score);
    range = sort([range(I),mean(range)]);
    m = min(score);
    imp = M-m;
    iter = iter+1;
    cont = (imp>stopCrit) & (iter<maxIter) & (range(1)<range(2));
end

k = round(mean(range));
if k<max(find(rank>px,1),find(rank<=px,1))
    k = max(find(rank>px,1),find(rank<=px,1));
end
[Xinds,Yinds] = stabilityScoreSelection(A,B,k);