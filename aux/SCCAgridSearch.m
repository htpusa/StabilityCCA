function cxy = SCCAgridSearch(X,Y,SCCA)

% Find optimal sparsity parameters using grid search
% SCCA is either "SCCAec" or "PMDCCA"

rounds = 10;
K = 3;
gridSize = 15;
xgrid = linspace(1,sqrt(size(X,2)),gridSize);
ygrid = linspace(1,sqrt(size(Y,2)),gridSize);
score = zeros(gridSize,rounds);
for ri=1:rounds
    part = crossvalind("Kfold",size(X,1),K);
    test = part==1;
    Xtrain = X(~test,:);
    Ytrain = Y(~test,:);
    Xtest = X(test,:);
    Ytest = Y(test,:);
    parfor i=1:gridSize
        for j=1:gridSize
            if SCCA=="SCCAec"
                [a,b] = SCCAec(Xtrain,Ytrain,'cxy',[xgrid(i),ygrid(j)]);
            else
                [a,b] = PMDCCA(Xtrain,Ytrain,'cxy',[xgrid(i),ygrid(j)]);
            end
            u = Xtest*a;
            v = Ytest*b;
            score(i,j,rounds) = corr(u,v);
        end
    end
end
score = mean(score,3);
[~,I] = max(score,[],"all");
[i,j] = ind2sub([gridSize gridSize],I);
cxy = [xgrid(i),ygrid(j)];