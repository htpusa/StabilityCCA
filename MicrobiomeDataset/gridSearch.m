function [cx cy] = gridSearch(X,Y)

gridSize = 15;
cx = logspace(log10(1),log10(sqrt(size(X,2))),gridSize);
cy = logspace(log10(1),log10(sqrt(size(Y,2))),gridSize);

rounds = 10;
K = 5;

score = zeros(gridSize,gridSize,rounds,K);

for r=1:rounds
   part = crossvalind("Kfold",size(X,1),K);
   for k=1:K
       test = part==k;
       Xtrain = X(~test,:);Ytrain = Y(~test,:);
       Xtest = X(test,:);Ytest = Y(test,:);
       for i=1:gridSize
           parfor j=1:gridSize
               [a b] = SCCAec(Xtrain,Ytrain,'cxy',[cx(i) cy(j)]);
               score(i,j,r,k) = corr(Xtest*a,Ytest*b);
           end
       end
   end
end

score = mean(score,3:4);
[~,I] = max(score,[],'all','linear');
[i,j] = ind2sub(size(score),I);
cx = cx(i);
cy = cy(j);