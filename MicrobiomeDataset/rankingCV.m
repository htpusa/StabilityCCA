function [stabCCAall,scca,sccaRank] = rankingCV(Xstruct,Ystruct)

X = Xstruct.data;
Y = Ystruct.data;

k = 5;
rounds = 50;
maxRank = 50;

px = size(X,2);
py = size(Y,2);

stabCCA.score = zeros(rounds, maxRank);
stabCCA.ranking = zeros(px+py,rounds);

% odds ratio
stabCCA2.score = zeros(rounds, maxRank);
stabCCA2.ranking = zeros(px+py,rounds);

% alternate, auc
stabCCA3.score = zeros(rounds, maxRank);
stabCCA3.ranking = zeros(px+py,rounds);

% alternate, odds ratio
stabCCA4.score = zeros(rounds, maxRank);
stabCCA4.ranking = zeros(px+py,rounds);

% first to thresh
th = 0.8;
stabCCA5.score = zeros(rounds, maxRank);
stabCCA5.ranking = zeros(px+py,rounds);

scca.score = zeros(rounds, 1);
scca.A = zeros(px,rounds);
scca.B = zeros(py,rounds);

sccaRank.score = zeros(rounds, maxRank);
sccaRank.ranking = zeros(px+py,rounds);

cx = 0;
cy = 0;

for r=1:rounds
    fprintf('\t Round %d of %d\n',r,rounds);
    test = crossvalind('Kfold',size(X,1),k) == 1;
    Xtrain = X(~test,:);
    Xtest = X(test,:);
    Ytrain = Y(~test,:);
    Ytest = Y(test,:);
    % SCCA
    % find best parameters once
    if r==1
        [cx,cy] = gridSearch(Xtrain,Ytrain);
    end
    [Atmp,Btmp] = SCCAec(Xtrain,Ytrain,'cxy',[cx,cy]);
    scca.score(r) = corr(Xtest*Atmp,Ytest*Btmp);
    scca.A(:,r) = Atmp; scca.B(:,r) = Btmp;

    [~, sccarank] = sort([abs(Atmp)/max(abs(Atmp)); abs(Btmp)/max(abs(Btmp))],'descend');
    sccaRank.ranking(:,r) = sccarank;
    % stability CCA
    [A,B] = stabilityCCA(Xtrain,Ytrain,'verbose',0);
    [~, aucrank] = sort([A.auc; B.auc],'descend');
    stabCCA.ranking(:,r) = aucrank;
    
    Ar = mean(A.probs ./ mean(A.probs,2),1)';
    Br = mean(B.probs ./ mean(B.probs,2),1)';
    [~, aucrank2] = sort([Ar; Br],'descend');
    stabCCA2.ranking(:,r) = aucrank2;

    [~, aucrankA] = sort(A.auc,'descend');
    [~, aucrankB] = sort(B.auc,'descend'); aucrankB=aucrankB+px;
    aucrankA1 = aucrankA(1:min(px,py));
    aucrankB1 = aucrankB(1:min(px,py));
    aucrank3 = reshape([aucrankA1';aucrankB1'],2*min(px,py),1);
    aucrank3 = [aucrank3; aucrankA(min(px,py)+1:end); aucrankB(min(px,py)+1:end)];
    stabCCA3.ranking(:,r) = aucrank3;

    [~, aucrankA] = sort(Ar,'descend');
    [~, aucrankB] = sort(Br,'descend'); aucrankB=aucrankB+px;
    aucrankA1 = aucrankA(1:min(px,py));
    aucrankB1 = aucrankB(1:min(px,py));
    aucrank4 = reshape([aucrankA1';aucrankB1'],2*min(px,py),1);
    aucrank4 = [aucrank4; aucrankA(min(px,py)+1:end); aucrankB(min(px,py)+1:end)];
    stabCCA4.ranking(:,r) = aucrank4;

    % first to thresh
    Afirst = 101*ones(px,1);
    for pxi=1:px
        ind = find(A.probs(:,pxi)>th,1);
        if ~isempty(ind)
            Afirst(pxi) = ind;
        end
    end
    Bfirst = 101*ones(py,1);
    for pyi=1:py
        ind = find(B.probs(:,pyi)>th,1);
        if ~isempty(ind)
            Bfirst(pyi) = ind;
        end
    end
    [~,aucrank5] = sort([Afirst;Bfirst],'ascend');
    stabCCA5.ranking(:,r) = aucrank5;

    for nSel=1:maxRank
        % scca ranking
        sel = sccarank(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            sccaRank.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            sccaRank.score(r,nSel) = NaN;
        end
        % stability CCA ranking
        sel = aucrank(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            stabCCA.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            stabCCA.score(r,nSel) = NaN;
        end
        % stability CCA ranking based on ratio to average
        sel = aucrank2(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            stabCCA2.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            stabCCA2.score(r,nSel) = NaN;
        end
        % stability CCA alternating ranking based on auc
        sel = aucrank3(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            stabCCA3.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            stabCCA3.score(r,nSel) = NaN;
        end
        % stability CCA alternating ranking based on odds ratio
        sel = aucrank4(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            stabCCA4.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            stabCCA4.score(r,nSel) = NaN;
        end
        % stability CCA ranking based on first to threshold
        sel = aucrank5(1:nSel);
        Ainds = sel(sel<=px);
        Binds = sel(sel>px)-px;
        if (~isempty(Ainds)) & (~isempty(Binds))
            [Atmp,Btmp] = canoncorr(Xtrain(:,Ainds), Ytrain(:,Binds));
            stabCCA5.score(r,nSel) = corr(Xtest(:,Ainds)*Atmp(:,1), ...
                                  Ytest(:,Binds)*Btmp(:,1));
        else
            stabCCA5.score(r,nSel) = NaN;
        end
    end
end

stabCCAall.stabCCA = stabCCA;
stabCCAall.stabCCA2 = stabCCA2;
stabCCAall.stabCCA3 = stabCCA3;
stabCCAall.stabCCA4 = stabCCA4;
stabCCAall.stabCCA5 = stabCCA5;

