function [Xinds,Yinds] = stabilityScoreSelection(A,B,parameter)

% Select variables based on stability scores, ie the area under stability
% path (A.auc and B.auc)
% A and B should be outputs of stabilityCCA
% If paramater is an integer, Top-k variables are selected
% If parameter is in (0,1), variables with auc above parameter are selected
% (threshold selection)

px = numel(A.auc);
auc = [A.auc;B.auc];

if parameter<1
    inds = find(auc>=parameter);
else
    [~,rank] = sort(auc,'descend');
    inds = rank(1:parameter);
end
Xinds = inds(inds<=px);
Yinds = inds(inds>px)-px;
