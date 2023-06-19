function [tpr,fpr,auc] = rocCurv(values,Y)

thr = [min(values)-1, linspace(min(values),max(values),100), max(values)+1];
tpr = zeros(numel(thr),1);
fpr = tpr;
for t=1:numel(thr)
    pred = values>thr(t);
    tpr(t) = sum(pred & Y)/sum(Y);
    fpr(t) = sum(pred & ~Y)/sum(~Y);
end

sorted = sortrows([fpr tpr],'ascend');
auc = sum(...
    diff([0;sorted(:,1)]) .* sorted(:,2) - ...
    (diff([0;sorted(:,1)]) .* diff([0;sorted(:,2)]) /2) ...
    );