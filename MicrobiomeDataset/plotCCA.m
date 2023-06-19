function plotCCA(Xvar,Yvar,XvarName,YvarName,nSel)

[X,Y] = prepareViews(Xvar,Yvar);

path = "res/" + Xvar + '_' + Yvar;
load(path + "/stableCCA.mat");

% thresh = 0.8;
% cont = 1;
% nSel = 2;
% [~,rank] = sort([A.auc; B.auc], 'descend');
% while cont
%     sel = rank(1:nSel);
%     Ainds = sel(sel<=numel(A.auc));
%     Binds = sel(sel>numel(A.auc))-numel(A.auc);
%     if numel(Ainds)>0 & (numel(Binds)>0)
%         [a,b,r,u,v] = canoncorr(X.data(:,Ainds),Y.data(:,Binds));
%         cont = r(1)<thresh;
%     end
%     nSel = nSel + 1;
% end

[~,rank] = sort([A.auc; B.auc], 'descend');
cancor = nan(50,1);
for i=2:50
    sel = rank(1:i);
    Ainds = sel(sel<=numel(A.auc));
    Binds = sel(sel>numel(A.auc))-numel(A.auc);
    if numel(Ainds)>0 & (numel(Binds)>0)
        [~,~,r] = canoncorr(X.data(:,Ainds),Y.data(:,Binds));
        cancor(i) = r(1);
    end
end
figure
    plot(1:50,cancor,'o-k','LineWidth',2)
    xlim([0 51]);ylim([0 1])
    title(XvarName + ' - ' + YvarName)
    xlabel('number of variables');ylabel('canonical correlation')
    exportgraphics(gcf,'plots/'+Xvar+'_'+Yvar+'_CCAcor.png', ...
        'Resolution',150)

sel = rank(1:nSel);
Ainds = sel(sel<=numel(A.auc));
Binds = sel(sel>numel(A.auc))-numel(A.auc);
[a,b,r,u,v] = canoncorr(X.data(:,Ainds),Y.data(:,Binds));

su = sign(mean(u(X.status~='Control'),1) - ...
    mean(u(X.status=='Control'),1));
sv = sign(mean(v(X.status~='Control'),1) - ...
    mean(v(X.status=='Control'),1));
a = su*a; b = sv*b; u = su*u; v = sv*v;

% Variables
figure
    gscatter(u(:,1),v(:,1),X.status,'rkb','os^',10)
    legend('Location','northwest')
    xlabel(XvarName + ' variable'); ylabel(YvarName + ' variable')
    % title(sprintf('Correlation=%.2f',r(1)))
    title(XvarName + ' - ' + YvarName)
set(gcf,'Position',[100,100,600,600])
fontsize(gcf,14,"points")
exportgraphics(gcf,'plots/'+Xvar+'_'+Yvar+'_CCAvars_top'+string(nSel)+'.png', ...
    'Resolution',300)

% Coefficients and pairwise contributions
Xnames = replace(X.names,'_','-');
Ynames = replace(Y.names,'_','-');

if Xvar=='enzymes'
    Xnames = X.shortnames;
end

pc = pairwiseContributions(X.data(:,Ainds),Y.data(:,Binds),a(:,1),b(:,1));
Xdist = pdist(pc);
Xtree = linkage(Xdist);
Xord = optimalleaforder(Xtree,Xdist);
Ydist = pdist(pc');
Ytree = linkage(Ydist);
Yord = optimalleaforder(Ytree,Ydist);
Xtmp = Xnames(Ainds);
Ytmp = Ynames(Binds);

tiledlayout(2,2,'TileSpacing','tight')
    nexttile
    limits = [0,max(pc,[],'all')];
    heatmap(Xtmp(Xord),Ytmp(Yord),round(pc(Xord,Yord),2)',...
        'ColorbarVisible',0,'ColorLimits',limits);
    ax = gca;
    ax = struct(ax);
    ax.Axes.XAxisLocation = 'top';
    ax.XAxis.TickLabelRotation=90;

    nexttile
    barh(b(fliplr(Yord),1),'r');
    ylim([0.6 numel(b(:,1))+0.4]);yticklabels([])
    %xticks([-0.4,0,0.4])
    %yticks(1:numel(Binds));yticklabels(Y.names(Binds))
    
    nexttile
    bar(a(Xord,1),'r');
    xlim([0.6 numel(a(:,1))+0.4])
    xticklabels([]); %yticks([-0.6,-0.3,0])
    %xticks(1:numel(Ainds));xticklabels(Xnames(Ainds));xtickangle(90)
    %ylabel('coefficient')

set(gcf,'Position',[100,100,900,900])
fontsize(gcf,14,"points")
exportgraphics(gcf,'plots/'+Xvar+'_'+Yvar+'_CCAcoeff_top'+string(nSel)+'.png', ...
    'Resolution',300)

% combine
% tiledlayout(4,4,'TileSpacing','tight');    
%     nexttile([2 4])
%     image(imread('plots/' + Xvar + '_' + Yvar + '_top10vars.png'));
%     xticks([]);yticks([])
%     nexttile(10,[2 2])
%     image(imread('plots/' + Xvar + '_' + Yvar + '_top10coeff.png'));
%     xticks([]);yticks([])
% exportgraphics(gcf, ...
%     'plots/'+Xvar+'_'+Yvar+'_top10.png','Resolution',600)

% ROC
figure
    [tpr,fpr,auc1] = rocCurv(u(:,1)+v(:,1),X.status~="Control");
    plot(fpr,tpr,'k','LineWidth',2)
    xlabel('FPR');ylabel('TPR')
    hold on
    cd = X.status~="UC";
    [tpr,fpr,auc2] = rocCurv(u(cd,1)+v(cd,1),X.status(cd)~="Control");
    plot(fpr,tpr,':r','LineWidth',2)
    legend({sprintf('IBD vs. Control: %.2f',auc1), ...
        sprintf('CD vs. Control: %.2f',auc2)},...
        'Location','southeast')
    title(XvarName + ' - ' + YvarName)

    exportgraphics(gcf,'plots/'+Xvar+'_'+Yvar+'_CCAroc_top'+string(nSel)+'.png', ...
    'Resolution',300)



