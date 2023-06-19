function plotStabPath(Xvar,Yvar,XvarName,YvarName,nSel)

[X,Y] = prepareViews(Xvar,Yvar);

path = "res/" + Xvar + '_' + Yvar;
load(path + "/stableCCA.mat");

[~,rank] = sort([A.auc; B.auc], 'descend');
sel = rank(1:nSel);
Ainds = sel(sel<=numel(A.auc));
Binds = sel(sel>numel(A.auc))-numel(A.auc);

% quartiles
Xq = quantile(A.probs',3)';
Yq = quantile(B.probs',3)';
% quantiles
Xqn = quantile(A.probs',9)';
Yqn = quantile(B.probs',9)';

figure
    subplot(2,1,1);plot(1:size(A.c,1),A.probs(:,Ainds),'LineWidth',2)
    xticks(floor(linspace(1,size(A.c,1),5)))
    %xticklabels(round(A.c(xticks,1),2))
    xticklabels("("+round(A.c(xticks,1),1)+","+round(B.c(xticks,1),1)+")")
    ylim([0 1.01]);xlabel('c_x');ylabel('probability');
    xlabel('(c_x,c_y)')
    title(XvarName); hold on;
    %plot(1:size(A.c,1),Xq(:,1),'--k','LineWidth',2)
    %plot(1:size(A.c,1),Xq(:,2),'k','LineWidth',2)
    %plot(1:size(A.c,1),Xq(:,3),'-.k','LineWidth',2)
    plot(1:size(A.c,1),Xqn(:,[1:4,6:9]),':k','LineWidth',1.5)
    plot(1:size(A.c,1),Xqn(:,5),'-.k','LineWidt',1.75)
    Anames = replace(X.names(Ainds),'_',' ');
    %legend([Anames;"Q_1";"Q_2";"Q_3";"90% quantile"],'Location','northeastoutside', ...
        %'AutoUpdate','off')
    legend([Anames;"Deciles"],'Location','northeastoutside', ...
        'AutoUpdate','off')
    tmp = A.probs;
    tmp(:,Ainds) = [];
    %plot(1:size(A.c,1),tmp,'Color',[0 0 0 0.05])
            
    subplot(2,1,2);plot(1:size(B.c,1),B.probs(:,Binds),'LineWidth',2)
    xticks(floor(linspace(1,size(B.c,1),5)))
    %xticklabels(round(B.c(xticks,1),2))
    xticklabels("("+round(A.c(xticks,1),1)+","+round(B.c(xticks,1),1)+")")
    ylim([0 1.01]);xlabel('c_y');ylabel('probability');
    xlabel('(c_x,c_y)')
    title(YvarName); hold on;
    %plot(1:size(B.c,1),Yq(:,1),'--k','LineWidth',2)
    %plot(1:size(B.c,1),Yq(:,2),'k','LineWidth',2)
    %plot(1:size(B.c,1),Yq(:,3),'-.k','LineWidth',2)
    plot(1:size(B.c,1),Yqn(:,[1:4,6:9]),':k','LineWidth',1.5)
    plot(1:size(B.c,1),Yqn(:,5),'-.k','LineWidth',1.75)
    Bnames = replace(Y.names(Binds),'_',' ');
    %legend([Bnames;"Q_1";"Q_2";"Q_3"],'Location','northeastoutside', ...
    legend([Bnames;"Deciles"],'Location','northeastoutside', ...
        'AutoUpdate','off')
    tmp = B.probs;
    tmp(:,Binds) = [];
    %plot(1:size(B.c,1),tmp,'Color',[0 0 0 0.05])

set(gcf,'Position',[100,100,1200,1000])
fontsize(gcf,14,"points")
sgtitle(XvarName + ' - ' + YvarName)  
exportgraphics(gcf,"plots/"+Xvar+'_'+Yvar + "_stabPath_top"+string(nSel)+ ...
    ".png", 'Resolution',300)