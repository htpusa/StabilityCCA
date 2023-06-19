function plotRanking(Xvar,Yvar,XvarName,YvarName)

path = "res/" + Xvar + '_' + Yvar;
load(path + "/ranking.mat");

% selection stats
px = size(scca.A,1);py = size(scca.B,1);
maxRank = size(sccaRank.score,2);
nRounds = size(scca.score,1);

scca.card = sum(scca.A~=0,1) + sum(scca.B~=0,1);

sccaRank.sel = zeros(px+py,maxRank);
sccaRank.prop = zeros(nRounds,maxRank);
stabCCA.stabCCA.sel = zeros(px+py,maxRank);
stabCCA.stabCCA.prop = zeros(nRounds,maxRank);
stabCCA.stabCCA2.sel = zeros(px+py,maxRank);
stabCCA.stabCCA2.prop = zeros(nRounds,maxRank);
stabCCA.stabCCA3.sel = zeros(px+py,maxRank);
stabCCA.stabCCA3.prop = zeros(nRounds,maxRank);
stabCCA.stabCCA4.sel = zeros(px+py,maxRank);
stabCCA.stabCCA4.prop = zeros(nRounds,maxRank);
stabCCA.stabCCA5.sel = zeros(px+py,maxRank);
stabCCA.stabCCA5.prop = zeros(nRounds,maxRank);
for r=1:nRounds
    for nSel=1:maxRank
        sccaRank.sel(sccaRank.ranking(1:nSel,r),nSel) = ...
            sccaRank.sel(sccaRank.ranking(1:nSel,r),nSel) + 1;
        sccaRank.prop(r,nSel) = sum(sccaRank.ranking(1:nSel,r)<=px)/nSel;
        stabCCA.stabCCA.sel(stabCCA.stabCCA.ranking(1:nSel,r),nSel) = ...
            stabCCA.stabCCA.sel(stabCCA.stabCCA.ranking(1:nSel,r),nSel) + 1;
        stabCCA.stabCCA.prop(r,nSel) = ...
            sum(stabCCA.stabCCA.ranking(1:nSel,r)<=px)/nSel;
        % 2
        stabCCA.stabCCA2.sel(stabCCA.stabCCA2.ranking(1:nSel,r),nSel) = ...
            stabCCA.stabCCA2.sel(stabCCA.stabCCA2.ranking(1:nSel,r),nSel) + 1;
        stabCCA.stabCCA2.prop(r,nSel) = ...
            sum(stabCCA.stabCCA2.ranking(1:nSel,r)<=px)/nSel;
        % 3
        stabCCA.stabCCA3.sel(stabCCA.stabCCA3.ranking(1:nSel,r),nSel) = ...
            stabCCA.stabCCA3.sel(stabCCA.stabCCA3.ranking(1:nSel,r),nSel) + 1;
        stabCCA.stabCCA3.prop(r,nSel) = ...
            sum(stabCCA.stabCCA3.ranking(1:nSel,r)<=px)/nSel;
        % 4
        stabCCA.stabCCA4.sel(stabCCA.stabCCA4.ranking(1:nSel,r),nSel) = ...
            stabCCA.stabCCA4.sel(stabCCA.stabCCA4.ranking(1:nSel,r),nSel) + 1;
        stabCCA.stabCCA4.prop(r,nSel) = ...
            sum(stabCCA.stabCCA4.ranking(1:nSel,r)<=px)/nSel;
        % 5
        stabCCA.stabCCA5.sel(stabCCA.stabCCA5.ranking(1:nSel,r),nSel) = ...
            stabCCA.stabCCA5.sel(stabCCA.stabCCA5.ranking(1:nSel,r),nSel) + 1;
        stabCCA.stabCCA5.prop(r,nSel) = ...
            sum(stabCCA.stabCCA5.ranking(1:nSel,r)<=px)/nSel;
    end
end
sccaRank.stab = nogueiraStability(sccaRank.sel,nRounds);
stabCCA.stabCCA.stab = nogueiraStability(stabCCA.stabCCA.sel,nRounds);
stabCCA.stabCCA2.stab = nogueiraStability(stabCCA.stabCCA2.sel,nRounds);
stabCCA.stabCCA3.stab = nogueiraStability(stabCCA.stabCCA3.sel,nRounds);
stabCCA.stabCCA4.stab = nogueiraStability(stabCCA.stabCCA4.sel,nRounds);
stabCCA.stabCCA5.stab = nogueiraStability(stabCCA.stabCCA5.sel,nRounds);

figure
% correlation
    subplot(1,3,1)
    hold on
    % scca ranking
    start = find(all(~isnan(sccaRank.score),1),1);
    plot(start:maxRank,mean(sccaRank.score(:,start:maxRank),1),'-or',...
        'LineWidth',2)
    errorbar(start:maxRank,mean(sccaRank.score(:,start:maxRank),1),...
        std(sccaRank.score(:,start:maxRank)),'r','LineStyle','none')
    ylim([0 1.01]);ylabel('correlation');yticks(0:0.2:1)
    xlim([0,maxRank+1]);xticks([0:5:maxRank]);xlabel('number of variables')
    % auc
    start = find(all(~isnan(stabCCA.stabCCA.score),1),1);
    plot(start:maxRank,mean(stabCCA.stabCCA.score(:,start:maxRank),1),'-db',...
        'LineWidth',2)
    errorbar(start:maxRank,mean(stabCCA.stabCCA.score(:,start:maxRank),1),...
        std(stabCCA.stabCCA.score(:,start:maxRank)),'b','LineStyle','none')
    % scca full model
    plot([1,maxRank],repmat(mean(scca.score),2,1),'--k', 'LineWidth',2)
    errorbar([2,maxRank-1],repmat(mean(scca.score),2,1),...
        repmat(std(scca.score),2,1),'k','LineStyle','none')
    legend({'SCCA-EC ranking','','StabilityCCA ranking','', ...
        sprintf('Full SCCA-EC model (%.1f variables)',mean(scca.card))},...
        'Location','southeast')
    hold off;title('Correlation');grid on
% stability
    subplot(1,3,2)
    hold on
    % scca ranking
    start = 1;
    plot(start:maxRank,sccaRank.stab(start:maxRank),'-or',...
        'LineWidth',2)
    minim = 1-(px+py)/(px+py-maxRank);
    ylim([minim 1]);ylabel('Stability estimator');yticks(0:0.2:1)
    xlim([0,maxRank+1]);xticks([0:5:maxRank]);xlabel('number of variables')
    % auc
    start = 1;
    plot(start:maxRank,stabCCA.stabCCA.stab(start:maxRank),'-db',...
        'LineWidth',2)
    %legend({'gradSCCA ranking','Stable CCA ranking'},...
    %    'Location','northeast')
    hold off;title('Selection stability');grid on
% proportion
    subplot(1,3,3)
    hold on
    % scca ranking
    start = 1;
    plot(start:maxRank,mean(sccaRank.prop(:,start:maxRank),1),'-or',...
        'LineWidth',2)
    errorbar(start:maxRank,mean(sccaRank.prop(:,start:maxRank),1),...
        std(sccaRank.prop(:,start:maxRank)),'r','LineStyle','none')
    ylim([0 1.01]);ylabel('correlation');yticks(0:0.2:1)
    xlim([0,maxRank+1]);xticks([0:5:maxRank]);xlabel('number of variables')
    % auc
    start = 1;
    plot(start:maxRank,mean(stabCCA.stabCCA.prop(:,start:maxRank),1),'-db',...
        'LineWidth',2)
    errorbar(start:maxRank,mean(stabCCA.stabCCA.prop(:,start:maxRank),1),...
        std(stabCCA.stabCCA.prop(:,start:maxRank)),'b','LineStyle','none')
    grid on;title('Proportion of X-variables')

   sgtitle(XvarName + ' - ' + YvarName)

   set(gcf,'Position',[100 100 1300 500])
   fontsize(gcf,12,"points")

   exportgraphics(gcf,"plots/" + Xvar + "_" + Yvar + "_ranking.png", ...
       'Resolution',300)

%% different rankings
% figure
% % correlation
%     subplot(1,2,1)
%     hold on
%     % 1
%     start = find(all(~isnan(stabCCA.stabCCA.score),1),1);
%     plot(start:maxRank,mean(stabCCA.stabCCA.score(:,start:maxRank),1),'-db',...
%         'LineWidth',2)
%     %errorbar(start:maxRank,mean(sccaRank.score(:,start:maxRank),1),...
%     %    std(sccaRank.score(:,start:maxRank)),'r','LineStyle','none')
%     ylim([0 1.01]);ylabel('correlation');yticks(0:0.2:1)
%     xlim([0,maxRank+1]);xticks([0:5:maxRank]);xlabel('number of variables')
% 
%     % 2
%     start = find(all(~isnan(stabCCA.stabCCA2.score),1),1);
%     plot(start:maxRank,mean(stabCCA.stabCCA2.score(:,start:maxRank),1),'-or',...
%         'LineWidth',2)
% 
%     % 3
%     start = find(all(~isnan(stabCCA.stabCCA3.score),1),1);
%     plot(start:maxRank,mean(stabCCA.stabCCA3.score(:,start:maxRank),1),'-sg',...
%         'LineWidth',2)
% 
%     % 4
%     start = find(all(~isnan(stabCCA.stabCCA4.score),1),1);
%     plot(start:maxRank,mean(stabCCA.stabCCA4.score(:,start:maxRank),1),'-^k',...
%         'LineWidth',2)
% 
%     % 5
%     start = find(all(~isnan(stabCCA.stabCCA5.score),1),1);
%     plot(start:maxRank,mean(stabCCA.stabCCA5.score(:,start:maxRank),1),'-+m',...
%         'LineWidth',2)
% 
%     legend({'AUC','Odds ratio','AUC, alternate','Odds ratio, alternate', ...
%         'First to 0.8'},'Location','southeast')
% 
%     hold off;title('Correlation');grid on
% % Stability
%     subplot(1,2,2)
%     hold on
%     % 1
%     start = 1;
%     plot(start:maxRank,stabCCA.stabCCA.stab(start:maxRank),'-db',...
%         'LineWidth',2)
%     minim = 1-(px+py)/(px+py-maxRank);
%     ylim([minim 1]);ylabel('Stability estimator');yticks(0:0.2:1)
%     xlim([0,maxRank+1]);xticks([0:5:maxRank]);xlabel('number of variables')
% 
%     % 2
%     plot(start:maxRank,stabCCA.stabCCA2.stab(start:maxRank),'-or',...
%         'LineWidth',2)
% 
%     % 3
%     plot(start:maxRank,stabCCA.stabCCA3.stab(start:maxRank),'-sg',...
%         'LineWidth',2)
% 
%     % 4
%     plot(start:maxRank,stabCCA.stabCCA4.stab(start:maxRank),'-^k',...
%         'LineWidth',2)
% 
%     % 5
%     plot(start:maxRank,stabCCA.stabCCA5.stab(start:maxRank),'-+m',...
%         'LineWidth',2)
% 
%     hold off;title('Stability');grid on;
%     sgtitle('Different ranking strategies')
% 
%     set(gcf,'Position',[100 100 1200 600])
% 
%     exportgraphics(gcf,"plots/" + Xvar + "_" + Yvar + "_ranking_comp.png", ...
%        'Resolution',300)






