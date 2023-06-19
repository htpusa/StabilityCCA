function plotCCApath(A,c,vars)

% plotCCApath(A,c,vars) Plot the output of `SCCAec` for one view
%   plotCCApath(A,c,vars) plots the regularisation path calculated by
% `SCCAec` for one view
%
%   INPUTS:
%   A       -   p-by-l matrix where every column is a canonical coefficient
%               vector
%   c       -   l-by-1 matrix of regularisation parameter values
%   vars    -   vector of indices <=p, variables to be highlighted
%               Empty matrix: highlight nothing
%
%   EXAMPLE:
%   [A,B,~,~,~,cxy] = SCCAec(X,Y);
%   plotCCApath(A,cxy(:,1),1:5)


%figure
    hold on
    remain = 1:size(A,1);
    remain(vars) = [];
    plot(1:numel(c),A(remain,:),':k')
    xticks(floor(linspace(1,numel(c),5)))
    xticklabels(round(c(xticks),2))
    ylim([min(A,[],'all')-0.01, max(A,[],'all')+0.01])
    xlabel('c');ylabel('canonical coefficient');
    if ~isempty(vars)
        plot(1:numel(c),A(vars,:),'LineWidth',2)
    end
    hold off