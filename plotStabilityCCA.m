function plotStabilityCCA(probs,c,vars)

% plotStabilityCCA(probs,c,vars) Plot a CCA stability path
%   plotStabilityCCA(probs,c,vars) creates a plot of the output of
%   stabilityCCA for one view and one canonical variable.
%
%   INPUTS:
%   probs       -   l-by-p matrix with probabilities of variables to be 
%                   selected as a function of l regularisation parameter
%                   values (one page from stabilityCCA output structure)
%   c           -   l-by-1 matrix of regularisation parameter values
%   vars        -   vector of indices <=p, variables to be highlighted
%                   Empty matrix: highlight nothing
%
%   EXAMPLE:
%   [A,B] = stabilityCCA(X,Y);
%   plotStabilityCCA(A.probs,A.c,1:5)


%figure
    hold on
    remain = 1:size(probs,2);
    remain(vars) = [];
    plot(1:numel(c),probs(:,remain),':k')
    xticks(floor(linspace(1,numel(c),5)))
    xticklabels(round(c(xticks),2))
    ylim([0 1.01]);xlabel('c');ylabel('probability');
    if ~isempty(vars)
        plot(1:numel(c),probs(:,vars),'LineWidth',2)
    end
    hold off
    grid on