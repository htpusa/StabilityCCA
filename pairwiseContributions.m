function pc = pairwiseContributions(X,Y,a,b)

% pairwiseContributions Pairwise contributions to canonical correlation
%   pc = pairwiseContributions(X,Y,a,b) calculates the relative
%   contributions of each pair of variables to the CCA objective maximised
%   by coefficient vectors 'a' and 'b'
%
%   INPUTS:
%   X           -   n-by-px data matrix
%   Y           -   n-by-py data matrix
%   a           -   px-by-1 coefficient vector
%   b           -   py-by-1 coefficient vector
%
%   OUTPUTS:
%   pc          -   px-by-py vector of contributions where the entry ij
%                   is the proportion of the canonical correlation due to
%                   the correlation between the projections of X-variable i
%                   and Y-variable j
%
%   EXAMPLE:
%      [a,b] = SCCAec(X,Y,'cxy',[2 2]);
%       pc = pairwiseContributions(X,Y,a,b);

pc = ((a .* (X'*Y))' .* b)' / ((X*a)'*(Y*b));