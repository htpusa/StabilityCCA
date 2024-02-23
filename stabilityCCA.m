function [A,B] = stabilityCCA(X,Y,varargin)

% stabilityCCA(X,Y) Canonical correlation stability path
%   [A,B] = stabilityCCA(X,Y) calculates CCA stability paths [1,2] for
%   the variables in the columns of X and Y.
%   Given a sequence of pairs of regularisation parameters [cx cy], X and Y
%   are concurrently split in two randomly 50 times, and a regularisation
%   path is calculated over [cx cy] using a sparse CCA function. The
%   probability of a variable to be selected at given [cx cy] is the
%   portion of subsets it has a non-zero (or above threshold) coefficient.
%
%   INPUTS:
%   X           -   n-by-px data matrix
%   Y           -   n-by-py data matrix
%   OPTIONAL INPUTS:
%   'SCCA'      -   which sparse CCA method to use:
%                       "SCCAec" (default)
%                       "PMDCCA"
%   'cxy'       -   l-by-2 matrix of regularisation parameters for X and Y
%                   views: cxy(:,1) for X, and cxy(:,2) for Y.
%                   Should be in ascending order (decreasing sparsity)
%                   default: default values used by SCCAec
%   'D'         -   how many canonical vectors are found (default: 1)
%   'eps'       -   threshold for selecting variables: variable i is
%                   considered selected in coefficient vector x if
%                   |x_i| > eps * max(|x|) (default: 0)
%   'verbose'   -   Boolean, print stuff or not (default: 1)
%
%   OUTPUTS:
%   A           -   a structure containing the results for the X-view:
%                   .probs  -   l-by-px-by-D matrix with probabilities of X
%                       variables to be selected as a function of l 
%                       regularisation parameter values
%                   .auc    -   px-by-D vector, average selection
%                       probability for each variable
%                   .numSel -   l-by-D vector, average number of variables
%                       selected at each regularisation parameter value
%                   .c     -   l-by-1 vector of regularisation parameter 
%                       values used
%   B           -   a structure containing the results for the Y-view
%
%   EXAMPLE:
%      load carbig;
%      data = [Displacement Horsepower Weight Acceleration MPG];
%      nans = sum(isnan(data),2) > 0;
%      X = data(~nans,1:3); Y = data(~nans,4:5);
%      [A,B] = stabilityCCA(X,Y);

%   References:
%       [1] Meinshausen, Nicolai, and Peter Bühlmann. "Stability 
%           selection." Journal of the Royal Statistical Society: Series B 
%           (Statistical Methodology) 72.4 (2010): 417-473.
%       [2] Shah, Rajen D., and Richard J. Samworth. "Variable selection 
%           with error control: another look at stability selection." 
%           Journal of the Royal Statistical Society: Series B (Statistical
%           Methodology) 75.1 (2013): 55-80.

verbose = 1;
eps = 0;
D = 1;
cxy = [];
fun = "SCCAec";

if ~isempty(varargin)
    if rem(size(varargin, 2), 2) ~= 0
		error('Check optional inputs.');
    else
        for i = 1:2:size(varargin, 2)
            switch varargin{1, i}
                case 'SCCA'
					fun = varargin{1, i+1};
                case 'D'
					D = varargin{1, i+1};
                case 'cxy'
					cxy = varargin{1, i+1};
                case 'eps'
					eps = varargin{1, i+1};
                case 'verbose'
					verbose = varargin{1, i+1};
                otherwise
					error(['Could not recognise optional input names.' ...
                        '\nNo input named "%s"'],...
						varargin{1,i});
            end
        end
    end
end

if ~ismember(fun,["SCCAec";"PMDCCA"])
    error('SCCA should be SCCAec or PMDCCA')
end

rounds = 100;
parts = false(size(X,1),rounds);
for i=1:2:rounds
    part = crossvalind('HoldOut',size(X,1));
    parts(:,i:i+1) = [part, ~part];
end
nParts = size(parts,2);

if isempty(cxy)
    [~,~,~,~,~,cxy] = SCCAec(X,Y);
end
l = size(cxy,1);

probsA = zeros(l,size(X,2),D);
probsB = zeros(l,size(Y,2),D);
numSelA = zeros(l,nParts,D);
numSelB = zeros(l,nParts,D);

parfor i=1:nParts
    verbose && fprintf('\t Subsample %d of %d\n',i,nParts);
    Xtmp = X(parts(:,i),:);
    Ytmp = Y(parts(:,i),:);
    if fun=="SCCAec"
        [Atmp,Btmp] = SCCAec(Xtmp,Ytmp,'cxy',cxy,'D',D);
    else
        [Atmp,Btmp] = PMDCCA(Xtmp,Ytmp,'cxy',cxy,'D',D);
    end
    Asel = abs(Atmp)>eps*max(abs(Atmp),[],1);
    probsA = probsA + pagetranspose(Asel);
    numSelA(:,i,:) = sum(Asel,1);
    Bsel = abs(Btmp)>eps*max(abs(Btmp),[],1);
    probsB = probsB + pagetranspose(Bsel);
    numSelB(:,i,:) = sum(Bsel,1);
end

probsA = probsA/nParts;
AUCA = squeeze(pagetranspose(sum(probsA,1)/rounds));
numSelA = squeeze(mean(numSelA,2));
A.probs = probsA;
A.auc = AUCA;
A.numSel = numSelA;

A.c = cxy(:,1);
probsB = probsB/nParts;
AUCB = squeeze(pagetranspose(sum(probsB,1)/rounds));
numSelB = squeeze(mean(numSelB,2));
B.probs = probsB;
B.auc = AUCB;
B.numSel = numSelB;
B.c = cxy(:,2);