function [A B r U V cxy] = SCCAec(X,Y,varargin)

% SCCAec Sparse canonical correlation analysis
%   [A B r U V] = SCCAec(X,Y) performs sparse canonical correlation
%   analysis using the alternating projected gradient algorithm from [1] 
%   and the L1-L2-penalty from [2].
%   Finds canonical coefficient vectors a such that |a|_1 <= c and
%   |a|_2 <= 1 simultaneously.
%   A regularisation path is calculated for a sequence of cx,cy pairs using
%   a warmstart algorithm that proceeds from less to more regularisation.
%   Multiple coefficient vectors can be found by using deflation from [3].
%
%   INPUTS:
%   X           -   n-by-px data matrix
%   Y           -   n-by-py data matrix
%   OPTIONAL INPUTS:
%   
%   'cxy'       -   l-by-2 matrix of regularisation parameters for X and Y
%                   views: cxy(:,1) for X, and cxy(:,2) for Y.
%                   Should be in ascending order (decreasing sparsity)
%                   default: log space with 100 steps from 1 (fully sparse)
%                   to sqrt(p) (no sparsity)
%   'D'         -   how many canonical vectors are found (default: 1)
%   'maxIter'   -   maximum number of iterations
%   'eps'       -   tolerance parameter
%
%   OUTPUTS:
%   A           -   px-by-l-by-D matrix with canonical coefficients for X
%                   in columns
%   B           -   py-by-l-by-D matrix with canonical coefficients for Y
%                   in columns
%   r           -   l-by-D vector with the sample canonical correlations
%   U           -   n-by-l-by-D matrix with canonical variables/scores for 
%                   X in columns
%   V           -   n-by-l-by-D matrix with canonical variables/scores for 
%                   Y in columns
%   cxy         -   l-by-2 matrix of regularisation parameters
%
%   EXAMPLE:
%      load carbig;
%      data = [Displacement Horsepower Weight Acceleration MPG];
%      nans = sum(isnan(data),2) > 0;
%      X = data(~nans,1:3); Y = data(~nans,4:5);
%      [A,B] = SCCAec(X,Y);

%   References:
%     [1] Uurtio, Viivi, Sahely Bhadra, and Juho Rousu. "Large-scale sparse 
%       kernel canonical correlation analysis." International Conference on
%       Machine Learning. PMLR, 2019.
%     [2] Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. 
%       "A penalized matrix decomposition, with applications to sparse 
%       principal components and canonical correlation analysis." 
%       Biostatistics 10.3 (2009): 515-534.
%     [3] Mackey, Lester. "Deflation methods for sparse PCA." Advances in 
%       neural information processing systems 21 (2008).

px = size(X,2);
py = size(Y,2);
nSteps = 100;
cxy = [logspace(log10(1),log10(sqrt(px)),nSteps)',...
        logspace(log10(1),log10(sqrt(py)),nSteps)'];
D = 1;
param.maxIter = 500;
param.eps = 1e-10;

if ~isempty(varargin)
    if rem(size(varargin, 2), 2) ~= 0
		error('Check optional inputs.');
    else
        for i = 1:2:size(varargin, 2)
            switch varargin{1, i}
                case 'D'
					D = varargin{1, i+1};
                case 'cxy'
                    cxy = varargin{1, i+1};
                case 'maxIter'
					param.maxIter = varargin{1, i+1};
                case 'eps'
					param.eps = varargin{1, i+1};
                otherwise
					error(['Could not recognise optional input names.' ...
                        '\nNo input named "%s"'],...
						varargin{1,i});
            end
        end
    end
end

if size(cxy,2) ~= 2
    error('cxy should have 2 columns');
elseif ~( isequal(cxy(:,1),sort(cxy(:,1))) && ...
        isequal(cxy(:,2),sort(cxy(:,2))) )
    warning('cxy not properly sorted, strange things may happen')
end

l = size(cxy,1);
nx = size(X,1);
ny = size(Y,1);

A = zeros(px,l,D);
B = zeros(py,l,D);
r = zeros(l,D);
U = zeros(nx,l,D);
V = zeros(ny,l,D);

X = X - mean(X,1);
Y = Y - mean(Y,1);

for d=1:D
    [aInit,~,bInit] = svd(X'*Y);
    aInit = aInit(:,1);
    bInit = bInit(:,1);
    [A(:,end,d),B(:,end,d)] = SCCAecFromInit(X,Y,...
                                cxy(end,1),cxy(end,2),...
                                aInit,bInit,...
                                param);
    U(:,end,d) = X*A(:,end,d);
    V(:,end,d) = Y*B(:,end,d);
    r(end,d) = SCCAecObjective(U(:,end,d),V(:,end,d));
    for i=l-1:-1:1
        [A(:,i,d),B(:,i,d)] = SCCAecFromInit(X,Y,...
                                cxy(i,1),cxy(i,2),...
                                A(:,i+1,d),B(:,i+1,d),...
                                param);
        U(:,i,d) = X*A(:,i,d);
        V(:,i,d) = Y*B(:,i,d);
        r(i,d) = SCCAecObjective(U(:,i,d),V(:,i,d));
    end
    % deflate data
    normA = A(:,end,d)/norm(A(:,end,d),2);
    normB = B(:,end,d)/norm(B(:,end,d),2);
    X = X - (normA*normA'*X')';
    Y = Y - (normB*normB'*Y')';
end



