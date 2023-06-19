function xProj = projectL1L2(x,c)

% L1L2-projection from [1].

%   [1] Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. 
%    "A penalized matrix decomposition, with applications to sparse 
%    principal components and canonical correlation analysis." 
%    Biostatistics 10.3 (2009): 515-534.

convergenceCrit = 1e-5;
maxIter = 100;

if c<=1
    [m i] = max(abs(x));
    xProj = zeros(numel(x),1);
    xProj(i) = c*sign(x(i));
else
    xProj = projectL2(x,1);
end

cont = norm(xProj,1) > c;
deltaRange = [0 max(abs(x))];
iter = 0;

while cont
    delta = mean(deltaRange);
    xProj = projectL2(softThresh(x,delta),1);
    diff = norm(xProj,1) - c;
    if diff>0
        deltaRange(1) = delta;
    else
        deltaRange(2) = delta;
    end
    cont = abs(diff) > convergenceCrit && iter<=maxIter;
    iter = iter+1;
end
