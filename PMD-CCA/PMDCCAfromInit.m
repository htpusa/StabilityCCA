function [a,b] = PMDCCAfromInit(X,Y,cx,cy,aInit,bInit,param)

% [a,b] = PMDCCAfromInit(X,Y,cx,cy,aInit,bInit,param)
% Performs the alternating covariance maximisation algorithm starting from
% an initial point.

maxIter = param.maxIter;
eps = param.eps;

a = projectL2(aInit,1);
b = projectL2(bInit,1);

obj = 0;
iter = 0;
improvement = 42;

while improvement>eps && iter<maxIter
    
    a = projectL1L2(X'*Y*b,cx);
    b = projectL1L2(Y'*X*a,cy);
    objNew = a'*X'*Y*b;
    improvement = (objNew-obj)/abs(obj);
    obj = objNew;
    iter = iter + 1;
end

if iter==maxIter
    warning('PMDCCAfromInit reached maximum number of iterations')
end