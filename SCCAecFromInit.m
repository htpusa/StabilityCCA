function [a,b] = SCCAecFromInit(X,Y,cx,cy,aInit,bInit,param)

% [a,b] = SCCAecFromInit(X,Y,cx,cy,aInit,bInit,param)
% Performs the alternating projected gradient algorithm for sparse 
% canonical correlation starting from an initial point.

maxIter = param.maxIter;
eps = param.eps;

a = projectL2(aInit,1);
b = projectL2(bInit,1);
a = projectL1L2(a,cx);
b = projectL1L2(b,cy);
Xa = X*a;
Yb = Y*b;

obj = SCCAecObjective(Xa,Yb);
objOld = obj;
iter = 0;
improvement = 42;

while improvement>eps && iter<maxIter

    % line search for a
    gradientA = SCCAecGrad(X,Xa,Yb);
    contA = 1;
    %gamma = norm(gradientA,2);
    gamma = norm(a,2);
    while contA && gamma>1e-10
        aNew = projectL1L2(a+gamma*gradientA,cx);
        XaNew = X*aNew;
        objNew = SCCAecObjective(XaNew,Yb);
        if objNew > obj + 1e-4*abs(obj)
            contA = 0;
            a = aNew;
            Xa = XaNew;
            obj = objNew;
        else
            gamma = gamma/2;
        end
    end
    % line search for b
    gradientB = SCCAecGrad(Y,Yb,Xa);
    contB = 1;
    %gamma = norm(gradientB,2);
    gamma = norm(b,2);
    while contB && gamma>1e-10
        bNew = projectL1L2(b+gamma*gradientB,cy);
        YbNew = Y*bNew;
        objNew = SCCAecObjective(Xa,YbNew);
        if objNew > obj + 1e-4*abs(obj)
            contB = 0;
            b = bNew;
            Yb = YbNew;
            obj = objNew;
        else
            gamma = gamma/2;
        end
    end

    improvement = (obj-objOld)/abs(obj+objOld);
    objOld = obj;
    iter = iter + 1;
end

if iter==maxIter
    warning('SCCAecFromInit reached maximum number of iterations')
end