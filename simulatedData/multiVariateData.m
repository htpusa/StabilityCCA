function [X,Y,Xt,Yt] = multiVariateData(n,a,b)

% Generate multivariate Gaussian data with ground truth vectors a and b

    alpha = 0.5;
    
    sigmaX = zeros(numel(a));
    sigmaX(find(a),find(a)) = alpha;
    sigmaX = sigmaX - diag(sigmaX).*eye(numel(a)) + eye(numel(a));

    sigmaY = zeros(numel(b));
    sigmaY(find(b),find(b)) = alpha;
    sigmaY = sigmaY - diag(sigmaY).*eye(numel(b)) + eye(numel(b));

    sigmaXY = a*b';
    sigmaYX = sigmaXY';
    sigma = [sigmaX sigmaXY;sigmaYX sigmaY];
    sigma = makePosDef(sigma);
    XY = mvnrnd(zeros(numel(a)+numel(b),1),sigma,n+100);
    X = XY(1:n,1:numel(a));
    Y = XY(1:n,numel(a)+1:end);
    Xt = XY(n+1:end,1:numel(a));
    Yt = XY(n+1:end,numel(a)+1:end);

    function Apd = makePosDef(A)
        if all(eig(A)>=0)
            Apd = A;
        else
            d = size(A,1);
            [V,D] = eig(A);
            D = diag(D);
            [D,ord] = sort(D,'descend');
            V = V(:,ord);
            tol = d*max(abs(D))*eps;
            delta = 2*tol;
            tau = max(0, delta-D);
            dA = V*diag(tau)*V';
            Apd = A+dA;
        end
    end

end