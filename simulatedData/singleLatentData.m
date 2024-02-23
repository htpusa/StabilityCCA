function [X,Y,Xt,Yt] = singleLatentData(n,a,b,noise)

% generate single latent factor data with ground truth vector a and b
    
mu = randn(n+100,1);

X = normrnd(mu(1:n)*a',noise);
Y = normrnd(mu(1:n)*b',noise);
Xt = normrnd(mu(n+1:end)*a',noise);
Yt = normrnd(mu(n+1:end)*b',noise);