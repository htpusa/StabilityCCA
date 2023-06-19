% Create simulated example data

clear

% data 1
n = 100;

a1 = [ones(5,1); -ones(5,1); zeros(90,1)];
a2 = [-ones(2,1); ones(2,1); -ones(2,1); ones(4,1); zeros(90,1)];
b1 = [zeros(90,1); -ones(5,1); ones(5,1)];
b2 = [zeros(90,1); ones(4,1); -ones(2,1); ones(2,1); -ones(2,1);];

Z = randn(n,2); W = orth(Z);
w1 = Z(:,1); w2 = Z(:,2);

X = normrnd(w1*a1' + w2*a2', 0.09); X = normalize(X);
Y = normrnd(w1*b1' + w2*b2', 0.09); Y = normalize(Y);

data1.X = X;
data1.Y = Y;
data1.A = [a1 a2];
data1.B = [b1 b2];

% data 2
n = 50;

a1 = [ones(5,1); -ones(5,1); zeros(490,1)];
b1 = [zeros(980,1); -ones(10,1); ones(10,1)];

w1 = randn(n,1);

X = normrnd(w1*a1', 0.25); X = normalize(X);
Y = normrnd(w1*b1', 0.25); Y = normalize(Y);

data2.X = X;
data2.Y = Y;
data2.A = a1;
data2.B = b1;

save('stabCCA_example','data1','data2');

