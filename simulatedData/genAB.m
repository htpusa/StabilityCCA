function [a,b] = genAB(conf)

% generate a ground truth
% conf should be [px_true py_true px_noise py_noise]

aT = rand(conf(1),1);
bT = rand(conf(2),1);

a = zeros(conf(1)+conf(3),1);
tmp = randperm(numel(a));
a(tmp(1:numel(aT))) = aT;

b = zeros(conf(2)+conf(4),1);
tmp = randperm(numel(b));
b(tmp(1:numel(bT))) = bT;