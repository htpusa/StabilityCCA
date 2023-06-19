function gx = glog(x,lambda)

% Generalised logarithm transformation
%
% gx = glog(x,lambda)
%
% EXAMPLE
% x = rand(5,1);
% gx = glog(x,1);

    gx = log(x + (x.^2 + lambda).^0.5);