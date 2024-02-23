function y = softThresh(x,c)

% y = softThresh(x,c)
% Soft thresholding operator
%
% EXAMPLE
% x = linspace(-1,1,100)
% y = softThresh(x,0.4);
% plot(x,y)
% ylim([-1,1])

y = abs(x)-c;
y(y<0) = 0;
y = sign(x) .* y;