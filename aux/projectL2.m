function xProj = projectL2(x,c)

% xProj = projectL2(x,c)
% L2-projection

xProj = c*x/norm(x,2);