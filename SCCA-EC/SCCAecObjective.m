function obj = SCCAecObjective(Xa,Yb)

% obj = SCCAecObjective(Xa,Yb)
% CCA objective ie correlation between projections X*a and Y*b

obj = Xa'*Yb/(sqrt(Xa'*Xa)*sqrt(Yb'*Yb));