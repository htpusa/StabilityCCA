mkdir('plots')

sp_emSel = 10;
en_emSel = 15;
sp_mSel = 20;
en_mSel = 20;

plotStabPath("species","idMetabolites","Species","Identified metabolites",sp_emSel)
plotStabPath("species","metabolites","Species","All metabolites",sp_mSel)
% plotStabPath("genus","idMetabolites","Genera","Identified metabolites")
% plotStabPath("genus","metabolites","Genera","All metabolites")
plotStabPath("enzymes","idMetabolites","Enzymes","Identified metabolites",en_emSel)
plotStabPath("enzymes","metabolites","Enzymes","All metabolites",en_mSel)

%plotRanking("species","idMetabolites","Species","Identified metabolites")
plotRanking("species","metabolites","Species","All metabolites")
% plotRanking("genus","idMetabolites","Genera","Identified metabolites")
% plotRanking("genus","metabolites","Genera","All metabolites")
%plotRanking("enzymes","idMetabolites","Enzymes","Identified metabolites")
plotRanking("enzymes","metabolites","Enzymes","All metabolites")

plotCCA("species","idMetabolites","Species","Identified metabolites",sp_emSel)
plotCCA("species","metabolites","Species","All metabolites",sp_mSel)
% plotCCA("genus","idMetabolites","Genera","Identified metabolites",10)
% plotCCA("genus","metabolites","Genera","All metabolites",15)
plotCCA("enzymes","idMetabolites","Enzymes","Identified metabolites",en_emSel)
plotCCA("enzymes","metabolites","Enzymes","All metabolites",en_mSel)