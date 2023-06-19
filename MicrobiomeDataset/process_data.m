clear

%% metabolites
metaboliteData = readtable("data/metabolites.xlsx");
metaboliteHeader = readmatrix("data/metabolites.xlsx","Range",'1:8', ...
    "OutputType","string");

sampleID = metaboliteHeader(2,2:end)';
rawData = table2array(metaboliteData(8:end,2:end))';
metaboliteID = string(metaboliteData.x_Feature_Sample(8:end));
age = metaboliteHeader(3,2:end)';
diagnosis = metaboliteHeader(4,2:end)';

metaData = table(sampleID, age, diagnosis);
metaData.Properties.VariableNames = {'sampleID','age','status'};
metabolites.metaData = metaData;
metabolites.rawData = rawData;
metabolites.metaboliteID = metaboliteID;
metabolites.data = glog(metabolites.rawData,1e-8);

% metabolite metadata
metaboliteMetaData = readmatrix("data/metabolite_metadata.xlsx", "OutputType",...
    "string");
[~, I] = ismember(metaboliteID, metaboliteMetaData(:,1));
metabolites.putativeClass = metaboliteMetaData(I,5);
metabolites.exactMatch = metaboliteMetaData(I,6);
metabolites.best = metabolites.metaboliteID;
metabolites.best(~ismissing(metabolites.putativeClass)) = ...
    metabolites.putativeClass(~ismissing(metabolites.putativeClass));
metabolites.best(~ismissing(metabolites.exactMatch)) = ...
    metabolites.exactMatch(~ismissing(metabolites.exactMatch));

%% species
speciesData = readtable("data/species.xlsx");
speciesHeader = readmatrix("data/species.xlsx","Range",'1:8', ...
    "OutputType","string");

sampleID = speciesHeader(2,2:end)';
rawData = table2array(speciesData(8:end,2:end))';
speciesID = string(speciesData.x_Feature_Sample(8:end));
age = speciesHeader(4,2:end)';
diagnosis = speciesHeader(5,2:end)';

metaData = table(sampleID, age, diagnosis);
metaData.Properties.VariableNames = {'sampleID','age','status'};
species.metaData = metaData;
species.rawData = rawData;
species.data = glog(rawData,1e-8);
species.speciesID = speciesID;

% genus
genusID = string([]);
for i=1:numel(speciesID)
    ID = split(speciesID(i),'_');
    ID = ID(1);
    if ~ismember(ID, genusID)
        genusID = [genusID; ID];
    end
end

genusRawData = zeros(numel(sampleID), numel(genusID));
for i=1:numel(genusID)
    genusRawData(:,i) = sum(rawData(:,startsWith(speciesID,genusID(i))),2);
end

metaData = table(sampleID, age, diagnosis);
metaData.Properties.VariableNames = {'sampleID','age','status'};
genus.metaData = metaData;
genus.rawData = genusRawData;
genus.data = glog(genusRawData,1e-8);
genus.genusID = genusID;

%% enzymes
enzymeData = readtable("data/enzymes.xlsx");
enzymeHeader = readmatrix("data/enzymes.xlsx","Range",'1:8', ...
    "OutputType","string");

sampleID = enzymeHeader(2,2:end)';
rawData = table2array(enzymeData(8:end,2:end))';
enzymeID = string(enzymeData.x_Feature_Sample(8:end));
age = enzymeHeader(4,2:end)';
diagnosis = enzymeHeader(5,2:end)';

metaData = table(sampleID, age, diagnosis);
metaData.Properties.VariableNames = {'sampleID','age','status'};
enzymes.metaData = metaData;
enzymes.rawData = rawData;
enzymes.data = glog(rawData,1e-8);
enzymes.enzymeID = enzymeID;

%% are the samples in the same order? Yes.
all(species.metaData.sampleID==metabolites.metaData.sampleID)
all(enzymes.metaData.sampleID==metabolites.metaData.sampleID)

%%
save('data','metabolites','species','genus','enzymes');















