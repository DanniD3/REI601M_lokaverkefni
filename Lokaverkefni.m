addpath('c:/cobra')
initCobraToolbox
%% Extract model
model = readCbModel('yeast_7.00_cobra.xml');
writeCbModel(model,'text','yeast_7.00.txt');
save('yeast_7.00.dat','-struct','model');
%% Load model
model = load('yeast_7.00.dat','-mat');
model.S = full(model.S);
%% UTIL: Get mets
metAndNames = [model.mets model.metNames];
disp(metAndNames);
%% UTIL: Get rxns
rxnAndNames = [model.rxns model.rxnNames];
disp(rxnAndNames);
%% Mets we need
h2o = model.mets{strmatch('H2O [cytoplasm]',model.metNames)};
o2 = model.mets{strmatch('oxygen [cytoplasm]',model.metNames)};
oleate = model.mets{strmatch('oleate [cytoplasm]',model.metNames)};
malonyl_coa = model.mets{strmatch('malonyl-CoA [cytoplasm]',model.metNames)};
co2 = model.mets{strmatch('carbon dioxide [cytoplasm]',model.metNames)};
coa = model.mets{strmatch('coenzyme A [cytoplasm]',model.metNames)};

%% Adding pathway to model for EPA
% oleic acid -> linoleic acid
model = addMetabolite(model,'LA','linoleic acid');
model = addReaction(model,'d12d',{oleate,o2,'LA',h2o}, [-1 -1 1 2],false);

% linoleic acid -> alpha-linoleic acid
model = addMetabolite(model,'aLA','alpha-linoleic acid');
model = addReaction(model,'d15d',{'LA',o2,'aLA',h2o}, [-1 -1 1 2],false);

% alpha-linoleic acid -> stearidonic acid
model = addMetabolite(model,'STA','stearidonic acid');
model = addReaction(model,'d6d',{'aLA',o2,'STA',h2o}, [-1 -1 1 2],false);

% stearidonic acid -> eicosatetraenoic acid
model = addMetabolite(model,'ETA','eicosatetraenoic acid');
model = addReaction(model,'c18e',{'STA',malonyl_coa,'ETA',co2,coa}, [-1 -1 1 1 1],false);

% eicosatetraenoic acid -> ecosapentaenoic acid
model = addMetabolite(model,'EPA','eicosapentaenoic acid');
model = addReaction(model,'d5d',{'ETA',o2,'EPA',h2o}, [-1 -1 1 2],false);

%% Add EX for linoleic acid
model = addReaction(model,'EX_la','LA');
model = changeObjective(model,'EX_la');
sol = optimizeCbModel(model,'max');
disp(sol);
%% Add EX for eicosapentaenoic acid
model = addReaction(model,'EX_epa','EPA');
model = changeObjective(model,'EX_epa');
sol = optimizeCbModel(model,'max');
disp(sol);
%% Robustness Analysis
bioRxn = 'r_2111';
model = changeObjective(model,bioRxn);
robustnessAnalysis(model,'EX_epa');
%% Additional pathway for DHA from EPA
% eicosapentaenoic acid -> docosapentaenoic acid
model = addMetabolite(model,'DPA','docosapentaenoic acid');
model = addReaction(model,'c20e',{'EPA',malonyl_coa,'DPA',co2,coa}, [-1 -1 1 1 1],false);

% docosapentaenoic acid -> docosahexaenoic acid
model = addMetabolite(model,'DHA','docosahexaenoic acid');
model = addReaction(model,'d4d',{'DPA',o2,'DHA',h2o}, [-1 -1 1 2],false);

%% Add EX for docosahexaenoic acid
model = addReaction(model,'EX_dha','DHA');
model = changeObjective(model,'EX_dha');
sol = optimizeCbModel(model,'max');

%% Robustness Analysis
bioRxn = 'r_2111';
model = changeObjective(model,bioRxn);
robustnessAnalysis(model,'EX_dha');

%% KEGG DB
excelImp = importdata('unidb.xlsx');
nums = num2cell(excelImp.data.CuratedDB(:,1));
dats = excelImp.textdata.CuratedDB(2:end,2:4);

KEGGDB = [nums dats];
