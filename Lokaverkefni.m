addpath('c:/cobra');
initCobraToolbox;
%% Extract model (of yeast_7.00.dat is not present)
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
disp(sol);

%% Robustness Analysis
bioRxn = 'r_2111'; % growth
model = changeObjective(model,bioRxn);
robustnessAnalysis(model,'EX_dha');

%% Adding KEGGID to yeast model
addpath('c:/jsonparser') % http://www.mathworks.com/matlabcentral/fileexchange/42236-parse-json-text/content/example/html/usage.html
fname = 'yeast_to_kegg.json';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

keggData = JSON.parse(str);
fields = fieldnames(keggData);
for i = 1:length(fields)
    field = fields(i);
    model.mets(strmatch(field{1},model.mets)) = cellstr(keggData.(field{1}));
    model.rxns(strmatch(field{1},model.rxns)) = cellstr(keggData.(field{1}));
end

%% UTIL: Search for met in KEGGDB
oleateRxn = {};
for i = 1:length(KEGGDB(:,4))
    metName = KEGGDB(i,4);
    if ~isempty(strfind(metName{1},'Oleate'))
        oleateRxn = [oleateRxn;KEGGDB(i,:)];
    end
end
disp(oleateRxn);

%% KEGG DB needs work
excelImp = importdata('unidb.xlsx');
nums = num2cell(excelImp.data.CuratedDB(:,1));
dats = excelImp.textdata.CuratedDB(2:end,2:4);

KEGGDB = [nums dats];
%% Using Ecoli for probPathwayConstruction
model = readCbModel('ecoli_core_model.xml');
%keggList = [model.mets model.metKEGGID];

[EcNum EcText EcRaw] = xlsread('ecoli_core_model.xls','metabolites');
KEGGIDs = EcText(2:end,8);

%This list is correct, 99% sure (checked first and last)
mets = [model.mets KEGGIDs];

%adding KeggID to model
for i = 1:length(KEGGIDs)
    model.mets{i} = KEGGIDs{i};
end

%% Running probPathwayConstruction
%KEGG IDs:
oleateKEGG = 'C00712';
succKEGG = 'C00042';
hKEGG = 'C00001';
%dhaKEGG = 'C06429';
model = probPathwayConstruction(succKEGG, model, KEGGDB);

%% Finding oleate producing reactions
oleateS = model.S(strmatch(oleate, model.mets),:);
oleateRxnIs = [];
for i = 1:length(oleateS)
    if oleateS(i) > 0
        oleateRxnIs = [oleateRxnIs i];
    end
end
oleateRxns = model.rxns(oleateRxnIs);
disp(oleateRxns);

%% Find bottleneck formula
%bn = 'r_1277'; % bn f. dha
bn = 'r_2189'; % bn f. oleate
bnFormula = model.S(:,strmatch(bn,model.rxns));
bnMets = [];
for i = 1:length(bnFormula)
    if bnFormula(i) ~= 0
       bnMets = [bnMets; [model.mets(i) bnFormula(i) model.metNames(i)]] ;
    end
end
disp(bnMets);

%% BOTTLENECK CHECKS

% Add EX for oleic acid
model = addReaction(model,'EX_oa',oleate);
model = changeObjective(model,'EX_oa');
sol = optimizeCbModel(model,'max');
disp(sol);

%load Ecoli_core_model
targetRxn='EX_oa';

npoints = 1000;  % Stær?slembiúrtaks (fjöldi flæðisvigra, 5000 er betra)
nsec = 120;      % 1800 fyrir minni genome-scale líkön, 4800 fyrir stærri líkön
[s,mixedFraction] = gpSampler(model, npoints, [], nsec);

% mixedFraction gefur til kynna hversu góð nálgun fékkst, gildi?er ?bilinu 0 til 1.
% Æskilegt er a?f?gildi nálægt 0.5, ef þa?er t.d. 0.75 eða hærra þarf a?hækka 'nsec'.
fprintf('Mixed fraction=%1.2f\n', mixedFraction)


% Skoða dreifingu flæðisgilda fyrir PGI hvarfi??glýkólýsu
% Sj?líka plotSampleHist ?COBRA
nbins=25;
hist(s.points(findRxnIDs(model,'d4d'),:), nbins) % Sjáum a?hvarfi?gengur ?báðar áttir, meiri líkur ?a?þa?gangi til hægri

% Ákvarða fylgnistuðla
R=corrcoef(s.points');
nbest=20;
idxTarget=findRxnIDs(model,targetRxn);
[~,idx]=sort(abs(R(idxTarget,:)),'descend'); % Ignore directionality, focus on magnitude
for i=1:nbest
    fprintf('%s\t%1.2f\t(%s)\n', model.rxns{idx(i)}, R(idxTarget,idx(i)), model.subSystems{idx(i)})
end
