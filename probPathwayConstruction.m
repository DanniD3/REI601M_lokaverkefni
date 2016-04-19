function model = probPathwayConstruction( targetMet, model, KEGGDB )
%probPathwayConstruction( targetMet, pathway )
%   Constructs the additional pathway of producing target metabolite
%   using probabilistic selection scheme
    
    model = addReaction(model,strcat('EX_',targetMet),targetMet);
    model = changeObjective(model,strcat('EX_',targetMet));
    orgModel = model; % store original model for comparison
    tmpModel = model; % store the model with the best current pathway
    
    native = 0; % check if targetMet is native
    if ~isempty(strmatch(targetMet,model.mets))
        sol = optimizeCbModel(model);
        native = sol.f;
    end
    
    % chain limit can be set in the first line in function contructPath
    maxTimeInSec = 10;
    
    tic;
    while toc < maxTimeInSec
        model = constructPath(targetMet, model, KEGGDB, 0);
        sol = optimizeCbModel(model);
        if sol.f <= native
            model = tmpModel; % return to current best if current is bad
        else
            tmpModel = model; % store current as current best if is better
        end
    end
    disp(strcat('Native :', num2str(native), ' Pathway :', num2str(sol.f)));
    if sol.f <= native
        model = orgModel;
        disp('Original model is better');
    end
end

function model = constructPath( met, model, KEGGDB, chainLength )
    limit = 6;
    disp(strcat(num2str(chainLength),'. reaction adding'));
    % return unmodified model if exceeded allowed chain length
    if chainLength > limit
        disp('Reaction Limit REACHED!!');
        return;
    end

    % Grab all reactions from KEGG that has met
    reactions = searchRxnWithKeggID(met,KEGGDB);
    
    % Size difference for 1 rxn or multiple rxns
    rSize = size(reactions);
    disp(strcat('Adding reaction for :', met));
    disp(strcat('Found :', num2str(rSize(1)), ' reactions for :', met));
    disp('');
    
    if rSize(1) > 1
        % Case multiple rxns from KEGG
        
        % Add the randomly selected reaction to pathway
        rAdded = 0;
        while ~rAdded
            rLength = size(reactions);
            % Check if reactions exist
            if isempty(reactions)
                throw(MException('NoReactionAvailable','No reaction has been found for the metabolite'));
            end
            i = randi(rLength(1));
            reaction = reactions(i,:);
            rxnName = reaction(2);
            rxnFormula = reaction(3);
            [model,rxnIDexists] = addReaction(model, rxnName{1}, rxnFormula{1});
            if isempty(rxnIDexists)
                % Added a rection if rxnID did not exist
                rAdded = 1;
            else
                % Remove existing reaction from the list of selectable
                % reaction
                reactions(i,:) = [];
            end
        end
        disp(rxnFormula);
        
        % Get all reactants of rSelected
        reactants = getReactants(reaction(3), met);
        for m = 1:length(reactants)
            if checkProduction(model,reactants(m))
                continue
            else
                r = reactants(m);
                constructPath(r{1},model,KEGGDB,chainLength+1);
            end
        end
    else
        % Case at most 1 rxn from KEGG
        if isempty(reactions)
            throw(MException('NoReactionAvailable','No reaction has been found for the metabolite'));
        end
        
        % add the reaction to model
        rxnName = reactions(2);
        rxnFormula = reactions(3);
        [model,rxnIDexists] = addReaction(model, rxnName{1}, rxnFormula{1});
        %disp(rxnFormula);
        
        % Get all reactants of rSelected
        reactants = getReactants(reactions(3), met);
        disp(reactants);
        for m = 1:length(reactants)
            if checkProduction(model,reactants(m))
                continue
            else
                r = reactants(m);
                constructPath(r{1},model,KEGGDB,chainLength+1);
            end
        end
    end
end

function rxns = searchRxnWithKeggID( mets, KEGGDB )
    rxns = [];
    for i = 1:length(KEGGDB(:,3))
        metName = KEGGDB(i,3);
        if ~isempty(strfind(metName{1},mets))
            rxns = [rxns;KEGGDB(i,:)];
        end
    end
%     rxns = [rxns;KEGGDB(231,:)]; % SUCC GOLDEN APPLE
end

function reactants = getReactants( reactionFormula, met )
    C = strsplit(reactionFormula{1},{' ','+'},'CollapseDelimiters',true);
    reactants = [];
    products = [];
    side = 0;
    targetSide = 0;
    for i = 1:length(C)
        % switch side when <=>
        if strcmp(C(i),'<=>')
           side = 1;
        else
            % add metabolite to correct side
            if ~side
                reactants = [reactants C(i)];
            else
                products = [products C(i)];
            end
            % set targetSide when target is found
            if strcmp(C(i),met)
                targetSide = side;
            end
        end
    end
    
    % return products if target is on RHS
    if ~targetSide
        reactants = products;
    end
    
    % check reactants are not multipliers
    for i = length(reactants):1:-1
        r = reactants(i);
        rStr = r{1};
        if ~strcmp(rStr(1),'C')
            reactants(i) = [];
        end
    end
end

function isProduced = checkProduction( model, reactant )
    isProduced = 0;
    productions = full(model.S(strmatch(reactant,model.mets),:));
    for i = 1:length(productions)
        if productions(i) > 0
            isProduced = 1;
            break;
        end
    end
end

