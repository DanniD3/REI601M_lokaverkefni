function model = probPathwayConstruction( targetMet, model, KEGGDB )
%probPathwayConstruction( targetMet, pathway )
%   Constructs the additional pathway of producing target metabolite
%   using probabilistic selection scheme
    orgModel = model;
    model = constructPath(targetMet, model, KEGGDB, 0);
    model = addReaction(model,strcat('EX_',targetMet),targetMet);
    model = changeObjective(model,strcat('EX_',targetMet));
    sol = optimizeCbModel(model);
    disp(sol);
    native = 0; % native does not produce targetMet
    if sol.f < native
        model = orgModel;
    end
end

function model = constructPath( met, model, KEGGDB, chainLength )
    limit = 6;
    % return unmodified model if exceeded allowed chain length
    if chainLength > limit
        return
    end

    % Grab all reactions from KEGG that has met
    reactions = searchRxnWithKeggID(met,KEGGDB);
    % Foreach reaction, if exist ignore, else add to set of reactions
    picks = [];
    
    % Size difference for 1 rxn or multiple rxns
    rSize = size(reactions);
    if rSize(1) > 1
        % Case multiple rxns from KEGG
        %   This is incompletely since only 1 rxn was found
        for i = 1:length(reactions)
            % if reaction i is not in model
            reaction = reactions(i);
            rxn = reaction(2);
            match = strmatch(rxn{1},model.rxns);
            if ~isempty(match)
                % add reaction i to model
                picks = [picks reactions(i)];
            end
        end
        %Add the randomly selected reaction to pathway
        rSelected = picks(randi(length(picks)));
        model = addReaction(model, rSelected(2), rSelected(3));
        
        % Get all reactants of rSelected
        reactants = [];
        for m = 1:length(reactants)
            if strmatch(reactants(m), model.mets) ~= cellstr()
                continue
            else
                constructPath(m,model,KEGGDB,chainLength+1)
            end
        end
    else
        % Case at most 1 rxn from KEGG
        if isempty(reactions)
            return
        end
        
        % add the reaction to model
        rxnName = reactions(2);
        rxnFormula = reactions(3);
        [model,rxnIDexists] = addReaction(model, rxnName{1}, rxnFormula{1});
        disp(rxnFormula);
        % Get all reactants of rSelected
        reactants = getReactants(reactions(3), met);
        for m = 1:length(reactants)
            % NEED TO CHECK IF REACTANTS ARE PRODUCED IN S
            if ~isempty(strmatch(reactants(m), model.mets))
                continue
            else
                constructPath(m,model,KEGGDB,chainLength+1)
            end
        end
    end
end

function rxns = searchRxnWithKeggID( mets, KEGGDB )
    rxns = [];
    for i = 1:length(KEGGDB(:,3))
        metName = KEGGDB(i,3);
        if ~isempty(strfind(metName{1},mets))
            rxns = [rxns KEGGDB(i,:)];
        end
    end
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
end
