function model = probPathwayConstruction( targetMet, model )
%probPathwayConstruction( targetMet, pathway )
%   Constructs the additional pathway of producing target metabolite
%   using probabilistic selection scheme
    orgModel = model;
    model = constructPath(targetMet, model, 0);
    sol = optimizeCbModel(model);
    native = 0; % native does not produce targetMet
    if sol.f < native
        model = orgModel;
    end
end

function model = constructPath( met, model, chainLength )
    limit = 6;
    % return unmodified model if exceeded allowed chain length
    if chainLength > limit
        return
    end

    % Grab all reactions from KEGG that has met
    reactions = KEGG(met);
    % Foreach reaction, if exist ignore, else add to set of reactions
    picks = [];
    for i = 1:length(reactions)
        % if reaction i is not in model
        if reactions(i) == model
            % add reaction i to model
            picks = [picks reactions(i)];
        end
    end
    
    %Add the randomly selected reaction to pathway
    rSelected = picks(randi(picks));
    model = addReaction(rSelected);
    
    % Get all reactants of rSelected
    reactants = [];
    for m = 1:length(reactants)
        if strmatch(reactants(m), model.mets) ~= cellstr()
            continue
        else
            constructPath(m,model,chainLength+1)
        end
    end
end

