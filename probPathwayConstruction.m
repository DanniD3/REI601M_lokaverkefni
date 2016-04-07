function pathway = probPathwayConstruction( targetMet, selectionScheme, pathway )
%probPathwayConstruction( targetMet, selectionScheme )
%   Constructs the additional pathway of producing target metabolite
%   using probabilistic selection scheme

    pathway = constructPath(targetMet, selectionScheme, pathway);
    
    % RUN Flux Balance Analysis on pathway (sol.x)
    sol = optimizeCbModel(model);
    
    if sol.f < sol(pathway).f
        pathway = [];
    end
end

function pathway = constructPath( met,  selectionScheme, pathway )
    if length(pathway) < limit
        pathway = [];
        return
    end
    
    % KEGG database
    for i = 1:KEGG(met)
        if i == pathway
            rWeighting = 0;
        else
            %rWeighting = (based on) selectionScheme;
        end
    end
    
    %Randomly select a reaction based on rWeighting
    %Add the selected reaction to pathway
    %for %m = 1:SelectedReaction
        %if m == pathway
            %continue
        %else
            constructPath(m,selectionScheme,pathway)
        %end
    %end
end

