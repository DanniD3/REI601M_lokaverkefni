function index = findBioIndex( model )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    c = model.c;
    
    for i = 1:length(c)
        if c(i) == 1
            index = i;
        end
    end
end

