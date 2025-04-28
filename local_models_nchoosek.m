function [A, modelType, tot_n_models] = local_models_nchoosek(gridnames,gridstructure)

    % Initialize
    n_gridnames = length(gridnames);
    A = {}; modelType = {};
    count = 0;

    % How many models in total?
    tot_n_models = 0;
    for k = 1:n_gridnames
        sz = size(nchoosek(gridnames,k));
        tot_n_models = tot_n_models + sz(1);
    end

    % Create all possible grid expansions
    for k = n_gridnames:-1:1 % Decrementing k's
        
        chosen = nchoosek(gridnames,k);
        sz_chosen = size(chosen);
        n_chosen = sz_chosen(1);
        
        boolType = zeros(n_chosen,n_gridnames);
        for m = 1:n_chosen
            boolType(m,:) = ismember(gridnames,chosen(m,:));
        end

        % Populate final model structures
        for m = 1:n_chosen
            count = count+1;
            modelType{count} = boolType(m,:);

            hold = [];
            for n = 1:length(modelType{count})
                if modelType{m}(n) == 1
                    hold = horzcat(hold, getfield(gridstructure,chosen(n)));
                end
            end
            A{count} = hold;
        end
    end
    
    % Finished!
    A=A';
    modelType = modelType';

return