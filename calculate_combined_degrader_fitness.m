function err = calculate_combined_degrader_fitness(k,model,options)
% Function for calculating model fitness to MRM timecourse data.


    % Setting parameters in model according to k.
    if strcmp(options.OptScale,'log')
        k = exp(k);
    end
    for ii = 1:length(options.ParamsToOpt)
        model.parameters.(options.ParamsToOpt(ii)) = k(ii);    
    end
    
    % Getting species indices that correspond to P, EP, EGP
    EP_ind = find(strcmp(model.species_names,"EP"),1);
    EGP_ind = find(strcmp(model.species_names,"EGP"),1);
    P_ind = find(strcmp(model.species_names,"P"),1);
    
    

    
    % Iterating through each of the conditions we want to optimize,
    % simulating model at parameters in k at the conditions in the data,
    % and calculating error.
    
    % Initalizing error
    err = 0;
    for ii = 1:length(options.ConditionsToOpt)
        err_tmp = 0;
        % Extracting data for the condition being simulated
        ConditionData = model.experimental_data{options.ConditionsToOpt(ii)};

        model.parameters.G_name = ConditionData{5};
        model.parameters.P_name = ConditionData{6};
        model.parameters.experiment_type = ConditionData{8};
        
        % Simulating model at timepoints specificed by t_exp and
        % calculating fitness
        if strcmp(ConditionData{7}, "timecourse")
            y_baseline = model.parameters.(append("Pt_",model.parameters.P_name,"_MRM"));
            t_exp = ConditionData{1};
            y_exp = ConditionData{2};
            y_err = ConditionData{3};
            model.parameters.(options.MG_Parameter) = ConditionData{4};
            ModelPredictions = model.evaluate(t_exp);
            if size(ModelPredictions,1) ~= length(t_exp)
                err_tmp = 1000000;
            else
                EP = ModelPredictions(:,EP_ind);
                EGP = ModelPredictions(:,EGP_ind);
                P = ModelPredictions(:,P_ind);
     
                y_predicted = (EP + EGP + P)./y_baseline;
        
        
                err_tmp = err_tmp + options.ErrorFunction(y_predicted,y_exp,y_err);
            end
        elseif strcmp(ConditionData{7},"titration")

            if strcmp(ConditionData{8},"MRM")
                y_baseline = model.parameters.(append("Pt_",model.parameters.P_name,"_MRM"));
                t_exp = ConditionData{1};
                y_exp = ConditionData{2};
                y_err = ConditionData{3}; 
                G_exp = ConditionData{4};
    
                for jj = 1:length(G_exp)
                    model.parameters.(options.MG_Parameter) = G_exp(jj);
                    ModelPredictions = model.evaluate([0 t_exp]);
                    EP = ModelPredictions(end,EP_ind);
                    EGP = ModelPredictions(end,EGP_ind);
                    P = ModelPredictions(end,P_ind);
                    y_predicted = (EP + EGP + P)/y_baseline;
                    err_tmp = err_tmp + options.ErrorFunction(y_predicted,y_exp(jj),y_err(jj));
                end
            elseif strcmp(ConditionData{8},"BRET")

                t_exp = ConditionData{1};
                y_exp = ConditionData{2};
                y_err = ConditionData{3}; 
                G_exp = ConditionData{4};
                
                % Simulating model in absence of any glue.
                model.parameters.(options.MG_Parameter) = 0;
                BaselinePredictions = model.evaluate([0 t_exp]);
                y_baseline = BaselinePredictions(end, EP_ind) + BaselinePredictions(end, EGP_ind);

                % Iterating through glue doses and calculating error
                for jj = 1:length(G_exp)
                    model.parameters.(options.MG_Parameter) = G_exp(jj);
                    ModelPredictions = model.evaluate([0 t_exp]);
                    EP = ModelPredictions(end,EP_ind);
                    EGP = ModelPredictions(end,EGP_ind);
                    y_predicted = (EP + EGP)/y_baseline;
                    err_tmp = err_tmp + options.ErrorFunction(y_predicted,y_exp(jj), y_err(jj));

                end


            elseif strcmp(ConditionData{8},"AS")
                t_exp = ConditionData{1};
                y_exp = ConditionData{2};
                y_err = ConditionData{3}; 
                G_exp = ConditionData{4};
                
                % Simulating model in absence of any glue.
                model.parameters.(options.MG_Parameter) = 0;
                BaselinePredictions = model.evaluate([0 t_exp]);
                y_baseline = BaselinePredictions(end, EP_ind) + BaselinePredictions(end, EGP_ind);

                % Iterating through glue doses and calculating error
                for jj = 1:length(G_exp)
                    model.parameters.(options.MG_Parameter) = G_exp(jj);
                    ModelPredictions = model.evaluate([0 t_exp]);
                    EP = ModelPredictions(end,EP_ind);
                    EGP = ModelPredictions(end,EGP_ind);
                    y_predicted = (EP + EGP)/y_baseline;
                    err_tmp = err_tmp + options.ErrorFunction(y_predicted,y_exp(jj), y_err(jj));
                end
            elseif strcmp(ConditionData{8},"HIBIT")
                y_baseline = model.parameters.(append("Pt_HIBIT"));
                t_exp = ConditionData{1};
                y_exp = ConditionData{2};
                y_err = ConditionData{3}; 
                G_exp = ConditionData{4};
    
                for jj = 1:length(G_exp)
                    model.parameters.(options.MG_Parameter) = G_exp(jj);
                    ModelPredictions = model.evaluate([0 t_exp]);
                    EP = ModelPredictions(end,EP_ind);
                    EGP = ModelPredictions(end,EGP_ind);
                    P = ModelPredictions(end,P_ind);
                    y_predicted = (EP + EGP + P)/y_baseline;
                    err_tmp = err_tmp + options.ErrorFunction(y_predicted,y_exp(jj),y_err(jj));
                end
            end
        end
       err = err + err_tmp*options.ConditionWeights(ii);
    end
end
