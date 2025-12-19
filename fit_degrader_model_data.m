function outstruct = fit_dose_response_data(model,options)
% Function for fitting dose response data of either BRET or MM1S data.
% Declaring arguments and setting default values of parameters

arguments
    % model
    model modelstructure
    % optimization method
    options.OptMethod string = "particleswarm"
    % Optimization Functions
    options.OptOptions = optimoptions("particleswarm")
    % Scale of Optimization
    options.OptScale = "log";
    % Fitness function
    options.FitFunction function_handle = @calculate_dose_response_fitness
    % Error function.
    options.ErrorFunction function_handle = @se
    % Timepoints to simulate
    options.tspan = [0 5];
    % Parameter defining glue concentration
    options.MG_Parameter string = "Gt"
    % Parameters to optimize6 7;
    options.ParamsToOpt string = {"KD1", "KD2"}
    % Conditions we are optimizing to
    options.ConditionsToOpt double = [];
    % Weights for each condition being optimized
    options.ConditionWeights double = [];
    % If using patternsearch, number of initial sets generated
    options.NSets double = 100;
    % If using patternsearch, best N sets to take and optimize
    options.BestNSets double = 10;
    % If using patternsearch, can provide starting points here.
    options.StartingPoints double = [];
    % Use parallel computing for patternsearch?
    options.UseParallel logical = false
    % If using parallel, number of pools to use
    options.NPools double = [];
    % savename
    options.SaveName string = "";
end

% If no parameters are given, selecting all for optimization
if isempty(options.ParamsToOpt)
    options.ParamsToOpt = model.parameter_names;
end

% If no conditions are given, selecting all
if isempty(options.ConditionsToOpt)
    options.ConditionsToOpt = 1:length(model.experimental_data);
end

% If no condition weights given, setting all to 1.
if isempty(options.ConditionWeights)
    options.ConditionWeights = ones(size(options.ConditionsToOpt));
end
% If number of pools is not specified, one for each set optimized.
if isempty(options.NPools)
    options.NPools = options.BestNSets;
end
% Optimizing data based on fitting method.
if strcmp(options.OptMethod,"particleswarm")

    FitFunction = @(k) options.FitFunction(k,model,options);
    bounds = reformat_boundaries(model.bounds,options.ParamsToOpt);
    if strcmp(options.OptScale,"log")
        ub = log(bounds(:,2)); lb = log(bounds(:,1));
    else
        ub = bounds(:,2);lb = bounds(:,1);
    end
    tic
    [k_opt,err] = particleswarm(FitFunction,length(options.ParamsToOpt),lb,ub,options.OptOptions);
    RunTime = toc;
    % Setting parameters in model to optimized values
    for ii = 1:length(k_opt)
        if strcmp(options.OptScale,"log")
            model.parameters.(options.ParamsToOpt(ii)) = exp(k_opt(ii));
        else
            model.parameters.(options.ParamsToOpt(ii)) = k_opt(ii);
        end

    end


elseif strcmp(options.OptMethod,"patternsearch")

    % Setting up Fitness Function
    FitFunction = @(k) options.FitFunction(k,model,options);
    % Generating ensemble of points to be optimized.
    bounds = reformat_boundaries(model.bounds,options.ParamsToOpt);
    if strcmp(options.OptScale,"log")
        bounds = log(bounds);
    end
    InitialSwarm = create_initial_swarm(options.StartingPoints,bounds,options.ParamsToOpt,options.NSets);
    % Calculating Fitness for Initial Points.
    ErrInitial = zeros(options.NSets,1);
    for ii = 1:options.NSets
        ErrInitial(ii) = FitFunction(InitialSwarm(:,ii));
    end
    % Sorting based on fitness.
    [~,ind] = sort(ErrInitial);
    InitialSwarm = InitialSwarm(:,ind);
    % iterating through each parameter set and locally optimizing using
    % patternsearch.


    err = zeros(options.BestNSets,1);
    OptimizedEnsemble = zeros(size(InitialSwarm));
    rng('shuffle')
    tic
    if options.UseParallel
        parfor (ii = 1:options.BestNSets,options.NPools)
            ParameterSet = InitialSwarm(:,ii);
            [tmp1,tmp2] = patternsearch(FitFunction,ParameterSet,[],[],[],[],bounds(:,1),bounds(:,2),options.OptOptions);
            OptimizedEnsemble(:,ii) = tmp1;
            err(ii) = tmp2;
            fprintf("Finished set %i of %i.\n",ii,options.BestNSets)
        end
    else
        for ii = 1:options.BestNSets

            ParameterSet = InitialSwarm(:,ii);
            [tmp1,tmp2] = patternsearch(FitFunction,ParameterSet,[],[],[],[],bounds(:,1),bounds(:,2),options.OptOptions);
            OptimizedEnsemble(:,ii) = tmp1;
            err(ii) = tmp2;
            fprintf("Finished set %i of %i.\n",ii,options.BestNSets)
        end
    end
    RunTime = toc;
    [err,ind] = sort(err);
    OptimizedEnsemble = OptimizedEnsemble(:,ind);
    if strcmp(options.OptScale,'log')
        OptimizedEnsemble = exp(OptimizedEnsemble);
    end


elseif strcmp(options.OptMethod,"particleswarm_ensemble")
    % Setting up Fitness Function
    FitFunction = @(k) options.FitFunction(k,model,options);
    bounds = reformat_boundaries(model.bounds,options.ParamsToOpt);
    if strcmp(options.OptScale,"log")
        ub = log(bounds(:,2)); lb = log(bounds(:,1));
    else
        ub = bounds(:,2);lb = bounds(:,1);
    end
    err = zeros(options.NSets,1);
    OptimizedEnsemble = zeros(length(options.ParamsToOpt),length(options.NSets));
    rng('shuffle')
    tic
    if options.UseParallel
        ParamsToOpt = options.ParamsToOpt;
        OptOptions = options.OptOptions;
        parfor (ii = 1:options.NSets,options.NPools)
            [tmp1,tmp2] = particleswarm(FitFunction,length(ParamsToOpt),lb,ub,OptOptions);
            OptimizedEnsemble(:,ii) = tmp1;
            err(ii) = tmp2;
            fprintf("Finished set %i of %i.\n",ii,options.NSets)
        end
    else
        for ii = 1:options.NSets
            [tmp1,tmp2] = particleswarm(FitFunction,length(options.ParamsToOpt),lb,ub,options.OptOptions);
            OptimizedEnsemble(:,ii) = tmp1;
            err(ii) = tmp2;
            fprintf("Finished set %i of %i.\n",ii,options.NSets)
        end
    end
    RunTime = toc;
    [err,ind] = sort(err);
    OptimizedEnsemble = OptimizedEnsemble(:,ind);
    if strcmp(options.OptScale,'log')
        OptimizedEnsemble = exp(OptimizedEnsemble);
    end



end


% Storing results in outstruct.
outstruct = struct();

% storing optimized model
outstruct.model = model;
% storing error
outstruct.fitness = err;
% storing runntime
outstruct.RunTime = RunTime;
outstruct.fitted_parameters = options.ParamsToOpt;
% If an ensemble was generated, storing ensemble here.
if exist("OptimizedEnsemble","var")
    outstruct.OptimizedEnsemble = OptimizedEnsemble;
end
if ~strcmp(options.SaveName,"")
    save(options.SaveName,'outstruct');
end
end

%% Helper Functions
function bounds = reformat_boundaries(old_bounds,params_to_opt)
for ii = 1:length(params_to_opt)
    tmp1 = append(params_to_opt{ii},"_bounds");
    tmp2 = append("bounds(",num2str(ii),",:",") = old_bounds.('",tmp1,"');");
    eval(tmp2);
end
end

function initial_swarm = create_initial_swarm(StartingPoints,bounds,params_to_opt,n_sets)

if isempty(StartingPoints)
    initial_swarm = bounds(:,1) + rand(length(params_to_opt),n_sets).*(bounds(:,2)...
        -bounds(:,1));
elseif size(StartingPoints,2) >= n_sets
    initial_swarm = StartingPoints(:,1:n_sets);
else
    tmp1 = bounds(:,1) + rand(length(params_to_opt),n_sets-size(StartingPoints,2).*(bounds(:,2)...
        -bounds(:,1)));
    initial_swarm = [StartingPoints, tmp1];
end


end

