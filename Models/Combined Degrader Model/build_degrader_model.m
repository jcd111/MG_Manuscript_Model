%% Function for building Combined MG Degrader Model

addpath(genpath('Modeling Functions'))
addpath(genpath('Models/Combined Degrader Model'))

% Loading parameter information
[parameters,bounds,parameter_names] = initialize_degrader_model_parameters;

% Leaving experimetnal data blank, setting when fitted
experimental_data = {};

% setting simulation functions
eval_function = @degrader_model_eval_function;
rate_laws = @degrader_model_rate_laws;

% Species names
species_names = ["P","E","G","EP","EG","EGP"];

% building and saving
model = modelstructure('parameters',parameters,'bounds',bounds,'parameter_names', ...
    parameter_names,'experimental_data',experimental_data,'eval_function', ...
    eval_function,'rate_laws',rate_laws,'species_names',species_names);

save('Models/Combined Degrader Model/degrader_model','model');

