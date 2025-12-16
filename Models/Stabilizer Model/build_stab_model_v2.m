%% Function for building sRB model

addpath(genpath('Modeling Functions'))
addpath(genpath('Models/Stabilizer Model'))

% Loading parameter information
[parameters,bounds,parameter_names] = initialize_stabilizer_model_parameters_v2;

% Leaving experimetnal data blank, setting when fitted
experimental_data = {};

% setting simulation functions
eval_function = @stabilizer_model_eval_function;
rate_laws = @stabilizer_model_rate_laws;

% Species names
species_names = ["E","G","P","EG","EGP","EGP_cov"];

% building and saving
model = modelstructure('parameters',parameters,'bounds',bounds,'parameter_names', ...
    parameter_names,'experimental_data',experimental_data,'eval_function', ...
    eval_function,'rate_laws',rate_laws,'species_names',species_names);

save('Models/Stabilizer Model/Stabilizer_model_v4','model');

%% Testing simulation.

model.parameters.Gt = 50;


tspan = linspace(0,1000,1000);

y = model.evaluate(tspan);
P = sum(y(:,[6]),2);

figure(1)
clf
hold on
plot(tspan,P,'LineWidth',1.5)

xlabel("Time (hours)")
ylabel("Total P")
l = legend("Thalidomide","Pomalidomide","Lenalidomide","Location","SouthWest");
title(l,"Drug")



