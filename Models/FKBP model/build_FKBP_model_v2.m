%% Initializing Parameters
clear;
addpath(genpath('Models/FKBP model'))
addpath(genpath('Modeling Functions'))

[parameters, bounds, parameter_names] = initialize_FKBP_parameters_v2();

%% Setting  eval function

eval_function = @(obj,tspan) FKBP_eval_function_v2(obj);

%% Setting up experimental data

experimental_data = {};

%% Constructing model.

model = modelstructure('parameters',parameters,'parameter_names',parameter_names,...
    'bounds',bounds,'experimental_data',experimental_data,'species_names',...
    ["EGP","EG","G","E","P"], 'eval_function',eval_function);

save('Models/FKBP model/FKBPmodel_v2','model')

%%
load('FKBPmodel_v1.mat')
model.parameters.G_name = "Cmpd7";
model.parameters.P_name = "FRB";

model.parameters.Pt=1e-4;
BaselinePredictions = model.evaluate;
EGP_Baseline = BaselinePredictions(1);
y_baseline = EGP_Baseline;
P_range = logspace(-4,4);
y_predicted = zeros(length(P_range),1);

for ii = 1:length(P_range)
    model.parameters.Pt = P_range(ii);
    tmp1 = model.evaluate;
    y_predicted(ii) = ( tmp1(1))/y_baseline;
    
end

figure(1)
clf
hold on
plot(P_range,y_predicted,'b','LineWidth',1.5)
set(gca,'Xscale','log')
